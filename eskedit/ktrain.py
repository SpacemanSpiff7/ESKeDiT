"""
The purpose of this module is to take an input VCF (and reference FASTA) to create mutation probabilities
"""
import sys
import datetime
import time
import multiprocessing as mp
import array
from collections import Counter, defaultdict
from signal import signal, SIGINT
from cyvcf2 import VCF
from pyfaidx import Fasta, FetchError
import pandas as pd
from eskedit import get_bed_regions
from eskedit.ktools import is_dash, is_quality_snv, kmer_search, get_methylation_probability, make_directory, is_cpg, \
    read_and_build_models
import os
from itertools import repeat


def ktrain_region_driver(bedregions: iter, vcfpath: str, fastapath: str, kmer_size: int, methylation_vcf_path: str):
    # may want to add this as a keyword argument later
    AC_cutoff = 1
    fasta = Fasta(fastapath, sequence_always_upper=True)
    gnomadVCF = VCF(vcfpath)
    reference_kmer_counts = Counter()
    rare_transitions_count = defaultdict(lambda: array.array('L', [0, 0, 0, 0]))
    # rare_transitions_ac = defaultdict(lambda: array.array('L', [0, 0, 0, 0]))
    # common_transitions_ac = defaultdict(lambda: array.array('L', [0, 0, 0, 0]))
    # common_transitions_count = defaultdict(lambda: array.array('L', [0, 0, 0, 0]))
    num_singletons, num_all_variants = 0, 0
    if methylation_vcf_path is not None:
        meth_vcf = VCF(methylation_vcf_path)
        lo_count, mid_count, hi_count = 0, 0, 0
        low_meth = defaultdict(lambda: array.array('L', [0, 0, 0, 0]))  # < 0.2
        mid_meth = defaultdict(lambda: array.array('L', [0, 0, 0, 0]))  # 0.2-0.6
        hi_meth = defaultdict(lambda: array.array('L', [0, 0, 0, 0]))  # # > 0.6

    nuc_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for region in bedregions:
        # bed indices will be same as start/stop stored by GRegion
        # Since gnomad variants are aligned to forawrd strand, do not consider strandedness
        try:
            seq = region.get_seq_from_fasta(fasta, kmer_size=kmer_size)
        except KeyError:
            print('Fasta record {} : {} - {} not found'.format(region.chrom, region.start, region.stop),
                  file=sys.stderr, flush=True)
            continue

        kmer_results = kmer_search(seq, kmer_size)
        # Add count of kmers to master dictionary
        reference_kmer_counts += kmer_results['count_kmers']
        for variant in gnomadVCF(region.gnomad_rep()):
            varAC = variant.INFO.get('AC')
            if is_quality_snv(variant):
                num_all_variants += 1
                if varAC <= AC_cutoff:
                    num_singletons += 1
                    # welcome to indexing hell
                    # VCF is 1-based[,], BED is 0-based [,)
                    # zero based idx = (VCFidx - 1) - (BED start
                    seq_idx = variant.POS - region.start + kmer_size // 2 - 1
                    seq_context = seq[seq_idx - kmer_size // 2: seq_idx + kmer_size // 2 + 1]
                    if 'N' in seq_context or len(seq_context) == 0:
                        continue

                    if len(seq_context) != kmer_size:
                        raise IndexError(
                            f'Simone, you idiot... kmer_size: {kmer_size} recovered seq len: {len(seq_context)}')
                    # check index
                    if variant.REF != seq_context[len(seq_context) // 2]:
                        print(
                            f'ERROR: Fasta REF {seq_context[len(seq_context) // 2]} and VCF REF {variant.REF} don\'t match at position {variant.POS} on {variant.CHROM}',
                            flush=True, file=sys.stderr)

                    # This is for counting allele information
                    # check methylation if C followed by G

                    if methylation_vcf_path is not None and is_cpg(seq_context):
                        # Despite the IDE's best efforts, 'meth_vcf' et al. are guaranteed to be defined here
                        methylation = get_methylation_probability(
                            meth_vcf(f'{variant.CHROM}:{variant.POS}-{variant.POS}'))
                        if methylation < 0.2:  # low/none
                            lo_count += 1
                            low_meth[seq_context][nuc_idx[variant.ALT[0]]] += 1
                        elif 0.2 <= methylation <= 0.6:
                            mid_count += 1
                            mid_meth[seq_context][nuc_idx[variant.ALT[0]]] += 1
                        elif 0.6 < methylation <= 1:
                            hi_count += 1
                            hi_meth[seq_context][nuc_idx[variant.ALT[0]]] += 1
                        else:
                            print('Irregular meth prob')

                    rare_transitions_count[seq_context][nuc_idx[variant.ALT[0]]] += 1

        # print summary of kmer search on GRegion fields to stdout
        region.add_field('num_variants', f'AC_lte_{AC_cutoff}={num_singletons},all_snvs={num_all_variants}')
        if methylation_vcf_path is not None:
            region.add_field('methylation', f'low={lo_count},intermediate={mid_count},high={hi_count}')
        print(str(region), flush=True)

    # Package return values

    def freeze_transitions(t_dict):
        newdict = defaultdict(tuple)
        for key, value in t_dict.items():
            newdict[key] = tuple(value)
        return frozenset(newdict.items())

    # return_dict = []
    # return_dict.append(frozenset({'ref_kcounts': frozenset(reference_kmer_counts.items())}.items()))
    # return_dict.append(frozenset({'rare_transitions': freeze_transitions(rare_transitions_count)}.items()))
    if methylation_vcf_path is not None:
        return frozenset(reference_kmer_counts.items()), \
               freeze_transitions(rare_transitions_count), \
               freeze_transitions(low_meth), \
               freeze_transitions(mid_meth), \
               freeze_transitions(hi_meth)
    else:
        return frozenset(reference_kmer_counts.items()), freeze_transitions(rare_transitions_count)


def ktrain(bedpath: str, vcfpath: str, fastapath: str, kmer_size: int, meth_vcf_path=None, header: bool = False,
           nprocs: int = 6):
    start_time = time.time()
    dirname = "{}_{}".format('.'.join(bedpath.split('.')[:-1]),
                             datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))
    try:
        nprocs = int(nprocs)
    except ValueError:
        print('ERROR: nprocs must be an integer', file=sys.stderr, flush=True)
        exit(1)

    # define functions for handling early termination
    def singint_handler(signal_received, frame):
        print('Kill signal received: {}\nFrame: {}'.format(signal_received, frame))
        exit(1)

    signal(SIGINT, singint_handler)

    region_chunks = get_bed_regions(bedpath, nchunks=nprocs, hasheader=header)
    driver_args = [(chunk, vcfpath, fastapath, kmer_size, meth_vcf_path) for chunk in region_chunks]

    pool = mp.Pool(nprocs)

    # https://stackoverflow.com/questions/56075837/parallel-execution-of-a-list-of-functions

    with mp.Pool(nprocs) as pool:
        futures = [pool.starmap_async(ktrain_region_driver, driver_args)]
        results = [fut.get() for fut in futures]
    # results = pool.starmap(ktrain_region_driver, driver_args)
    # pool.close()
    # pool.join()

    cumulative_results = []
    for ridx, result in enumerate(results[0]):
        for idx, table in enumerate(result):
            if ridx == 0:
                if idx == 0:
                    cumulative_results.append(Counter(dict(table)))
                else:
                    cumulative_results.append(
                        defaultdict(lambda: array.array('L', [0, 0, 0, 0]),
                                    {a: array.array('L', b) for a, b in dict(table).items()}))
                continue
            else:
                if idx == 0:
                    cumulative_results[0] += Counter(dict(table))
                else:
                    # idx_nuc = dict(zip(list('ACGT'), range(4)))
                    # cumulative_results[idx] = sum_transitions(cumulative_results[idx],
                    #                                           defaultdict(lambda: array.array('L', [0, 0, 0, 0]),
                    #                                                       {a: array.array('L', b) for a, b in
                    #                                                        dict(table).items()}))

                    for kmer, count_array in defaultdict(lambda: array.array('L', [0, 0, 0, 0]),
                                                         {a: array.array('L', b) for a, b in
                                                          dict(table).items()}).items():
                        for alt_idx, count in enumerate(count_array):
                            cumulative_results[idx][kmer][alt_idx] += count

    results_index_order = ['ref_kmer_counts', 'rare_transitions', 'low_methylation',
                           'intermediate_methylation', 'high_methylation']
    dirpath = make_directory(dirname=dirname)
    for i, res in enumerate(cumulative_results):
        df = pd.DataFrame.from_dict(res, orient='index').sort_index()
        if i > 0:
            df.columns = list('ACGT')
        df.to_csv(os.path.join(dirpath, '.'.join([results_index_order[i], 'csv'])))

    models = read_and_build_models(dirpath)
    model_dir_name = "{}_MultinomialModel_{}".format('.'.join(bedpath.split('.')[:-1]),
                                                     datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))
    model_dir_path = make_directory(dirname=model_dir_name)
    for name, model_df in models.items():
        model_df.to_csv(os.path.join(model_dir_path, '.'.join([f'{name}_probability', 'csv'])))
    print(f'Done in {time.time() - start_time}')

    # Fasta indexing is inclusive start and stop (idx starts at 1)
    # VCF indexing works same as Fasta
    # BED indexing starts at 0 and does not include end coordinate
    # for region in regions:
    #     # shift start left by half ksize to capture nucleotide level mutability
    #     fa_idx_start = max(region.start - kmer_size // 2, 1)
    #     fa_idx_stop = region.stop + kmer_size//2 + 1
    #     # bed indices will be same as start/stop stored by GRegion
    #     if region.fields['strand'] is None:
    #         seq = fasta[region.chrom][region.start:region.stop].seq
    #     else:
    #         if is_dash(region.fields['strand']):
    #             seq = fasta[region.chrom][region.start:region.stop].complement.seq
    #         else:
    #             seq = fasta[region.chrom][region.start:region.stop].seq
    #
    #     # Count k-mers
    #
    #     # Local kmer count
    #
    #     # global kmer count
