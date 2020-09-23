# Steps:
# 1. load model - multinom or freq?
#     a. use raw count data
#     b. option for probability
# 2. sliding window
# 3. input bed, vcf, fasta
# for region in bed
#     count kmers
#     output TFIDF for each kmer (computed across all regions)
import array
from collections import defaultdict
import sys
from cyvcf2 import VCF
from pyfaidx import Fasta, FetchError
from eskedit import base_methylation_probability
from eskedit.kclasses import MethylationModel, GRegion
from eskedit import is_dash
from eskedit import kmer_search
from eskedit.ktools import count_kmers, count_orfs, count_start_codons, count_stop_codons, count_strong_start_codons, \
    count_strong_orfs, count_g4s
from eskedit import is_quality_snv
from eskedit import is_cpg
import numpy as np
import multiprocessing as mp

print_lock = mp.Semaphore(value=1)


def calculate_expectation_matrix(region: GRegion, meth_vcf: VCF, models: MethylationModel, sequence: str,
                                 kmer_size: int = 7) -> tuple:
    base_methylation_probability_cache = base_methylation_probability
    cpg_check = is_cpg

    ct_freq_sum = 0.0

    raw_freq_sum = np.zeros((4, 4), dtype=np.float)
    adj_freq_sum = np.zeros((4, 4), dtype=np.float)

    nuc_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    for start_idx in range(len(sequence) - kmer_size + 1):
        forward_kmer = sequence[start_idx: start_idx + kmer_size]
        reverse_complement = ''.join(
            [complement.get(base, 'N') for base in list(sequence[start_idx: start_idx + kmer_size][::-1])])
        next_kmers = [forward_kmer, reverse_complement]

        # analysis for forward kmer
        from_idx = nuc_idx.get(forward_kmer[kmer_size // 2], None)
        if from_idx is None:
            continue
        new_expected = models.get_frequency('rare_transitions', forward_kmer)
        raw_freq_sum[from_idx] = raw_freq_sum[from_idx] + new_expected
        if cpg_check(forward_kmer):
            vcf_idx = region.start + start_idx
            meth_prob = base_methylation_probability_cache(meth_vcf(f'{region.chrom}:{vcf_idx}-{vcf_idx}'))
            ct_freq_meth = models.get_methylation_frequency(forward_kmer, meth_prob)[3]  # should be column index for T
            ct_freq_sum += ct_freq_meth
            new_expected[3] = ct_freq_meth

        adj_freq_sum[from_idx] = adj_freq_sum[from_idx] + new_expected

        # analysis for reverse kmer
        from_idx = nuc_idx.get(reverse_complement[kmer_size // 2], None)
        if from_idx is None:
            continue
        new_expected = models.get_frequency('rare_transitions', reverse_complement)
        raw_freq_sum[from_idx] = raw_freq_sum[from_idx] + new_expected
        if cpg_check(reverse_complement):
            vcf_idx = region.start + start_idx
            meth_prob = base_methylation_probability_cache(meth_vcf(f'{region.chrom}:{vcf_idx}-{vcf_idx}'))
            ct_freq_meth = models.get_methylation_frequency(reverse_complement, meth_prob)[3]  # should be column index for T
            ct_freq_sum += ct_freq_meth
            new_expected[3] = ct_freq_meth

        adj_freq_sum[from_idx] = adj_freq_sum[from_idx] + new_expected

        # for next_kmer in next_kmers:
        #     from_idx = nuc_idx.get(next_kmer[kmer_size // 2], None)
        #
        #     if from_idx is None:
        #         continue
        #     new_expected = models.get_frequency('rare_transitions', next_kmer)
        #     raw_freq_sum[from_idx] = raw_freq_sum[from_idx] + new_expected
        #     if cpg_check(next_kmer):
        #         vcf_idx = region.start + start_idx
        #         meth_prob = base_methylation_probability_cache(meth_vcf(f'{region.chrom}:{vcf_idx}-{vcf_idx}'))
        #         ct_freq_meth = models.get_methylation_frequency(next_kmer, meth_prob)[3]  # should be column index for T
        #         ct_freq_sum += ct_freq_meth
        #         new_expected[3] = ct_freq_meth
        #
        #     adj_freq_sum[from_idx] = adj_freq_sum[from_idx] + new_expected
    return raw_freq_sum, adj_freq_sum, ct_freq_sum


def calculate_expected(region: GRegion, meth_vcf: VCF, models: MethylationModel, sequence: str,
                       kmer_size: int = 7) -> tuple:
    # cach functions for better performance
    base_methylation_probability_cache = base_methylation_probability
    cpg_check = is_cpg
    raw_frequency_sum = np.zeros(4, dtype=np.float)
    adj_freq_sum = np.zeros(4, dtype=np.float)
    ct_freq = 0
    for start_idx in range(len(sequence) - kmer_size + 1):
        next_kmer = sequence[start_idx: start_idx + kmer_size]
        new_expected = models.get_frequency('rare_transitions', next_kmer)
        raw_frequency_sum = raw_frequency_sum + new_expected
        if cpg_check(next_kmer):
            vcf_idx = region.start + start_idx
            # vcf_idx = region.start - kmer_size//2 + start_idx + 2
            meth_prob = base_methylation_probability_cache(meth_vcf(f'{region.chrom}:{vcf_idx}-{vcf_idx}'))
            ct_freq_meth = models.get_methylation_frequency(next_kmer, meth_prob)[3]  # should be column index for T
            ct_freq += ct_freq_meth
            new_expected[3] = ct_freq_meth

        adj_freq_sum = adj_freq_sum + new_expected

    return raw_frequency_sum, adj_freq_sum, ct_freq


def query_regions(region_chunk: iter, models: MethylationModel, vcf_path: str, fasta_path: str, methylation_vcf: str):
    quality_snv = is_quality_snv
    cpg_check = is_cpg
    isdash = is_dash
    search_kmers = kmer_search
    base_methylation_probability_cache = base_methylation_probability

    gnomad_vcf = VCF(vcf_path)
    meth_vcf = VCF(methylation_vcf)
    fasta = Fasta(fasta_path, sequence_always_upper=True)
    kmer_size = models.kmer_size
    shift = kmer_size // 2
    # ref_kcounts = Counter()
    if kmer_size == 0:
        raise ValueError(f'Kmer size of {kmer_size} is not allowed. Revise model format.')
    for region in region_chunk:
        num_all_snvs = 0
        try:
            sequence = fasta.get_seq(region.chrom, region.start - shift, region.stop + shift).seq
            # sequence = region.get_seq_from_fasta(fasta, kmer_size=kmer_size)
            # if region.strand is not None:
            #     if isdash(region.strand):
            #         sequence = fasta.get_seq(region.chrom, region.start - shift,
            #                                  region.stop + shift).complement.seq[::-1]
            #     else:
            #         sequence = fasta.get_seq(region.chrom, region.start - shift, region.stop + shift).seq
            # else:
            #     sequence = fasta.get_seq(region.chrom, region.start - shift, region.stop + shift).seq
            # exp = window.calculate_expected(sequence)  # this does account for strandedness
            # AF, AC, AN, singletons, count = count_regional_alleles(vcf(str(region)))

        except (KeyError, FetchError):
            # TODO: meaningful error message
            continue

        if 'N' in sequence or sequence is None:
            continue
        kmer_counting = search_kmers(sequence, kmer_size, additional_functions=[count_kmers,
                                                                                count_orfs,
                                                                                count_strong_orfs,
                                                                                count_start_codons,
                                                                                count_strong_start_codons,
                                                                                count_stop_codons,
                                                                                count_g4s])
        ref_kcounts = kmer_counting['count_kmers']
        # frequency_array, adj_freq, ct_freq = calculate_expected(region, meth_vcf, models, sequence)
        freq_matrix, adj_freq_matrix, ct_freq = calculate_expectation_matrix(region, meth_vcf, models, sequence)
        nuc_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        no_count = 0
        num_ct = 0
        observed_context = defaultdict(lambda: array.array('L', [0, 0, 0, 0]))
        observed_transitions = np.zeros((4, 4), dtype=np.int)
        # measure tfidf and scale observed accordingly. this may address the high O/E
        for variant in gnomad_vcf(region.gnomad_rep()):
            varAC = variant.INFO.get('AC')
            if quality_snv(variant):
                num_all_snvs += 1
                from_idx = nuc_idx.get(variant.REF)
                to_idx = nuc_idx.get(variant.ALT[0])
                from_idx_comp = nuc_idx.get(complement.get(variant.REF, 'N'))
                to_idx_comp = nuc_idx.get(complement.get(variant.ALT[0], 'N'))
                observed_transitions[from_idx, to_idx] += 1
                observed_transitions[from_idx_comp, to_idx_comp] += 1
                # seq_idx = variant.POS - region.start + kmer_size // 2 - 1
                # seq_context = sequence[seq_idx - kmer_size // 2: seq_idx + kmer_size // 2 + 1]
                # if 'N' in seq_context or len(seq_context) == 0:
                #     continue
                # if cpg_check(seq_context) and variant.ALT[0] == 'T':
                #     num_ct += 1

        freq_string = '|'.join([str(num) for num in freq_matrix.reshape((16,))])
        adj_freq_str = '|'.join([str(num) for num in adj_freq_matrix.reshape((16,))])
        observed_str = '|'.join([str(num) for num in observed_transitions.reshape((16,))])
        region.add_field('expectation',
                         f'raw_expectation={freq_string};adjusted_expectation={adj_freq_str};expected_ct={ct_freq}')
        region.add_field('observed',
                         f'num_all_snvs={num_all_snvs};observed_transitions={observed_str};num_ct={num_ct}')
        region.add_field('predicted_structures',
                         f'num_start_codons={kmer_counting["count_start_codons"]};num_stop_codons={kmer_counting["count_stop_codons"]};num_orfs={kmer_counting["count_orfs"]};num_strong_start_codons={kmer_counting["count_strong_start_codons"]};num_strong_orfs={kmer_counting["count_strong_orfs"]};num_g4s={kmer_counting["count_g4s"][0]};g4_mfe={kmer_counting["count_g4s"][1]}')

        # print_lock.acquire()
        # print(str(region), flush=True)
        # print_lock.release()
        with print_lock:
            print(str(region), flush=True)
    pass


def kquery(bedpath: str, vcfpath: str, fastapath: str, kmer_size: int, methylation: str = None, header: bool = False,
           nprocs: int = 12, raw_counts: str = None, model_path: str = None):
    from eskedit import read_and_build_frequencies
    from eskedit.kclasses import MethylationModel
    from eskedit.ktools import read_models, get_bed_regions

    # read in model data
    if raw_counts is not None:
        models = read_and_build_frequencies(raw_counts)
    elif model_path is not None:
        models = read_models(model_path)
    else:
        raise ValueError("KQuery requires either raw_counts or model_path to not be \'None\'")
    model_operator = MethylationModel(models)
    # stores model names as a list of strings and adds all the models to the model operation
    model_names = model_operator.names

    # read input bed file, returns list of GRegions
    region_chunks = get_bed_regions(bedpath, nchunks=nprocs)

    args = []
    for region_chunk in region_chunks:
        # query_regions(region_chunk, model_operator, vcfpath, fastapath, methylation)
        args.append((region_chunk, model_operator, vcfpath, fastapath, methylation))

    with mp.Pool(nprocs) as pool:
        pool.starmap(query_regions, args)

    pass
