import re
import itertools
import random
import multiprocessing as mp
from collections import Counter
from pyfaidx import Fasta
import numpy as np
from eskedit import GRegion
import os
import datetime
import pandas as pd
from os import listdir
from os.path import isfile, join, splitext
from scipy.stats import multinomial
import math
import RNA
from eskedit.kclasses import Model, ModelFreq


def rowmnd(row):
    ns = [r for r in row.iloc[:4] if r > 0]
    ns.append(row.iloc[4] - sum(ns))
    alphas = [n / row.iloc[4] for n in ns]
    return multinomial.pmf(ns, n=row.iloc[4], p=alphas)


def generate_frequencies(transitions_path, counts_path):
    ts = pd.read_csv(transitions_path, index_col=0).sort_index()
    cts = pd.read_csv(counts_path, index_col=0).sort_index()
    cts.columns = ['counts']
    # merge to align indices
    ts = ts.merge(cts, left_on=None, right_on=None, left_index=True, right_index=True)

    # commented out becuase doesn't consider all kmers that mutate
    # ts = ts.iloc[:, :4].div(ts.iloc[:, 4], axis=0)
    def row_freq_calc(row: pd.Series):
        row = row.astype(np.float128)
        for i, v in row.iteritems():
            row[i] = v / max((row[-1] - v), 0.1)
        return row

    ts = ts.apply(row_freq_calc, axis=1)
    # ts = ts.iloc[:, :4].div((ts.iloc[:, 4] - ts.iloc[:, :4].apply(sum, axis=1)), axis=0)
    return ts.iloc[:, :4]


def generate_multinomial_probabilities(transitions_path, counts_path, ksize=None):
    ts = pd.read_csv(transitions_path, index_col=0).sort_index()
    cts = pd.read_csv(counts_path, index_col=0).sort_index()
    cts.columns = ['counts']
    ts = ts.merge(cts, left_on=None, right_on=None, left_index=True, right_index=True)
    if ksize is None:
        ksize = int(math.log(len(cts), 4))
    ts['probability'] = ts.apply(rowmnd, axis=1)
    return ts[['probability']].fillna(0)


def read_and_build_models(directory_path: str) -> list:
    # https://stackoverflow.com/questions/3207219/how-do-i-list-all-files-of-a-directory
    expected_keys = ['ref_kmer_counts', 'singleton_transitions', 'hi_methylation', 'intermediate_methylation',
                     'low_methylation']
    files = {}
    for f in listdir(directory_path):
        if isfile(join(directory_path, f)):
            files.update({splitext(f)[0]: join(directory_path, f)})
    counts = {}
    ref_counts_path = files[expected_keys[0]]
    for name, path in files.items():
        if name != expected_keys[0]:
            counts.update({name: generate_multinomial_probabilities(path, ref_counts_path)})
    return df_to_kmodel(counts)


def read_and_build_frequencies(directory_path: str) -> list:
    # https://stackoverflow.com/questions/3207219/how-do-i-list-all-files-of-a-directory
    expected_keys = ['ref_kmer_counts', 'singleton_transitions', 'hi_methylation', 'intermediate_methylation',
                     'low_methylation']
    files = {}
    for f in listdir(directory_path):
        if isfile(join(directory_path, f)):
            files.update({splitext(f)[0]: join(directory_path, f)})
    counts = {}
    ref_counts_path = files[expected_keys[0]]
    for name, path in files.items():
        if name != expected_keys[0]:
            counts.update({name: generate_frequencies(path, ref_counts_path)})
    return df_to_freq_model(counts)


def read_models(directory_path: str) -> list:
    files = {}
    for f in listdir(directory_path):
        if isfile(join(directory_path, f)):
            files.update({splitext(f)[0]: join(directory_path, f)})
    counts = {}
    for name, path in files.items():
        counts.update({name: pd.read_csv(path, index_col=0)})
    return df_to_kmodel(counts)


def df_to_kmodel(counts: dict) -> list:
    models = []
    for k, v in counts.items():
        models.append(Model(k, v))
    return models


def df_to_freq_model(counts: dict) -> list:
    models = []
    for k, v in counts.items():
        models.append(ModelFreq(k, v))
    return models


def is_cpg_ct(seq: str, alt: str) -> bool:
    k = len(seq)
    if k < 3:
        if k == 2:
            return seq[0] == 'C' and seq[1] == 'G' and alt == 'T'
        else:
            return False
    elif seq[k // 2] == 'C' and seq[k // 2 + 1] == 'G' and alt == 'T':
        return True
    else:
        return False


def is_cpg(seq: str) -> bool:
    if seq is None:
        return False
    k = len(seq)
    if k < 3:
        if k == 2:
            return seq[0] == 'C' and seq[1] == 'G'
        else:
            return False
    else:
        return seq[k // 2] == 'C' and seq[k // 2 + 1] == 'G' and 'N' not in seq


def avg_seq_methylation_probability(vcf_iter, seq_len: int = None) -> float:
    meth_probs = []
    for variant in vcf_iter:
        # float(v.format('methylation')[0])
        meth_probs.append(float(variant.format('methylation')[0]) / 100)

    if seq_len is not None:
        len_meth_prob = seq_len
    else:
        len_meth_prob = len(meth_probs)
    if len(meth_probs) == 0:
        # perhaps change this return value to 'None' to distinguish from bases with data
        return 0.0
    else:
        return sum(meth_probs) / len_meth_prob


def base_methylation_probability(vcf_iter, name=None) -> float:
    if name is None:
        name = 'methyl_mean'
    meth_prob = None
    for variant in vcf_iter:
        # float(v.format('methylation')[0])
        # print(f'REF: {variant.REF}\tALT: {variant.ALT}\t{variant.POS}')
        meth_prob = float(variant.format('methylation')[0]) / 100
    if meth_prob is None:
        # perhaps change this return value to 'None' to distinguish from bases with data
        return -1.0
    else:
        return meth_prob


def is_quality_snv(variant) -> bool:
    return variant.FILTER is None and len(variant.ALT) == 1 and len(variant.REF) == 1 and len(variant.ALT[0]) == 1


def count_kmers(sequence: str, kmer_length: int) -> Counter:
    counts = Counter()
    for i in range(len(sequence) - (kmer_length - 1)):
        next_seq = sequence[i:(i + kmer_length)]
        if not ('N' in next_seq):
            counts[next_seq.upper()] += 1
    return counts


def count_start_codons(sequence: str, kmer_length: int) -> int:
    start_codon = re.compile(r'ATG')
    matches = start_codon.findall(sequence)
    return len(matches)


def count_stop_codons(sequence: str, kmer_length: int) -> int:
    regex = re.compile(r'TAA|TAG|TGA')
    matches = regex.findall(sequence)
    return len(matches)


def get_MFE(seq: str):
    ss, mfe = RNA.fold(seq)
    return mfe


def g4s_sliding_window(seq):
    num_consecutive_gs = 0
    num_interg_nucs = 0
    in_pg4 = False
    for idx in range(len(seq)):
        if seq[idx] == 'G':
            num_consecutive_gs += 1
            if num_consecutive_gs >= 2:
                in_pg4 = True
            else:
                in_pg4 = False
            continue
        else:
            num_consecutive_gs = 0
        if in_pg4:
            num_interg_nucs += 1
            if num_interg_nucs > 7:
                in_pg4 = False


def count_g4s(sequence: str, kmer_length: int) -> tuple:
    g4_regex = re.compile('([G]{2,5}[ACGT]{1,7}){3}[G]{2,5}')
    g4_regex = re.compile(r'([G]{2,5}\w{1,7}){3}[G]{2,5}', re.UNICODE)
    count = 0
    pg4s = []
    for match in g4_regex.finditer(sequence):
        count += 1
        start, stop = match.span()
        pg4 = sequence[start:stop]
        # print(pg4, get_MFE(pg4))
        pg4s.append(':'.join([pg4, str(get_MFE(pg4))]))

    # if count == 0:
    #     return 0, 'none'
    return count, '|'.join(pg4s)


def _kozak_strength(sequence: str) -> int:
    # 3: Strong = (A/G)NNAUGG
    # 2: Moderate = (A/G)NNAUGN or NNNAUGG
    # 1: Weak = NNNAUGN
    # 0: Nothing = NNNNNNN
    if len(sequence) != 7:
        return 0
    else:
        score = 0
        putative_start_codon = sequence[3:6]
        if putative_start_codon == 'ATG':
            score += 1
        else:
            return 0
        if sequence[0] == 'A' or sequence[0] == 'G':
            score += 1
        if sequence[6] == 'G':
            score += 1
        return score


def count_strong_start_codons(sequence: str, kmer_length: int) -> int:
    # start_codon = re.compile(r'ATG')
    count = 0
    kozak_str = _kozak_strength

    for cdn_start in range(len(sequence)):
        next_codon = sequence[cdn_start:cdn_start + 3]
        if next_codon == 'ATG':
            if cdn_start - 3 < 0 or cdn_start + 4 > len(sequence) - 1:
                continue
            else:
                context = sequence[cdn_start - 3:cdn_start + 4]
                if kozak_str(context) == 3:
                    count += 1
    return count


def count_orfs(sequence: str, kmer_length: int) -> int:
    # regex = r'(?<=ATG)(.*)(?=TAA|TAG|TGA)'
    # TODO: noncanonical start codons and all reading frames
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    in_orf = False
    orf_start = None
    # kozak = 0
    orfs = []
    for idx, base in enumerate(sequence):
        if not in_orf:
            putative_codon = sequence[idx:idx + 3]
            if putative_codon == start_codon:
                in_orf = True
                # kozak = _kozak_strength(sequence[idx - 3:idx + 4])
                orf_start = idx
                idx += 2
        else:
            next_codon = sequence[idx:idx + 3]
            if next_codon in stop_codons:
                orfs.append((orf_start, idx + 3))
                in_orf = False
                # kozak = 0
                orf_start = idx + 3
            idx += 2
    # if in_orf:
    #     orfs.append((orf_start, None))

    return len(orfs)


def count_strong_orfs(sequence: str, kmer_length: int) -> int:
    # regex = r'(?<=ATG)(.*)(?=TAA|TAG|TGA)'
    # TODO: noncanonical start codons
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    kozak = 0
    orfs = []
    check_kozak = _kozak_strength
    # split seq into 3 reading frames
    reading_frames = [sequence[i:] for i in range(3)]
    for seq in reading_frames:
        in_orf = False
        orf_start = None

        for codon_idx in range(0, len(seq), 3):
            next_codon = seq[codon_idx:codon_idx + 3]
            if not in_orf:
                if next_codon == start_codon:
                    kozak = check_kozak(seq[codon_idx - 3:codon_idx + 4])
                    if kozak == 3:
                        in_orf = True
                        orf_start = codon_idx
            else:
                if next_codon in stop_codons:
                    orfs.append((orf_start, codon_idx + 3))
                    in_orf = False
    return len(orfs)
    # for idx, base in enumerate(sequence):
    #     if not in_orf:
    #         putative_codon = sequence[idx:idx + 3]
    #         if putative_codon == start_codon:
    #             in_orf = True
    #             kozak = _kozak_strength(sequence[idx - 3:idx + 4])
    #             orf_start = idx
    #             idx += 2
    #     else:
    #         next_codon = sequence[idx:idx + 3]
    #         if next_codon in stop_codons:
    #             orfs.append((orf_start, idx + 3))
    #             in_orf = False
    #             # kozak = 0
    #             orf_start = idx + 3
    #         idx += 2
    # # if in_orf:
    #     orfs.append((orf_start, None))


def kmer_search(sequence: str, kmer_length: int, additional_functions: iter = None) -> dict:
    """
    Driver for get_kmer_count
    :param additional_functions: an iterable consisting of functions that accept sequence and kmer_length (in that order) as positional arguments
    :param sequence:
    :param kmer_length:
    :return: a dictionary containing the results of the analysis (will always count_kmers by default). Maps function name to result
    """
    if additional_functions is None:
        additional_functions = []
    if count_kmers not in additional_functions:
        additional_functions.append(count_kmers)
    results = {}

    if additional_functions is not None:
        for custom_function in additional_functions:
            results.update({custom_function.__name__: custom_function(sequence, kmer_length)})

    return results
    # counts = Counter()
    # for i in range(len(sequence) - (kmer_length - 1)):
    #     next_seq = sequence[i:(i + kmer_length)]
    #     if not ('N' in next_seq or 'n' in next_seq):
    #         counts[next_seq.upper()] += 1
    # if count_gc and count_n:
    #     nucleotides = Counter(sequence)
    #     gc_content = (nucleotides['G'] + nucleotides['C']) / max(
    #         (nucleotides['A'] + nucleotides['T'] + nucleotides['G'] + nucleotides['C']), 1)
    #     return counts, gc_content, nucleotides['N']
    # elif count_gc and not count_n:
    #     nucleotides = Counter(sequence)
    #     gc_content = (nucleotides['G'] + nucleotides['C']) / max((
    #             nucleotides['A'] + nucleotides['T'] + nucleotides['G'] + nucleotides['C']), 1)
    #     return counts, gc_content
    # elif not count_gc and count_n:
    #     nucleotides = Counter(sequence)
    #     return counts, nucleotides['N']
    # else:
    #     return counts


def get_bed_regions(bedpath, hasheader: bool = False, nchunks: int = None, return_extra=False):
    headernames = None
    numfields = 0
    regions = []
    chroms = set()
    # the number of fields less than or equal to  6  must be positional in concordance with bed format
    # additional fields are stored  as keyword arguments
    with open(bedpath, 'r') as bed:
        if hasheader:
            # even though not used, we still need to read a line if a header is present, so it stays! (for now)
            headernames = bed.readline().split('\t')
        for line in bed.readlines():
            fields = [s.strip() for s in line.split('\t') if len(s) > 0]
            if len(fields) > 6:  # 6 is the maximumum number of positional arguments for GRegion
                kwargs = {}
                addcount = 1
                for addlfield in fields[6:]:
                    addstring = 'field{}'.format(addcount)
                    kwargs.update({addstring: addlfield})
                    addcount += 1
                new_region = GRegion(*fields[:6], **kwargs)
                regions.append(new_region)
            else:
                new_region = GRegion(*fields)
                regions.append(new_region)
            chroms.add(new_region.chrom)
            numfields = max(numfields, new_region.num_fields())

        if headernames is None:
            headernames = [str(i) for i in range(1, numfields + 1)]

    if nchunks is None:
        return regions
    else:
        nchunks = int(nchunks)
        regions = np.array_split(regions, nchunks)

    if return_extra:
        return regions, chroms, headernames
    else:
        return regions


def ran_seq(seq_len):
    """
    Driver method for gen_random_sequence
    :param seq_len:
    :return: random lsequence of length seq_len
    """
    sequence = ""
    for i in range(seq_len):
        sequence += random.choice('ACTG')
    return sequence


def generate_kmers(k):
    """Generates a list of all possible DNA sequences of length k. E.g. generate_kmers(2) will return
    [ AA, AC, AT, AG, CA, CC, CT, CG, TA, TC, TT, TG, GA, GC, GT, GG ] """
    len_k = int(k)
    if len_k < 0:
        raise ValueError("Must be a positive integer")
    combos = list(itertools.product('ACTG', repeat=len_k))
    seqs = []
    for seq in combos:
        seqs.append(''.join(seq))
    return seqs


def gen_random_sequence(length):
    """
    :param length: (integer) number of nucleotides in a random sequence
    :return: a random DNA sequence of length 'length'
    """
    if length > 5_000_000:
        pool = mp.Pool()
        nprocs = mp.cpu_count()
        chunk_size = int(length / nprocs) + 1
        args = []
        tot = chunk_size
        while tot <= length:
            args.append(chunk_size)
            tot += chunk_size
        new_tot = np.array(args).sum()
        if new_tot < length:
            args.append(length - new_tot)
        results = [funccall.get() for funccall in [pool.map_async(ran_seq, args)]]
        pool.close()
        random_seq = ""
        for seq in results[0]:
            random_seq += seq
    else:
        random_seq = ran_seq(length)
    return random_seq


def ref_genome_as_string(ref_fasta, keys=None):
    ref_genome = Fasta(ref_fasta)
    if keys is None:
        keys = ref_genome.keys()
    ref_seq = ""
    for chrom in keys:
        ref_seq += str(ref_genome[chrom])
    return ref_seq


def is_dash(pdash):
    regex = '[\u002D\u058A\u05BE\u1400\u1806\u2010-\u2015\u2E17\u2E1A\u2E3A\u2E3B\u2E40\u301C\u3030\u30A0\uFE31\uFE32\uFE58\uFE63\uFF0D]'
    if re.match(regex, pdash) is not None:
        return True


def make_directory(dirname='Model_{}'.format(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))):
    # https://stackoverflow.com/questions/14115254/creating-a-folder-with-timestamp/14115286
    mydir = os.path.join(os.getcwd(), dirname)
    tempdir = mydir
    count = 1
    while True:
        try:
            os.makedirs(mydir)
            break
        except FileExistsError:
            mydir = tempdir + '_' + str(count)
            count += 1
            continue
    return mydir
