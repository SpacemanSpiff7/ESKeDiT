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

def get_methylation_probability(vcf_iter, name=None):
    if name is None:
        name = 'methyl_mean'
    meth_prob = None
    for variant in vcf_iter:
        # float(v.format('methylation')[0])
        meth_prob = float(variant.format('methylation')[0])
    if meth_prob is None:
        # perhaps change this return value to 'None' to distinguish from bases with data
        return 0.0
    else:
        return meth_prob


def is_quality_snv(variant):
    return variant.FILTER is None and len(variant.ALT) == 1 and len(variant.REF) == 1 and len(variant.ALT[0]) == 1


def count_kmers(sequence, kmer_length):
    counts = Counter()
    for i in range(len(sequence) - (kmer_length - 1)):
        next_seq = sequence[i:(i + kmer_length)]
        if not ('N' in next_seq):
            counts[next_seq.upper()] += 1
    return counts


def count_start_codons(sequence, kmer_length):
    start_codon = 'ATG'
    matches = re.findall(start_codon, sequence)
    return len(matches)


def count_stop_codons(sequence, kmer_length):
    regex = r'TAA|TAG|TGA'
    matches = re.findall(regex, sequence)
    return len(matches)


def count_g4s():
    # TODO
    pass


def count_orfs(sequence, kmer_length):
    # regex = r'(?<=ATG)(.*)(?=TAA|TAG|TGA)'
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    in_orf = False
    orf_start = None
    orfs = []
    for idx, base in enumerate(sequence):
        if not in_orf:
            putative_codon = sequence[idx:idx + 3]
            if putative_codon == start_codon:
                in_orf = True
                orf_start = idx
                idx += 2
        else:
            next_codon = sequence[idx:idx + 3]
            if next_codon in stop_codons:
                orfs.append((orf_start, idx + 3))
                in_orf = False
                orf_start = idx + 3
            idx += 2
    # if in_orf:
    #     orfs.append((orf_start, None))

    return len(orfs)


def kmer_search(sequence, kmer_length, additional_functions=None):
    """
    Driver for get_kmer_count
    :param additional_functions: an iterable consisting of functions that accept sequence and kmer_length (in that order) as positional arguments
    :param sequence:
    :param kmer_length:
    :param count_gc: Return GC content (default = False)
    :param count_n: Reurn count number of unreferenced nucleotides (default = False)
    :return: a list containing the results of the analysis (will always count_kmers by default)
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


def get_bed_regions(bedpath, hasheader=False, nchunks=None, verbose=False):
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
        if not isinstance(nchunks, int):
            return regions
        else:
            nchunks = int(nchunks)
            regions = np.array_split(regions, nchunks)

    if verbose:
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
