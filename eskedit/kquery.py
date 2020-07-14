# Steps:
# 1. load model - multinom or freq?
#     a. use raw count data
#     b. option for probability
# 2. sliding window
# 3. input bed, vcf, fasta
# for region in bed
#     count kmers
#     output TFIDF for each kmer (computed across all regions)
from cyvcf2 import VCF
from pyfaidx import Fasta, FetchError

from eskedit.kclasses import ModelOps


def query_regions(region_chunk: iter, models: ModelOps, vcf_path: str, fasta_path: str):
    from eskedit import is_dash
    vcf = VCF(vcf_path)
    fasta = Fasta(fasta_path, sequence_always_upper=True)
    kmer_size = models.kmer_size
    shift = kmer_size // 2
    if kmer_size == 0:
        raise ValueError(f'Kmer size of {kmer_size} is not allowed. Revise model format.')
    for region in region_chunk:
        try:
            if region.strand is not None:
                if is_dash(region.strand):
                    sequence = fasta.get_seq(region.chrom, region.start - shift,
                                             region.stop + shift).complement.seq.upper()
                else:
                    sequence = fasta.get_seq(region.chrom, region.start - shift, region.stop + shift).seq.upper()
            else:
                sequence = fasta.get_seq(region.chrom, region.start - shift, region.stop + shift).seq.upper()
            # exp = window.calculate_expected(sequence)  # this does account for strandedness
            # AF, AC, AN, singletons, count = count_regional_alleles(vcf(str(region)))

        except (KeyError, FetchError):
            # TODO: meaningful error message
            continue
        continue
    pass


def kquery(bedpath: str, vcfpath: str, fastapath: str, kmer_size: int, methylation: bool = True, header: bool = False,
           nprocs: int = 6, raw_counts: str = None, model_path: str = None):
    from eskedit import read_and_build_models
    from eskedit.kclasses import ModelOps
    from eskedit.ktools import read_models, get_bed_regions

    # read in model data
    if raw_counts is not None:
        models = read_and_build_models(raw_counts)
    elif model_path is not None:
        models = read_models(model_path)
    else:
        raise ValueError("KQuery requires either raw_counts or model_path to not be \'None\'")
    model_operator = ModelOps()
    # stores model names as a list of strings and adds all the models to the model operation
    model_names = model_operator.add_models(models)

    # read input bed file, returns list of GRegions
    region_chunks = get_bed_regions(bedpath, nchunks=nprocs)

    for region_chunk in region_chunks:
        query_regions(region_chunk, model_operator, vcfpath, fastapath)
        continue
    pass
