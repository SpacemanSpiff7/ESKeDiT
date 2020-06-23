"""
ESKeDiT

Created by Simone Longo at the University of Utah 2019
"""
from eskedit import runtests
import argparse

if __name__ == "__main__":
    main_parser = argparse.ArgumentParser(description='ESKeDiT version 2.0.0')
    main_parser.add_argument('-b', '--bed', action='store', dest='bed_file_path', help='Input BED file path')
    main_parser.add_argument('-v', '--vcf', action='store', dest='vcf_file_path',
                             help='Input file path for variant VCF')
    main_parser.add_argument('-m', '--methylation_vcf', action='store', dest='meth_vcf_path',
                             help='Input path for methylation VCF (optional)')
    main_parser.add_argument('-f', '--fasta', action='store', dest='fasta_path',
                             help='Input path the reference fasta file')
    main_parser.add_argument('-@', '--nprocs', action='store', dest='nprocs', help='Number of available CPUs')
    main_parser.add_argument('-k', '--kmer_size', action='store', dest='kmer_size', help='K-mer size to use')
    main_parser.add_argument('-t', action='store_true', dest='run_test')
    main_parser.add_argument('-u', action='store_true', dest='unit_test')
    cmd_args = main_parser.parse_args()
    from eskedit import ktrain

    # This block is for testing purposes
    if cmd_args.run_test:
        nprocs = 6
        bedpath = '/Users/simonelongo/Documents/QuinlanLabFiles/ESKeDiT/notebooks/notebook_resources/pc_exon_complement_22june2020.bed'
        vcf_path = "/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomad.genomes.r3.0.sites.vcf.bgz"
        fasta_path = "/Users/simonelongo/too_big_for_icloud/ref_genome/hg38/hg38.fa"
        meth_vcf_path = "/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomadv3_methylation_2.vcf.bgz"
        ktrain(bedpath, vcf_path, fasta_path, 7, meth_vcf_path=meth_vcf_path, nprocs=nprocs)
        exit(0)
    if cmd_args.unit_test:
        runtests()
        exit(0)

    if cmd_args.meth_vcf_path is None:
        # run without methylation
        ktrain(cmd_args.bed_file_path, cmd_args.vcf_file_path, cmd_args.fasta_path, 7, meth_vcf_path=None, nprocs=nprocs)
        pass
    else:  # run with metylation
	ktrain(cmd_args.bed_file_path, cmd_args.vcf_file_path, cmd_args.fasta_path, 7, meth_vcf_path=cmd_args.meth_vcf_path, nprocs=cmd_args.nprocs)
        pass
