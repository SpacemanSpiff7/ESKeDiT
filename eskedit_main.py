"""
ESKeDiT

Created by Simone Longo at the University of Utah 2019
"""
import argparse


def ktrain_driver(args):
    if args.vcf_file_path is None or \
            args.bed_file_path is None or \
            args.meth_vcf_path is None or \
            args.fasta_path is None or \
            args.kmer_size is None:
        raise ValueError('All file paths must be provided.')
    from eskedit import ktrain
    ktrain(args.bed_file_path,
           args.vcf_file_path,
           args.fasta_path,
           int(args.kmer_size),
           meth_vcf_path=args.meth_vcf_path,
           nprocs=args.nprocs)
    pass


def kquery_driver(args):
    if args.vcf_file_path is None or \
            args.bed_file_path is None or \
            args.meth_vcf_path is None or \
            args.fasta_path is None or \
            args.kmer_size is None or \
            args.counts_directory_path is None:
        raise ValueError('All file paths must be provided.')
    from eskedit import kquery
    kquery(args.bed_file_path,
           args.vcf_file_path,
           args.fasta_path,
           args.kmer_size,
           methylation=args.meth_vcf_path,
           raw_counts=args.counts_directory_path,
           nprocs=args.nprocs)
    return


def test_driver(args):
    from eskedit import runtests
    runtests()
    return


if __name__ == "__main__":
    from eskedit import ktrain, kquery

    function_map = {'train': ktrain_driver,
                    'query': kquery_driver,
                    'test': test_driver}

    main_parser = argparse.ArgumentParser(description='ESKeDiT version 2.0.0')
    main_parser.add_argument('command', choices=function_map.keys())
    main_parser.add_argument('-b', '--bed', action='store', dest='bed_file_path', help='Input BED file path')
    main_parser.add_argument('-v', '--vcf', action='store', dest='vcf_file_path',
                             help='Input file path for variant VCF')
    main_parser.add_argument('-m', '--methylation_vcf', action='store', dest='meth_vcf_path',
                             help='Input path for methylation VCF (optional)')
    main_parser.add_argument('-f', '--fasta', action='store', dest='fasta_path',
                             help='Input path the reference fasta file')
    main_parser.add_argument('-o', '--output', action='store', dest='output_path',
                             help='Directory path for outpile files.')
    main_parser.add_argument('-@', '--nprocs', action='store', dest='nprocs', help='Number of available CPUs',
                             default=1)
    main_parser.add_argument('-k', '--kmer_size', action='store', dest='kmer_size', help='K-mer size to use', default=7)
    # main_parser.add_argument('-t', action='store_true', dest='run_test')
    # main_parser.add_argument('-u', action='store_true', dest='unit_test')
    main_parser.add_argument('--counts', action='store', dest='counts_directory_path')
    cmd_args = main_parser.parse_args()

    function = function_map.get(cmd_args.command, None)
    if function is not None:
        function(cmd_args)

    # This block is for testing purposes
    # if cmd_args.run_test:
    #     nprocs = 6
    #     # bedpath = '/Users/simonelongo/Documents/QuinlanLabFiles/ESKeDiT/notebooks/notebook_resources/pc_exon_complement_22june2020.bed'
    #     bedpath = "/Users/simonelongo/Documents/QuinlanLabFiles/ESKeDiT/resources/testfiles/variant_rich.bed"
    #
    #     vcf_path = "/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomad.genomes.r3.0.sites.vcf.bgz"
    #     fasta_path = "/Users/simonelongo/too_big_for_icloud/ref_genome/hg38/hg38.fa"
    #     meth_vcf_path = "/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomadv3_methylation_2.vcf.bgz"
    #     ktrain(bedpath, vcf_path, fasta_path, 7, meth_vcf_path=meth_vcf_path, nprocs=nprocs)
    #     exit(0)
    # if cmd_args.unit_test:
    #     from eskedit import runtests
    #
    #     runtests()
    #     exit(0)

    # if cmd_args.meth_vcf_path is None:
    #     # run without methylation
    #     ktrain(cmd_args.bed_file_path, cmd_args.vcf_file_path, cmd_args.fasta_path, 7, meth_vcf_path=None,
    #            nprocs=nprocs)
    #     pass
    # else:  # run with metylation
    #     ktrain(cmd_args.bed_file_path, cmd_args.vcf_file_path, cmd_args.fasta_path, 7,
    #            meth_vcf_path=cmd_args.meth_vcf_path, nprocs=cmd_args.nprocs)
    #     pass
