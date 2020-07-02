import time
from collections import Counter

from cyvcf2 import VCF


def test_ktrain():
    start = time.time()
    from eskedit import GRegion, ktrain
    #test_region = GRegion('chr1', 1361264.0, 1361777.0)
    gnomad_vcf = VCF("/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomad.genomes.r3.0.sites.vcf.bgz")
    vcf_path = "/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomad.genomes.r3.0.sites.vcf.bgz"
    fasta_path = "/Users/simonelongo/too_big_for_icloud/ref_genome/hg38/hg38.fa"
    # bed_path = '/Users/simonelongo/Documents/QuinlanLabFiles/ESKeDiT/resources/testfiles/test_ENSEMBL.bed'
    bed_path = "/Users/simonelongo/Documents/QuinlanLabFiles/ESKeDiT/resources/testfiles/variant_rich.bed"
    meth_vcf_path = "/Users/simonelongo/too_big_for_icloud/gnomAD_v3/gnomadv3_methylation_2.vcf.bgz"
    results = ktrain(bed_path, vcf_path, fasta_path, 7, meth_vcf_path)
    print(results)
    #var_counts = Counter()
    # for variant in gnomad_vcf(test_region.gnomad_rep()):
    #     if variant.FILTER is None:
    #         var_counts[variant.INFO.get('AC')] += 1
    # print(var_counts)
    return 'ktrain test done in {}'.format(time.time() - start)


def test_ktools():
    start = time.time()
    # generate kmers
    for i in range(1, 10):
        from eskedit import ktools
        assert len(ktools.generate_kmers(i)) == (4 ** i), "Generated wrong number of {}-mers".format(i)
    print('Generate kmers test passed')

    testbedpath = '/Users/simonelongo/Documents/QuinlanLabFiles/ESKeDiT/resources/testfiles/small_test.bed'
    from eskedit import get_bed_regions
    print(get_bed_regions(testbedpath))

    t1 = time.time()
    test3fields = get_bed_regions(
        '/Users/simonelongo/Documents/QuinlanLabFiles/ESKeDiT/resources/testfiles/test3fields.bed')
    time_t1 = time.time() - t1
    assert len(test3fields) == 100

    t1 = time.time()
    testENSEMBL = get_bed_regions(
        '/Users/simonelongo/Documents/QuinlanLabFiles/ESKeDiT/resources/testfiles/test_ENSEMBL.bed')
    time_t2 = time.time() - t1
    assert len(testENSEMBL) == 32

    t1 = time.time()
    testBIG = get_bed_regions(
        '/Users/simonelongo/Documents/QuinlanLabFiles/ESKeDiT/resources/testfiles/big_test.bed')
    time_t3 = time.time() - t1
    assert len(testBIG) == 19027

    # t1 = time.time()
    testnchunks = 12
    testBIG = get_bed_regions(
        '/Users/simonelongo/Documents/QuinlanLabFiles/ESKeDiT/resources/testfiles/big_test.bed', nchunks=testnchunks)
    # time_t3 = time.time() - t1
    assert len(testBIG) == testnchunks
    totallen = sum([len(reg) for reg in testBIG])
    # print('Length testBIG: {} Total Length: {}'.format(len(testBIG), totallen))
    assert totallen == 19027

    print(
        'get_bed_regions passed. time for reading 100 lines 3 fields: {}\n100 lines with n fields: {}\n19,027 lines with 5 fields: {}'.format(
            time_t1, time_t2, time_t3))
    print()

    # testing kmer_search
    from eskedit.ktools import gen_random_sequence
    sequence = gen_random_sequence(10000)
    s1 = 'GAGAGATGCGCGCGCGCGTAAGCGCGGCATGCGGCGGTAGGGGATGGGGGGGG'
    s2 = 'GAGAGATGGCGCGCGCGCGTAAGCGCGGCATGCGGCGGTAGGGGATGGGGGGGG'
    from eskedit.ktools import kmer_search
    from eskedit.ktools import count_kmers
    from eskedit.ktools import count_start_codons, count_stop_codons, count_orfs

    ktime = time.time()
    results1 = kmer_search(s1, 3, additional_functions=[count_kmers, count_start_codons, count_stop_codons, count_orfs])
    assert results1['count_orfs'] == 2
    assert results1['count_start_codons'] == 3
    assert results1['count_stop_codons'] == 2
    # for k, v in results1.items():
    #     print(k, str(v))
    print(f'kmer_search wall time: {time.time() - ktime}')
    return 'ktools test done in {}'.format(time.time() - start)


def test_kquery():
    start = time.time()
    print()
    return 'kquery test done in {}'.format(time.time() - start)


def test_kclasses():
    start = time.time()

    from eskedit import GRegion
    test_region = GRegion('chr1', 5646, 987987, 'bleh', 'greh', '+', info='blablabla')
    print(str(test_region))
    test_region = GRegion('chr1', 5646, 987987, info='blablabla')
    print(str(test_region))
    test_region = GRegion('chr1', 5646, 987987, 'bleh', 'greh')
    print(str(test_region))
    test_region = GRegion('chr1', 5646, 987987)
    print(str(test_region))
    test_region = GRegion('chr1', 5646, 987987, info='blablabla')
    print(str(test_region))

    print('GRegion constructor runs! Check output.')
    print()
    return 'kclasses test done in {}'.format(time.time() - start)


def runtests():
    """
    This class runs tests listed in the 'results' declaration and prints the results
    :return:
    """
    start = time.time()
    results = [test_ktrain(),
               test_ktools(),
               test_kquery(),
               test_kclasses()]

    for result in results:
        print(str(result))

    from eskedit import constants
    print("\n{}\nAll tests done in {}".format(constants.ESKEDIT_STRING, (time.time() - start)))
