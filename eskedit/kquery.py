"""
Steps:
1. load model - multinom or freq?
    a. use raw count data
    b. option for probability
2. sliding window
3. input bed, vcf, fasta
for region in bed
    count kmers
    output TFIDF for each kmer (computed across all regions)
"""