#!/bin/bash 

conda activate kmer
conda install bioconda::meryl

# count kmers
meryl count k=21 output PharDB.meryl ONT_reads_pc.fastq.gz

# generate kmer histogram
meryl histogram PharDB.meryl > PharDB.meryl.hist

# # GENOMESCOPE # #
Rscript genomescope.R PharDB.meryl.hist 21 150 Phar_genomescope_21
