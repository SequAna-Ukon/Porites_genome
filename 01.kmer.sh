#!/bin/bash 

conda activate kmer
conda install bioconda::meryl

# count kmers
meryl count k=21 output Phar_reads_porechopped.meryl Phar_reads_porechopped.fastq.gz

# generate kmer histogram
meryl histogram Phar_reads_porechopped.meryl > Phar_meryldb.hist

# # GENOMESCOPE # #
Rscript genomescope.R Phar_meryldb.hist 21 150 Phar_genomescope_21
