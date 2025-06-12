#!/bin/bash 

conda activate kmer
conda install bioconda::meryl

# count kmers
meryl count k=21 output Phar_porechopped_reads.meryl Phar_porechopped_reads.fastq.gz

# generate kmer histogram
meryl histogram Phar_porechopped_reads.meryl > Phar_porechopped_reads.meryl.hist

# # GENOMESCOPE # #
Rscript genomescope.R Phar_porechopped_reads.meryl.hist 21 150 Phar_genomescope_21
