#!/bin/bash 

conda activate kmer
conda install -c bioconda jellyfish

# jellyfish count 
jellyfish count -C -m 21 -s 600M -t 32 *.fastq -o Phar_reads.jf

# generate kmer histogram
jellyfish histo -t 32 Phar_reads.jf > Phar_reads.histo

# # GENOMESCOPE # #
Rscript genomescope.R Phar_reads.histo 21 150 Phar_genomescope_21

# # SMUDGEPLOT # # 
L=$(smudgeplot.py cutoff Phar_03_reads.histo L)
U=$(smudgeplot.py cutoff Phar_03_reads.histo U)

echo $L $U # these need to be sane values like 30 800 or so

jellyfish dump -c -L $L -U $U Phar_reads.jf | smudgeplot.py hetkmers -o kmer_pairs

smudgeplot.py plot kmer_pairs_coverages.tsv -o Phar_smudgeplot
