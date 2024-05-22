#!/bin/bash
# ONT reads preprocessing
# -----------------------

# # ENVIRONMENTS & WORKDIRECTORY # #
workdir="/PATH/TO/READS"
ref_genomes="/PATH/TO/SYMBIONT_GENOME_REF"

mamba activate reads_env
mamba install -c bioconda porechop
mamba install -c bioconda fastqc
mamba install -c bioconda multiqc
mamba install -c bioconda chopper

# # ADAPTER REMOVAL USING PORECHOP # #
porechop --input $workdir/ONT_raw_reads.fastq.gz -o ONT_reads_pc.fastq.gz --discard_middle

# # QUALITY CONTROL # #
fastqc -t 8 ONT_reads_pc.fastq.gz
multiqc ONT_reads_pc_fastqc.zip -o multiqc --interactive
NanoPlot -t 10 --fastq ONT_reads_pc.fastq.gz --plots dot --legacy hex --N50

# # FILTER & SPLIT READS # #

# quality trim into longer assembly reads (minimum average quality 3, minimum length 1,000 bp) & remove symbiodiniacaea reads
gunzip -c ONT_reads_pc.fastq.gz | chopper -q 3 -l 1000 --contam $ref_genomes/all_symbiodiniaceae_v2.fna.gz | gzip > assembly_reads.fastq.gz

# quality trim into shorter polishing reads (minimum average quality 5, minimum length 500 bp) & remove symbiodiniacaea reads
gunzip -c ONT_reads_pc.fastq.gz | chopper -q 5 -l 500 --contam $ref_genomes/all_symbiodiniaceae_v2.fna.gz | gzip > polishing_reads.fastq.gz
