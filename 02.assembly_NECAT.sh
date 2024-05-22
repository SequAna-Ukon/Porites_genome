#!/bin/bash
# Script for genome assembly using NECAT

# --------------------------------------
# # ENVIRONMENT & INSTALLATION # #
mamba activate assembly_env
necat_path="/PATH/TO/NECAT_INSTALLATION/"

# install from source code
git clone https://github.com/xiaochuanle/NECAT.git
cd NECAT/src/
make
cd ../Linux-amd64/bin
export PATH=$PATH:$(pwd)

# create config file
$necat_path/Linux-amd64/bin/necat.pl config phar_config.txt

# make a read list with where the reads to be processed are stored:
nano phar_read_list.txt

# edit the config file:
nano phar_config.txt

# PROJECT=phar_necat
# ONT_READ_LIST=/PATH/TO/phar_read_list.txt
# GENOME_SIZE=599000000
# THREADS=10
# MIN_READ_LENGTH=1000
# PREP_OUTPUT_COVERAGE=40
# OVLP_FAST_OPTIONS=-n 500 -z 20 -b 2000 -e 0.5 -j 0 -u 1 -a 1000
# OVLP_SENSITIVE_OPTIONS=-n 500 -z 10 -e 0.5 -j 0 -u 1 -a 1000
# CNS_FAST_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
# CNS_SENSITIVE_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
# TRIM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 1 -a 400
# ASM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400
# NUM_ITER=2
# CNS_OUTPUT_COVERAGE=30
# CLEANUP=1
# USE_GRID=false
# GRID_NODE=0
# GRID_OPTIONS=
# SMALL_MEMORY=0
# FSA_OL_FILTER_OPTIONS=
# FSA_ASSEMBLE_OPTIONS=
# FSA_CTG_BRIDGE_OPTIONS=
# POLISH_CONTIGS=true

# -----------------------------------------
# ASSEMBLY # #

# 1. corrections step:
$necat_path/Linux-amd64/bin/necat.pl correct phar_config.txt

# 2. assembly step:
$necat_path/Linux-amd64/bin/necat.pl assemble phar_config.txt

# 3. briding step: bridge the contigs (basically a scaffolding step)
$necat_path/Linux-amd64/bin/necat.pl bridge phar_config.txt

# rename assembly file
mv /PATH/TO/polished_contigs.fasta /PATH/TO/Phar_assembly_NECAT.fasta

# ------------------------------------------
# # QC OF THE ASSMEBLY # #
mamba activate qc_env

mamba install -c bioconda gfastats
mamba install -c bioconda busco

gfastats /PATH/TO/Phar_assembly_NECAT.fasta

busco -i /PATH/TO/Phar_assembly_NECAT.fasta -m geno -l eukaryota_odb10 -c 30 -o busco_euk

busco -i /PATH/TO/Phar_assembly_NECAT.fasta -m geno -l metazoa_odb10 -c 30 -o busco_meta

# ------------------------------------------
# # GENOME COVERAGE # #
mamba install -c bioconda bwa
mamba install -c bioconda samtools
mamba install -c conda-forge ncbi-datasets-cli

# PORITES REFERENCE FROM NCBI

# Download whole taxon reference genomes record and collect in single fasta file -- from NCBI database
# Example for "Porites" - Taxonomy ID: 46719

# # Download all the reference genome associated with the 46719 taxon
datasets download genome taxon 46719

# # Unzip the files
unzip ncbi_dataset.zip

# # Join together in single fasta file ready for indexing
cat ncbi_dataset/data/*/*.fna > ref_porites.fna

# Index the reference with bwa index
bwa index -a bwtsw /PATH/TO/ref_porites.fasta

# MAP NECAT ASSEMBLY TO PORITES REFERENCE
bwa mem -t 32 /PATH/TO/ref_porites.fasta /PATH/TO/Phar_assembly_NECAT.fasta | samtools view -S -b > phar.necat.against.porites.bam
samtools flagstats phar.necat.against.porites.bam
