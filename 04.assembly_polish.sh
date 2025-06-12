#!/bin/bash
# polish Phar_assembly_NECAT.fasta using Racon & Medaka

#############
# # RACON # #
#############

# directories
polishing_reads="/path/to/raw_reads/phar_polishing_reads.fastq.gz"
phar_assembly="/path/to/necat/Phar_assembly_NECAT.fasta"

# align polishing reads to assembly
minimap2 -a -x map-ont -t 32 $phar_assembly $polishing_reads -o phar_aligned_racon.sam

# polish with alignment
/path/to/racon -u -m 3 -x -5 -g -4 -w 500 -t 32 $polishing_reads phar_aligned_racon.sam $phar_assembly > Phar_assembly_NECAT_racon.fasta

##############
# # MEDAKA # #
##############

conda activate polish
conda install -c conda-forge -c bioconda medaka

# align polishing reads to assembly
polishing_reads="/path/to/raw_reads/phar_polishing_reads.fastq.gz"
phar_assembly="/path/to/necat/Phar_assembly_NECAT_racon.fasta"

minimap2 -a -x map-ont -t 32 $phar_assembly $polishing_reads -o phar_aligned_medaka.sam

# prepare files from sam to bam
samtools view -Sb -@ 32 phar_aligned_medaka.sam > phar_aligned_medaka.bam
samtools sort -@ 32 -O BAM phar_aligned_medaka.bam > phar_aligned_medaka.sort.bam
samtools index -b -@ 32 phar_aligned_medaka.sort.bam

# MEDAKA
# use model: r1041_e82_400bps_hac_v4.3.0 - high accuracy basecaller (hac) and latest R10.4.1 technology
medaka consensus phar_aligned_medaka.sort.bam phar_medaka.hdf --model  r1041_e82_400bps_hac_v4.3.0 --threads 48
medaka stitch --threads 48 phar_medaka.hdf $phar_assembly Phar_assembly_NECAT_racon_medaka.fasta

mv Phar_assembly_NECAT_racon_medaka.fasta Phar_assembly_NECAT_polished.fasta

######################
# # CLEAN ASSEMBLY # #
######################

conda activate anno_env

funannotate clean -i Phar_assembly_NECAT_polished.fasta -m 200 -o Phar_clean.fasta
funannotate sort -i Phar_clean.fasta -o PAG_Phar_UKon_1.1.fasta -b Phar_scaff

###########################################
# final assembly: PAG_Phar_UKon_1.1.fasta #
###########################################

#############
# # BUSCO # #
#############

busco -i PAG_Phar_UKon_1.1.fasta -m geno -l eukaryota_odb10 -c 30 -o busco_euk
busco -i PAG_Phar_UKon_1.1.fasta -m geno -l metazoa_odb10 -c 30 -o busco_meta
