#!/bin/bash

#########################################################################################
# MAP NANOPORE PORECHOPPED READS TO MITOCHONDRIAL REFERENCE GENOME OF PORITES HARRISONI #
#########################################################################################

# Q10 quality of assembly reads
gunzip -c phar_assembly_reads.fastq.gz | chopper -q 10 | gzip > phar_assembly_reads_q10.fastq.gz

# Define input files
READS="phar_assembly_reads.fastq.gz"
REF="Phar_mitogenome_ref.fasta"
OUTPUT_PREFIX="mapped_mito_reads"

# Index the reference genome (if not already indexed)
minimap2 -d ${REF}.mmi ${REF}

# Align reads to the mitochondrial genome and output in BAM format
minimap2 -ax map-ont ${REF}.mmi ${READS} | samtools view -bS - | samtools sort -o ${OUTPUT_PREFIX}.sorted.bam

# Index the sorted BAM file
samtools index ${OUTPUT_PREFIX}.sorted.bam

# Extract mapped reads
samtools view -b -F 4 ${OUTPUT_PREFIX}.sorted.bam > ${OUTPUT_PREFIX}.mapped.bam

# Convert extracted reads back to FASTQ
samtools fastq ${OUTPUT_PREFIX}.mapped.bam > ${OUTPUT_PREFIX}.fastq

# Compress the output FASTQ file
gzip ${OUTPUT_PREFIX}.fastq

# THEN ASSEMBLE USING CANU
# CANU v2.3
canu -d CANU -p Phar genomeSize=18k -nanopore -trimmed mapped_mito_reads.fastq.gz gridOptions="--cpus=50"

# assembled mitochondrial genome will be Phar.contigs.fasta

##########################
# CIRCULARIZE MITOGENOME #
##########################

docker pull sangerpathogens/circlator
docker run -v /home/fiesingera/proj/pag_longread_genomes/Phar_Genome_2024/mitogenome:/data -it sangerpathogens/circlator bash

# inside docker container:
circlator all --assembler canu --merge_min_id 85 --merge_breaklen 1000 /data/Phar.contigs.trim.fasta /data/mapped_mito_reads.fastq.gz /data/circlator_out_CANU_trim_V3

gfastats 06.fixstart.fasta
# # scaffolds: 1
# Total scaffold length: 18645

# PUBLISHED MITOGENOME PHAR (Terraneo et al. 2018): 18,630bp

#############################
# # POLISHING USING RACON # #
#############################

mkdir POLISH
cd POLISH

scp circlator_out_CANU_trim_V3/06.fixstart.fasta PAG_Phar_mitogenome_v1.fasta

ln -s /home/fiesingera/proj/pag_longread_genomes/host_dna_ont/raw_data/phar_polishing_reads.fastq.gz .

# align polishing reads to assembly
minimap2 -a -x map-ont -t 32 PAG_Phar_mitogenome_v1.fasta phar_polishing_reads.fastq.gz -o racon.sam

# polish with alignment
/home/fiesingera/tools/racon/build/bin/racon -u -m 3 -x -5 -g -4 -w 500 -t 32 phar_polishing_reads.fastq.gz racon.sam PAG_Phar_mitogenome_v1.fasta > PAG_Phar_mitogenome_v2.polished.fasta

gfastats PAG_Phar_mitogenome.polished.fasta

# # scaffolds: 1
# Total scaffold length: 18640
