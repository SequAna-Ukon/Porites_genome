# _Porites harrisoni_ reference genome assembly & annotation
Scripts used for the _Porites harrisoni_ referencegenome assembly and annotation using Oxford Nanopore longreads. The _Porites harrisoni_ reference genome is published under NCBI accession number: JBDLLT000000000 and the workflow is described in the following paper: xxx

## Assembly
* **00.reads_preprocessing.sh**

ONT reads were preprocessed using PoreChop (https://github.com/rrwick/Porechop) and quality control was done using FastQC (https://github.com/s-andrews/FastQC) and NanoPlot (https://github.com/wdecoster/NanoPlot). The reads were filtered and split into assembly (longer; min. length 1000 bp) & polishing (shorter; min. length 500 bp) reads using chopper (https://github.com/wdecoster/chopper).

* **01.kmer.sh**
  
Kmer profiling was done using meryl (https://github.com/marbl/meryl) and GenomeScope 2.0 (https://github.com/tbenavi1/genomescope2.0).

* **02.assembly_NECAT.sh**
  
The genome was assembled with NECAT (https://github.com/xiaochuanle/NECAT) using the longer assembly reads (see above) and the assembly was assessed using gfastats (https://github.com/vgl-hub/gfastats) and BUSCO (https://github.com/metashot/busco). The genome coverage was assessed by mapping the assembly to a Porites reference database from NCBI with bwa (https://github.com/lh3/bwa).

* **03.assembly_filtering.sh**
  
The assembly was filtered using BlobToolKit (https://blobtoolkit.genomehubs.org/).

* **04.assembly_polish.sh**
  
The assembled and filtered genome was polished using the shorter polishing reads with racon (https://github.com/lbcb-sci/racon) and Medaka (https://github.com/nanoporetech/medaka).

## Annotation
* **05.repeats.sh**
  
Repeats in the _Porites harrisoni_ genome were identified using EDTA (https://github.com/oushujun/EDTA) & RepeatModeler (https://github.com/Dfam-consortium/RepeatModeler) and soft-masked using RepeatMasker (https://github.com/rmhubley/RepeatMasker).

* **06.structural_annotation.sh**
### tRNA prediction

````bash
tRNAscan-SE -E -I -H --detail --thread 50 -o trnascan-se.out -f trnascan-se.tbl -m trnascan-se.log  PAG_UKon_Phar.fasta.masked

EukHighConfidenceFilter -i trnascan-se.out -s trnascan-se.tbl -o eukconf -p filt
````
- I used the latest Metazoa-specific protein set from [OrthoDB](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/).
- RNASeq reads are REVERSE-stranded according to Rach. Note then, the following commands have been used:
````bash
#mapping
STAR --runThreadN 50 --runMode genomeGenerate --genomeDir Phar_index --genomeFastaFiles PAG_UKon_Phar.fasta --genomeSAindexNbases 10

for i in `ls *1P.fq.gz|sed 's/_1P.fq.gz//g'`;do STAR --runThreadN 50 --genomeDir Phar_index --readFilesIn ${i}_1P.fq.gz ${i}_2P.fq.gz --readFilesCommand "zcat" --outSAMtype  BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFileNamePrefix ${i}_ --limitBAMsortRAM 10000000000; done

samtools merge -@ 50 Phar_RNASeqAll.STAR.bam *.sortedByCoord.out.bam

#Split Stranded RNA-Seq BAM 
#Plus strand
samtools view -h Phar_RNASeqAll.STAR.bam | awk 'BEGIN{OFS="\t"} /^@/ || ($2==99 || $2==83)' | samtools view -b -o Phar_plus_strand.bam

#Minus strand
samtools view -h Phar_RNASeqAll.STAR.bam | awk 'BEGIN{OFS="\t"} /^@/ || ($2==147 || $2==163)' | samtools view -b -o Phar_minus_strand.bam
````
- I used [BRAKER3](https://github.com/Gaius-Augustus/BRAKER)

````bash
sudo docker run --user 1000:100  -v $(pwd):/home/jovyan/work  teambraker/braker3:latest braker.pl --species=Porites_harrisoni --genome=work/PAG_UKon_Phar.clean.fasta.masked --bam=work/Phar_plus_strand.bam,work/Phar_minus_strand.bam --stranded=+,- --threads 50 --prot_seq=work/Metazoa.fa --busco_lineage=metazoa_odb10
````
- Decorate a braker.gtf with UTRs from a stringtie.gff file.

````bash
pip install intervaltree
python3.8 ../stringtie2utr.py -g braker.gtf -s GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff -o braker_with_utrs.gtf
````
- The script currently resides here: https://github.com/Gaius-Augustus/BRAKER/blob/utr_from_stringtie/scripts/stringtie2utr.py 

- Implement the tRNA prediction 

````bash
#covert tRNA to gff after removing non-high confident
perl convert_tRNAScanSE_to_gff3.pl --input=filter.out > trna_annotation.gff
#gtf to gff
cat braker/braker_with_utrs.gtf |gtf2gff.pl --gff3 -o braker.gff3

#merge gff files
agat_sp_merge_annotations.pl --gff braker.gff3 --gff trna_annotation.gff --out merged.gff
#export protein sequences to proceed with functional annotation
gffread merged.gff -g PAG_UKon_Phar.clean.fasta -y Phar.braker.prot.fasta
````

- As a new rule, i'm checking for overlapping genes using agat and validating the gff file using gt

````bash
agat_sp_fix_overlaping_genes.pl -f merged.gff -o Porites_harrisoni.gff3
gt gff3validator Porites_harrisoni.gff3
````
  

* **07.functional_annotation.sh**
tbd
