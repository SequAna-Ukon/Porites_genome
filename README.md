# Porites_genome
Porites Genome assembly and annotation using ONT reads

## Adapter Removal using PoreChop:

````bash
porechop –-input reads.fastq -o reads_porechopped.fastq --discard_middle
````
## Quality Control:

Assess the quality of the raw nanopore reads using tools like FastQC or NanoPlot to identify any issues that may need attention.

https://github.com/wdecoster/NanoPlot

````bash
NanoPlot -t 2 --fastq reads1.fastq.gz reads2.fastq.gz --maxlength 40000 --plots dot --legacy hex
````
- Not very informative but we can give it a try

## Genome assembly
using Flye, CANU, and MaSuRCA .....

### 1. Flye

````bash
flye --nano-hq -g 0.421g --input [input.fastq] --out-dir [output_directory] --scaffold -t 50
````
### 2. CANU

````bash
canu -p [output_prefix] -d [output_directory] genomeSize=0.421g stopOnLowCoverage=5 -nanopore-raw [input.fastq]
````
### 3.MaSuRCA

````bash
runCanu.sh nanopore-[read_type] [config_file]
````
- Configuration file for MaSuRCA

````bash
# Configuration file for MaSuRCA

DATA
  PE = 
  JUMP = 
  OTHER = nanopore-[read_type] raw_reads.fastq
  # Add additional libraries as needed

PARAMETERS
  # Specify assembly parameters here, such as genome size estimate, k-mer size, etc.

````
## Kmer profiling
Usually, this is for short-reads or high-accurate long reads as "HiFi technology" but we could try
