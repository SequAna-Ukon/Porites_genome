#!/bin/bash

conda activate repeats
conda install -c bioconda -c conda-forge barrnap
conda install -c bioconda vsearch

# ------------------------- REPEAT IDENTIFICATION & MASKING ---------------

assembly="/path/to/necat/PAG_Phar_UKon_1.1.fasta"

# # Identify rRNA # #
barrnap -q -k euk $assembly --threads 50 --outseq Phar_rrna.fasta > Phar_rrna.gff 

# # Pull docker images for RepeatModeler and EDTA # #
docker pull dfam/tetools:latest
docker pull quay.io/biocontainers/edta:2.2.2--hdfd78af_1

############## REPEAT MODELER ###############

docker run -v $PWD:/in -w /in dfam/tetools:latest BuildDatabase -name Phar_genome $assembly
docker run -v $PWD:/in -w /in dfam/tetools:latest RepeatModeler -database Phar_genome -LTRStruct -threads 40

# The results have been saved to:
#   /in/
#     Phar_genome-families.fa  - Consensus sequences for each family identified.
#     Phar_genome-families.stk - Seed alignments for each family identified.
#     Phar_genome-rmod.log     - Execution log.  Useful for reproducing results.

################## EDTA #####################
docker run -v $PWD:/in -w /in quay.io/biocontainers/edta:2.2.2--hdfd78af_1 EDTA.pl --genome $assembly --sensitive 1 --anno 1 -t 32 --overwrite 1 --force 1
   
# TE annotation using the EDTA library has finished! Check out:
#                 Whole-genome TE annotation (total TE: 43.97%): .EDTA.TEanno.gff3
#                 Whole-genome TE annotation summary: .EDTA.TEanno.sum
#                 Whole-genome TE divergence plot: _divergence_plot.pdf
#                 Whole-genome TE density plot: .EDTA.TEanno.density_plots.pdf

############ COMBINE REPEAT MODELER & EDTA ############### 

cat *-families.fa *.mod.EDTA.TElib.fa > Phar_RE_DB.fa

# #  Remove duplicates from repeats db # #
vsearch --derep_fulllength Phar_RE_DB.fa --output Phar_RE_DB_dedup.fasta

# # Get Repeat distribution # #
grep '>' Phar_RE_DB_dedup.fasta | sed -r 's/.+#//' | sed -r 's/\s+.+//' | sort | uniq -c

# # Repeats masking using Repeats DB from previous step # #
docker run -v $PWD:/in -w /in dfam/tetools:latest RepeatMasker $assembly -lib Phar_RE_DB_dedup.fasta -pa 8 -norna -xsmall

# The full report of the repeats percentage will be outputted after masking in PAG_Phar_UKon_1.1.fasta.tbl file

