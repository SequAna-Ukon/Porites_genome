# _Porites harrisoni_ reference genome assembly & annotation
Scripts used for the _Porites harrisoni_ referencegenome assembly and annotation using Oxford Nanopore longreads. The _Porites harrisoni_ reference genome is published under NCBI accession number: JBDLLT000000000 and the workflow is described in the following paper: xxx

## Assembly
* **00.reads_preprocessing.sh**

ONT reads were preprocessed using PoreChop (https://github.com/rrwick/Porechop) and quality control was done using FastQC (https://github.com/s-andrews/FastQC) and NanoPlot (https://github.com/wdecoster/NanoPlot). The reads were filtered and split into assembly reads (longer; min. length 1000 bp) & polishing reads (shorter; min. length 500 bp) using chopper (https://github.com/wdecoster/chopper).

* **01.kmer.sh**
  
Kmer profiling was done using meryl (https://github.com/marbl/meryl) and GenomeScope 2.0 (https://github.com/tbenavi1/genomescope2.0).

* **02.assembly_NECAT.sh**
  
The genome was assembled with NECAT (https://github.com/xiaochuanle/NECAT) using the longer assembly reads (see above) and the assembly was assessed using gfastats (https://github.com/vgl-hub/gfastats) and BUSCO (https://github.com/metashot/busco). The genome coverage was assessed by mapping the assembly to a Porites reference database from NCBI with bwa (https://github.com/lh3/bwa).

* **03.assembly_filtering.sh**
  
The assembly was filtered using BlobToolKit (https://blobtoolkit.genomehubs.org/).

* **04.assembly_polish.sh**
  
The assembled and filtered genome was polished using the shorter polishing reads with racon (https://github.com/lbcb-sci/racon) and Medaka (https://github.com/nanoporetech/medaka). The polished assembly was cleaned and headers were renamed and sorted using funannotate clean and funannotate sort (https://github.com/nextgenusfs/funannotate).

## Annotation
* **05.repeats.sh**
  
Repeats in the _Porites harrisoni_ genome were identified using EDTA (https://github.com/oushujun/EDTA) & RepeatModeler (https://github.com/Dfam-consortium/RepeatModeler) and soft-masked using RepeatMasker (https://github.com/rmhubley/RepeatMasker).

* **06.structural_annotation.sh**

Structural annotation was done using BRAKER3 (https://github.com/Gaius-Augustus/BRAKER). First, tRNAs were identified using tRNAscan-SE (https://github.com/UCSC-LoweLab/tRNAscan-SE), which were subsequently filtered for high-confidence tRNAs using EukHighConfidenceFilter implemented in tRNAscan-SE. Second, transcript evidence for the structural annotation was prepared from RNASeq data, which were trimmed using Trimmomatic (https://github.com/usadellab/Trimmomatic) and mapped to the _Porites harrisoni_ assembly using STAR (https://github.com/alexdobin/STAR). The resulting bam files were merged, and strand-specific RNA was extracted from the merged bam file. Using the script stringtie2utr.py (https://github.com/Gaius-Augustus/BRAKER/blob/utr_from_stringtie/scripts/stringtie2utr.py) from the BRAKER3 suite, untranslated regions (UTRs) were added to the gtf file output by BRAKER3. The high-confidence set of tRNAs was merged with the BRAKER3 structural predictions. Then, the gff3 file was checked for overlapping genes using AGAT (https://agat.readthedocs.io/en/latest/index.html) and validated using GenomeTools (https://github.com/genometools/genometools). Finally, GffRead (https://github.com/gpertea/gffread) was used to extract the predicted protein sequences from the merged file to use in the functional annotation.
  
* **07.functional_annotation.sh**

The predicted genes were annotated using InterProScan (https://github.com/ebi-pf-team/interproscan),  EggNOG-mapper (https://github.com/eggnogdb/eggnog-mapper) and Phobius (https://phobius.sbc.su.se/). The respective annotation files were then fed into funannotate annotate (https://github.com/nextgenusfs/funannotate) with the predicted genes in gff3 file format to synthesize all annotations. The final annotation was assessed using BUSCO (https://busco.ezlab.org/busco_userguide.html).




