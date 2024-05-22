# _Porites harrisoni_ reference genome assembly & annotation
Scripts used for the _Porites harrisoni_ referencegenome assembly and annotation using Oxford Nanopore longreads. The _Porites harrisoni_ reference genome is published under NCBI accession number: JBDLLT000000000 and the workflow is described in the following paper: xxx

## Assembly
ONT reads were preprocessed using PoreChop (https://github.com/rrwick/Porechop) and quality control was done using FastQC (https://github.com/s-andrews/FastQC) and NanoPlot (https://github.com/wdecoster/NanoPlot). The reads were filtered and split into assembly (longer; min. length 1000 bp) & polishing (shorter; min. length 500 bp) reads using chopper (https://github.com/wdecoster/chopper)
* **00.reads_preprocessing.sh**
  
Kmer profiling was done using meryl (https://github.com/marbl/meryl) and GenomeScope 2.0 (https://github.com/tbenavi1/genomescope2.0)
* **01.kmer.sh**

The genome was assembled using NECAT (https://github.com/xiaochuanle/NECAT) and the longer assembly reads
* **02.assembly_NECAT.sh**

The assembly was filtered using BlobToolKit (https://blobtoolkit.genomehubs.org/)
* **03.assembly_filtering.sh**

The assembled and filtered genome was polished using the shorter polishing reads with racon (https://github.com/lbcb-sci/racon) and Medaka (https://github.com/nanoporetech/medaka)
* **04.assembly_polish.sh**

## Annotation
Repeats in the _Porites harrisoni_ genome were identified using EDTA (https://github.com/oushujun/EDTA) & RepeatModeler (https://github.com/Dfam-consortium/RepeatModeler) and soft-masked using RepeatMasker (https://github.com/rmhubley/RepeatMasker) 
* **05.repeats.sh**

Transcript and protein evidence were generated for the structural annotation using STAR (https://github.com/alexdobin/STAR), StringTie (https://github.com/gpertea/stringtie) and UniProt Database (https://www.uniprot.org/)
* **06.structural_annotation.sh**

Functional annotation was generated using funannotate (https://github.com/nextgenusfs/funannotate)
* **07.functional_annotation.sh**
