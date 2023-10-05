# CAGE-Seq-snakemake
The snakemake pipeline for CAGE-Seq analysis

The CAGE-Seq-snakemake is a bioinformatics analysis pipeline used for CAGE-seq sequencing data.

The pipeline is built using snakemake framework and takes raw fastq-files as input and includes steps for linker and artefact trimming (cutadapt), rRNA removal (SortMeRNA, alignment to a reference genome (BOWTIE2) and CAGE tag counting and clustering (paraclu). Additionally, several quality control steps (FastQC, RSeQC, MultiQC) are included to allow for easy verification of the results after a run.
