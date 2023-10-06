### CAGE-Seq-snakemake
The snakemake pipeline for CAGE-Seq analysis

The CAGE-Seq-snakemake is a bioinformatics analysis pipeline used for CAGE-seq sequencing data.

The pipeline is built using snakemake framework and takes raw fastq-files as input and includes steps for linker and artefact trimming (cutadapt), rRNA removal (SortMeRNA, alignment to a reference genome (HISAT2) and CAGE tag counting and clustering (BEDtools and DPI). Additionally, several quality control steps (FastQC, RSeQC) are included to allow for easy verification of the results after a run.

#### Software dependencies
1. Input read QC ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Filter and adapter trimming ([Cutadapt](https://github.com/marcelm/cutadapt))
3. rRNA removing ([SortMeRNA](https://github.com/sortmerna/sortmerna))
4. Reference alignment ([HISAT2](https://github.com/DaehwanKimLab/hisat2))
5. Sort bam ([SAMtools](http://www.htslib.org/))
6. Mark duplicates ([Picard](https://broadinstitute.github.io/picard/))
7. BAM to BED and CTSS ([BEDtools](https://github.com/arq5x/bedtools2))
8. Peak identification ([DPI](https://github.com/hkawaji/dpi1))
9. CAGE tag clustering QC ([RSeQC](https://rseqc.sourceforge.net/))

#### Step 1: Conda environment
```bash
conda create -n snakemake python=3.10.6
conda activate snakemake
conda install -c bioconda fastqc=0.12.1 cutadapt=3.1 sortmerna=4.3.4 hisat2=2.2.1 samtools=1.3.1 bedtools=v2.31.0 rseqc=5.0.1
```
#### Step 2: config.yaml
First, he config.yaml should be edited. If you don't know YAML ,the link is useful [yaml](https://www.cloudbees.com/blog/yaml-tutorial-everything-you-need-get-started).

```bash
cat config.yaml

workdir:
  "/data/zhusitao/project/04.learning/snakemake/03.CAGE-Seq"

adapter:
  forward: AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG
  backward: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA

software:
  fastqc: /home/zhusitao/anaconda3/envs/snakemake/bin/fastqc
  cutadapt: /home/zhusitao/anaconda3/bin/cutadapt
  sortmerna: /home/zhusitao/software/sortmerna-4.3.4-Linux/bin/sortmerna
  bowtie2: /home/zhusitao/anaconda3/bin/bowtie2
  hisat2: /home/zhusitao/anaconda3/bin/hisat2
  samtools: /home/zhusitao/anaconda3/bin/samtools
  Rscript: /home/zhusitao/anaconda3/envs/R/bin/Rscript

USED_GTF: /path/GTF/ath.gtf
```

#### Step 3: snakemake
Define an edited snakemake

#### Step 4: execute
```bash
nohup snakemake --snakefile Snakefile --configfile config.yaml --cores 50 --latency-wait 100 &
```
