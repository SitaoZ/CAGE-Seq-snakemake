configfile: "config.yaml"

WORKDIR = config["workdir"]
FILES = config["samples"]
SAMPLE = FILES.keys()
USED_GTF = config["GTF"]
USED_GFF = config["GFF"]


ADAPTER_FORWARD = config['adapter']['forward']
ADAPTER_BACKWARD = config['adapter']['backward']

CUTADAPT = config['software']['cutadapt']
BOWTIE2 = config['software']['bowtie2']
HISAT2 = config['software']['hisat2']
SAMTOOLS = config['software']['samtools']
RSCRIPT = config['software']['Rscript']

rule all:
    input:
        "10.dpi/outPooled/tc.bed.gz",
        "11.visualize/10.density/peak_transcript_density_comp.csv"

def get_input_fastqs(wildcards):
    sample_name = wildcards.sample
    return FILES[sample_name]

rule fastqc:
    input:
        "data/samples/{sample}.fastq.gz"
    output:
        out1="01.fastqc/{sample}_fastqc.html",
        out2="01.fastqc/{sample}_fastqc.zip"
    threads: 16
    resources:
        tmpdir="tmp"
    message: "fastqc {input}: {threads} threads"
    benchmark:
        "benchmarks/fastqc.{sample}.csv"
    log:
        stdout="logs/fastqc.{sample}.stdout",
        stderr="logs/fastqc.{sample}.stderr"
    shell:
        "fastqc -t {threads} -o 01.fastqc {input} > {log.stdout} 2>{log.stderr}"

rule cutadapt:
    input:
        qc_html="01.fastqc/{sample}_fastqc.html",
        qc_zip="01.fastqc/{sample}_fastqc.zip",
        read="data/samples/{sample}.fastq.gz"
    output:
        "02.cutadapt/{sample}.clean.fq.gz",
    threads: 16
    resources:
        tmpdir="tmp"
    message: "cutadapt {input.read}: {threads} threads"
    benchmark:
        "benchmarks/cutadapt.{sample}.csv"
    log:
        stdout="logs/cutadapt.{sample}.stdout",
        stderr="logs/cutadapt.{sample}.stderr"
    params:
        "-m 20 -q 20 --max-n 0.2 -e 0.08 -Z --trim-n"
    shell:
        "{CUTADAPT} {params} -j {threads} -a {ADAPTER_FORWARD} -o {output} {input.read}> {log.stdout} 2>{log.stderr}"

rule sortmerna:
    input:
        read="02.cutadapt/{sample}.clean.fq.gz"
    output:
        out="03.sortmerna/{sample}/{sample}_non_rRNA.fq.gz",
    threads: 16
    resources:
        tmpdir="tmp"
    message: "sortmerna {input.read} : {threads} threads"
    benchmark:
        "benchmarks/sortmerna.{sample}.csv"
    log:
        stdout="logs/sortmerna.{sample}.stdout",
        stderr="logs/sortmerna.{sample}.stderr"
    params:
        fixed = "--sam --num_alignments 1 --fastx",
        aligned = "03.sortmerna/{sample}/{sample}_rRNA",
        other = "03.sortmerna/{sample}/{sample}_non_rRNA",
        workdir = "03.sortmerna/{sample}",
        idx = "03.sortmerna/{sample}/idx",
        kvdb = "03.sortmerna/{sample}/kvdb",
        readb = "03.sortmerna/{sample}/readb"
    shell:
        "sortmerna {params.fixed} --threads {threads} --workdir {params.workdir} \
         --ref /home/zhusitao/software/sortmerna/data/rRNA_databases/silva-euk-18s-id95.fasta \
         --ref /home/zhusitao/software/sortmerna/data/rRNA_databases/silva-euk-28s-id98.fasta \
         --ref /home/zhusitao/software/sortmerna/data/rRNA_databases/rfam-5.8s-database-id98.fasta \
         --ref /home/zhusitao/software/sortmerna/data/rRNA_databases/rfam-5s-database-id98.fasta \
         --reads {input.read} --aligned {params.aligned} --other {params.other} -v > {log.stdout} 2>{log.stderr};"
         "rm -rf {params.idx} {params.kvdb} {params.readb}"


rule hisat2:
    input:
        read="03.sortmerna/{sample}/{sample}_non_rRNA.fq.gz",
    output:
        "04.hisat2/{sample}/{sample}.bam"
    threads: 16 
    resources:
        tmpdir="tmp"
    message: "hisat2 {input.read} : {threads} threads"
    benchmark:
        "benchmarks/hisat2.{sample}.csv"
    log:
        stdout="logs/hisat2.{sample}.stdout",
        stderr="logs/hisat2.{sample}.stderr"
    params:
        "--phred33 --sensitive"
    shell:
        "{HISAT2} {params} -p {threads} -x /home/zhusitao/database/plant/ath/tair10/hisat2_index/TAIR10 \
         -U {input.read} 2>{log.stderr} | {SAMTOOLS} view -b -o {output} -"

rule sort_bam:
    input:
        "04.hisat2/{sample}/{sample}.bam"
    output:
        "05.sortbam/{sample}/{sample}.sort.bam"
    threads: 16
    resources:
        tmpdir="tmp"
    message: "samtools sort {input} : {threads} threads"
    benchmark:
        "benchmarks/sort_bam.{sample}.csv"
    log:
        stdout="logs/sort_bam.{sample}.stdout",
        stderr="logs/sort_bam.{sample}.stderr"
    threads: 16
    shell:
        """
        export TMPDIR=05.sortbam/{wildcards.sample}/tmp
        mkdir -p ${{TMPDIR}}
        {SAMTOOLS} sort -@ {threads} -T ${{TMPDIR}} -o {output} {input}
        """

rule mark_dup:
    input:
        "05.sortbam/{sample}/{sample}.sort.bam"
    output:
        "06.markdup/{sample}/{sample}.rmdup.sort.bam"
    threads: 1
    resources:
        tmpdir="tmp"
    message: "picards {input} : {threads} threads"
    benchmark:
        "benchmarks/mark_dup.{sample}.csv"
    log:
        stdout="logs/markdup.{sample}.stdout",
        stderr="logs/markdup.{sample}.stderr"
    shell:
        "java -jar /home/zhusitao/software/picard/picard.jar MarkDuplicates  --REMOVE_DUPLICATES true  --INPUT {input}  --OUTPUT {output} --METRICS_FILE 06.markdup/{wildcards.sample}/{wildcards.sample}.sort.matrix > {log.stdout} 2> {log.stderr}"


rule ctss:
    input:
        "06.markdup/{sample}/{sample}.rmdup.sort.bam"
    output:
        "07.ctss/{sample}/ctssAll.gz"
    threads: 1
    message: "ctss {input} : {threads} threads"
    benchmark:
        "benchmarks/ctss.{sample}.csv"
    params:
        qual="20",
        tpm="5"
    shell:
        "sample_name=`ls {input} | sed 's/\.bam//g'`; output_name=`echo $sample_name | sed 's/06.markdup/07.ctss/g'`; sh scripts/make_ctss.sh -q {params.qual} -i $sample_name -n $output_name; sh scripts/process_ctss.sh -t {params.tpm} -o 07.ctss/{wildcards.sample} {WORKDIR}/07.ctss/{wildcards.sample}/{wildcards.sample}.rmdup.sort.ctss.bed; gzip 07.ctss/{wildcards.sample}/ctssAll"

rule rseqc:
    input:
        "06.markdup/{sample}/{sample}.rmdup.sort.bam"
    output:
        "08.rseqc/{sample}/tags_distribution.txt"
    threads: 1
    message: "rseqc {input} : {threads} threads"
    benchmark:
        "benchmarks/rseqc.{sample}.csv"
    log:
        stdout="logs/rseqc.{sample}.stdout",
        stderr="logs/rseqc.{sample}.stderr"
    shell:
        """read_distribution.py -i {input} -r /home/zhusitao/database/plant/ath/All_isoforms_GeneID/All_isoforms_V1.0.bed | awk 'NR>=5 && NR<=15' | awk '{{print $1","$2","$3","$4}}' > {output}"""


rule pre_dpi:
    input:
        rseqc_result="08.rseqc/{sample}/tags_distribution.txt",
        ctss_gz="07.ctss/{sample}/ctssAll.gz"
    output:
        "09.dpi_data/{sample}.ctss.gz"
    threads: 1
    message: "cp {input} : {threads} threads"
    benchmark:
        "benchmarks/pre_dpi.{sample}.csv"
    log:
        stdout="logs/dpi_data.{sample}.stdout",
        stderr="logs/dpi_data.{sample}.stderr"
    shell:
        "cp {input.ctss_gz} 09.dpi_data/{wildcards.sample}.ctss.gz > {log.stdout} 2>{log.stderr}"

rule dpi:
    input:
        expand("09.dpi_data/{sample}.ctss.gz", sample=SAMPLE)
    output:
        bed="10.dpi/outPooled/tc.bed.gz",
        permissive="10.dpi/outPooled/tc.spi_merged.ctssMaxCounts3.bed.gz",
        robust="10.dpi/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.bed.gz"
    threads: 1
    message: "DPI {input} : {threads} threads"
    benchmark:
        "benchmarks/dpi.csv"
    log:
        stdout="logs/dpi.stdout",
        stderr="logs/dpi.stderr"
    params:
        dir="10.dpi"
    shell:
        "sh {WORKDIR}/scripts/dpi1-master/identify_tss_peaks.sh -g {WORKDIR}/scripts/chrom.sizes -d N -i '{WORKDIR}/09.dpi_data/*ctss.gz' -d N -o {params.dir} > {log.stdout} 2> {log.stderr}"


rule visualize:
    input:
        bed="10.dpi/outPooled/tc.bed.gz",
        permissive="10.dpi/outPooled/tc.spi_merged.ctssMaxCounts3.bed.gz",
        robust="10.dpi/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.bed.gz"
    output:
        "11.visualize/10.density/peak_transcript_density_comp.csv"
    threads: 1
    message: "visualize {input} : {threads} threads"
    benchmark:
        "benchmarks/visualize.csv"
    log:
        stdout="logs/visualize.stdout",
        stderr="logs/visualize.stderr"
    shell:
        "sh scripts/stat.sh -r {WORKDIR}/{input.robust} -p {WORKDIR}/{input.permissive} -g {USED_GFF} -o {WORKDIR}/11.visualize"
