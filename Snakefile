shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")


import csv
import os
import json

# configfile: "config.yaml"

localrules: all

FILES = json.load(open('samples.json'))

BWA_INDEX = '/home/yx157/data/Index/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa'  ## update later
SAMPLES = sorted(FILES.keys())

## list all samples by sample_name and sample_type
MARK_SAMPLES = []
for sample in SAMPLES:
    for sample_type in FILES[sample].keys():
        MARK_SAMPLES.append(sample + "_" + sample_type)

# which sample_type is used as control for calling peaks: e.g. Input, IgG...
# CONTROL = config["control"]
CONTROL = 'INPUT'
CONTROLS = [sample for sample in MARK_SAMPLES if CONTROL in sample] ## mark_sample is for the macs2 call
CASES = [sample for sample in MARK_SAMPLES if CONTROL not in sample]

## multiple samples may use the same control input/IgG files
CONTROLS_UNIQUE = list(set(CONTROLS))

## list BAM files
CONTROL_BAM = expand("03aln/{sample}.sorted.bam", sample=CONTROLS_UNIQUE) ## bam for each sample 
CASE_BAM = expand("03aln/{sample}.sorted.bam", sample=CASES)

## peaks and bigwigs
ALL_PEAKS = []
ALL_inputSubtract_BIGWIG = []

print( CASES)
for case in CASES:
    control = "_".join(case.split("_")[0:-1]) + "_" + CONTROL
    print (control) 
    print( CONTROLS)
    if control in CONTROLS:
        ALL_PEAKS.append("05peak_macs2/{}_vs_{}_macs2_peaks.xls".format(case, control))
        ALL_inputSubtract_BIGWIG.append("06bigwig_inputSubtract/{}_subtract_{}.bw".format(case, control))
        # ALL_SUPER.append("11superEnhancer/{}_vs_{}-super/".format(case, control))
ALL_SAMPLES = CASES + CONTROLS_UNIQUE
ALL_BAM     = CONTROL_BAM + CASE_BAM


# ALL_FIRST_QC  = expand("02fqc/{sample}_{read}.fastqc.zip", sample = ALL_SAMPLES, read = ["R1", "R2"])
ALL_TRIMMED_FASTQ_1 = expand("01_trim_seq/{sample}_1.fastq.gz", sample = ALL_SAMPLES)
ALL_TRIMMED_FASTQ_2 = expand("01_trim_seq/{sample}_2.fastq.gz", sample = ALL_SAMPLES)
ALL_FASTQC  = expand("02fqc/{sample}_1_fastqc.html", sample = ALL_SAMPLES)
ALL_INDEX = expand("03aln/{sample}.sorted.bam.bai", sample = ALL_SAMPLES)
ALL_BAM = CONTROL_BAM + CASE_BAM
ALL_FLAGSTAT = expand("03aln/{sample}.sorted.bam.flagstat", sample = ALL_SAMPLES)
ALL_BIGWIG = expand("04bigwig/{sample}.bw", sample = ALL_SAMPLES)
# ALL_QC = ["10multiQC/multiQC_log.html"]

TARGETS = []
TARGETS.extend(ALL_FASTQC)
TARGETS.extend(ALL_TRIMMED_FASTQ_1)
TARGETS.extend(ALL_TRIMMED_FASTQ_2)
TARGETS.extend(ALL_BAM)
TARGETS.extend(ALL_INDEX)
TARGETS.extend(ALL_PEAKS)
TARGETS.extend(ALL_BIGWIG)
TARGETS.extend(ALL_inputSubtract_BIGWIG)
TARGETS.extend(ALL_FLAGSTAT)
# TARGETS.extend(ALL_QC)

localrules: all
rule all:
    input: TARGETS


## get a list of fastq.gz files for the same mark, same sample
def get_fastq_r1(wildcards):

    s_id = "_".join(wildcards.sample.split("_")[0:-1])  ## extract the last _ element for the type
    mark = wildcards.sample.split("_")[-1]
    return FILES[s_id][mark]['R1']

def get_fastq_r2(wildcards):
    # print(wildcards)
    s_id = "_".join(wildcards.sample.split("_")[0:-1])  ## extract the last _ element for the type
    mark = wildcards.sample.split("_")[-1]
    return FILES[s_id][mark]['R2']

rule trim_adapter:
    input:  get_fastq_r1, get_fastq_r2 
    output: "01_trim_seq/{sample}_1.fastq.gz" , "01_trim_seq/{sample}_2.fastq.gz"
    params : jobname = "{sample}"
    threads : 4
    shell: 
        """
        NGmerge  -a  -n 4 -1 {input[0]} -2 {input[1]}  -o 01_trim_seq/{params.jobname}
        """

rule fastqc:
    input:  "01_trim_seq/{sample}_1.fastq.gz" , "01_trim_seq/{sample}_2.fastq.gz"
    output: "02fqc/{sample}_1_fastqc.html" , "02fqc/{sample}_2_fastqc.html"
    log:    "00_log/{sample}_fastqc"
    threads: 1
    params : jobname = "{sample}"
    message: "fastqc {input}: {threads}"
    shell:
        """
        module load fastqc
        fastqc -o 02fqc -f fastq --noextract {input}  2> {log}
        """
# get the duplicates marked sorted bam, remove unmapped reads by samtools view -F 4 and dupliated reads by samblaster -r
# samblaster should run before samtools sort

rule bwa_align:
    input:  "01_trim_seq/{sample}_1.fastq.gz" , "01_trim_seq/{sample}_2.fastq.gz"
    output: "03aln/{sample}.sorted.bam"
    threads: 6
    message: "aligning {input}: {threads} threads"
    log:
        markdup = "00log/{sample}.markdup"
    shell:
        """
        module load bwa
        module load samtools
        bwa mem  -t {threads} {BWA_INDEX} {input[0]} {input[1]}  \
        | samblaster --removeDups \
	    | samtools view -Sb -F 4 - \
	    | samtools sort -m 6G -@ {threads} -T {output[0]}.tmp -o {output[0]} 2> {log.markdup}
        """

rule index_bam:
    input:  "03aln/{sample}.sorted.bam"
    output: "03aln/{sample}.sorted.bam.bai"
    log:    "00log/{sample}.index_bam"
    threads: 1
    params: jobname = "{sample}"
    message: "index_bam {input}: {threads} threads"
    shell:
        """
        samtools index {input} 2> {log}
        """

# check number of reads mapped by samtools flagstat, the output will be used for downsampling
rule flagstat_bam:
    input:  "03aln/{sample}.sorted.bam"
    output: "03aln/{sample}.sorted.bam.flagstat"
    log:    "00log/{sample}.flagstat_bam"
    threads: 1
    params: jobname = "{sample}"
    message: "flagstat_bam {input}: {threads} threads"
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """

rule make_inputSubtract_bigwigs:
    input : "03aln/{control}.sorted.bam", "03aln/{case}.sorted.bam" #, "03aln/{control}.sorted.bam.bai", "03aln/{case}.sorted.bam.bai"
    output:  "06bigwig_inputSubtract/{case}_subtract_{control}.bw"
    log: "00log/{case}_{control}_inputSubtract.makebw"
    threads: 5
    params: jobname = "{case}"
    message: "making input subtracted bigwig for {input}"
    shell:
        """
        module load deepTools
        bamCompare --bamfile1 {input[1]} --bamfile2 {input[0]} \
        --normalizeUsing RPKM \
        --binSize 30 \
        --smoothLength 300 \
        -p {threads} \
         --extendReads 200 -o {output} 2> {log}
        """

rule make_bigwigs:
    input : "03aln/{sample}.sorted.bam", "03aln/{sample}.sorted.bam.bai"
    output: "04bigwig/{sample}.bw"
    log: "00log/{sample}.makebw"
    threads: 5
    params: jobname = "{sample}"
    message: "making bigwig for {input}"
    shell:
        """
        module load deepTools
        bamCoverage -b {input[0]} --normalizeUsing RPKM --binSize 30 --smoothLength 300 -p 5 --extendReads 200 -o {output} 2> {log}
        """


rule call_peaks_macs2:
    input: control = "03aln/{control}.sorted.bam", case="03aln/{case}.sorted.bam"
    output: bed = "05peak_macs2/{case}_vs_{control}_macs2_peaks.xls"  ## case control is defined by the output 
    log: "00log/{case}_vs_{control}_call_peaks_macs2.log"
    params:
        name = "{case}_vs_{control}_macs2",
        jobname = "{case}"
    message: "call_peaks macs2 {input}: {threads} threads"
    shell:
        """
       module load macs2
       ## for macs2, when nomodel is set, --extsize is default to 200bp, this is the same as 2 * shift-size in macs14.
        macs2 callpeak -t {input.case} \
            -c {input.control} --keep-dup all -f BAM -g hs \
            --outdir 05peak_macs2 -n {params.name} -p 1e-5 --broad --broad-cutoff 1e-5 --nomodel &> {log}
        """


rule multiQC:
    input :
        # expand("00log/{sample}.align", sample = ALL_SAMPLES), ## the sample
        expand("03aln/{sample}.sorted.bam.flagstat", sample = ALL_SAMPLES),
        expand("02fqc/{sample}_1_fastqc.html", sample = ALL_SAMPLES)
    output: "10multiQC/multiQC_log.html"
    log: "00log/multiqc.log"
    message: "multiqc for all logs"
    shell:
        """
        multiqc 02fqc 03aln 00log -o 10multiQC -d -f -v -n multiQC_log 2> {log}
        """

