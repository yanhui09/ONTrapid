# Allow users to fix the underlying OS via singularity.
container: "docker://snakemake/snakemake:latest"

import pandas as pd
#-----------------
configfile: "config.yaml"

#-----------------
WORK_DIR = config["work_dir"]
OUT_DIR = WORK_DIR + "/results"

# load sample table
def load_sample_table(sample_table=WORK_DIR + "/samples.tsv"):
    sampleTable = pd.read_csv(sample_table, index_col=0, sep="\t")
    return sampleTable

sampleTable = load_sample_table()
SAMPLES = sampleTable.index.values

#---------------------------------
rule all:
    input: 
        expand(OUT_DIR + "/{sample}/{qc}_summary.tsv", 
        sample=SAMPLES,qc=["busco", "quast"])

# intialize input     
def get_input(wildcards):
    fqs_dir = sampleTable.loc[wildcards.sample, 'fq_dir']
    return expand(fqs_dir + "/{f}.fastq.gz", 
    f = glob_wildcards(fqs_dir + "/{f}.fastq.gz").f)

rule init:
    input: get_input
    output: OUT_DIR + "/init/{sample}.fastq"
    message: "Initialize input [{wildcards.sample}]"
    log: OUT_DIR + "/logs/init/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/init/{sample}.txt"
    shell:
        "zcat {input} > {output} 2> {log}"

include: "rules/qc.smk"
include: "rules/assembly.smk"
include: "rules/variants.smk"