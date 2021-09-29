import pandas as pd
#-----------------
configfile: "config.yaml"

#-----------------
WORK_DIR = config["work_dir"]

# load sample table
def load_sample_table(sample_table=WORK_DIR + "/samples.tsv"):
    sampleTable = pd.read_csv(sample_table, index_col=0, sep="\t")
    return sampleTable

sampleTable = load_sample_table()

SAMPLES = sampleTable.index.values
def get_fqs_dir(sample, header):
    dirs = sampleTable.loc[sample, header]
    return list(dirs)

#---------------------------------
rull all:
     input: WORK_DIR + "/results/Finished"

# intialize input     
def get_input(wildcards):
    fqs_dir = get_fqs_dir(wildcards.sample, "fq_dir")
    return expand(fq_dir + "/{f}.fastq.gz", 
    f = glob_wildcards(fq_dir + "/{f}.fastq").f)

rule init:
    input: get_input
    output: WORK_DIR + "/results/init/{sample}.fastq"
    message: "Initialize input"
    log: WORK_DIR + "/results/logs/init/{sample}.log"
    benchmark: WORK_DIR + "/results/bechmarks/init/{sample}.txt"
    shell:
        "zcat {input} > {output} 2> {log}"

include: "rules/qc.smk"
include: "rules/assembly.smk"