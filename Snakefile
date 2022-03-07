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

include: "rules/qc.smk"
include: "rules/assembly.smk"
include: "rules/pangenomics.smk"
include: "rules/variants.smk"
#---------------------------------
rule all:
    input: 
        expand(OUT_DIR + "/{sample}/{qc}_summary.tsv", 
        sample=SAMPLES, qc=["busco", "quast"])

rule pangenome:
    input:
        rules.pan_genome.output.db,
        rules.calc_genome_similarity.output,
        rules.anvi_sum.output,
    output: touch(OUT_DIR + "/.pangenome_DONE")
