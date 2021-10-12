rule rasusa:
    input: rules.init.output,
    output: OUT_DIR + "/{sample}/subsampled_%sx.fastq" % config["rasusa"]["coverage"], 
    message: "Subsampling at estimated depth {params.coverage}x [{wildcards.sample}]"
    params:
        coverage=config["rasusa"]["coverage"],
        genome_size=config["rasusa"]["genome_size"],
    conda: "../envs/rasusa.yaml"
    log: OUT_DIR + "/logs/rasusa/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/rasusa/{sample}.txt"
    shell:
        "rasusa --input {input} --coverage {params.coverage}"
        " --genome-size {params.genome_size} > {output}"
        " 2> {log}"

# decide subsampling first
def get_qc_input(x):
    if x is True:
        return rules.rasusa.output
    elif x is False:
        return rules.init.output
    else:
        raise Exception('Subsampling only allows bool type [TRUE/FALSE].\n{} is used in the config file'.format(x))

SUBS = config["rasusa"]["subsampling"]

rule qc:
    input: get_qc_input(SUBS)
    output: directory(OUT_DIR + "/{sample}/qc_raw"),
    message: "NanoQC for raw fastq [{wildcards.sample}]"
    conda: "../envs/nanoqc.yaml"
    log: OUT_DIR + "/logs/qc_raw/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/qc_raw/{sample}.txt"
    shell:
        "nanoQC -o {output} {input} > {log} 2>&1"

rule porechop:
    input:
        rules.qc.output,
        raw = get_qc_input(SUBS),
    output: OUT_DIR + "/{sample}/porechopped.fastq",
    message: "Trimming adapters with Porechop [{wildcards.sample}]"
    conda: "../envs/porechop.yaml"
    log:  OUT_DIR + "/logs/porechop/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/porechop/{sample}.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        porechop -i {input.raw} -o {output} \
        --threads {threads} > {log} 2>&1
        """

def use_porchop(x, y):
    if y is True:
        return rules.porechop.output
    elif y is False:
        if x is False:
            return rules.init.output
        else: # exception will be first raised at get_qc_input()
            return rules.rasusa.output
    else:
        raise Exception('Use_porechop only allows bool type [TRUE/FALSE].\n{} is used in the config file'.format(x))
        
rule nanofilt:
    input: use_porchop(SUBS, config["use_porchop"])
    output: OUT_DIR + "/{sample}/nanofilted.fastq",
    message: "Filter low-quality reads [{wildcards.sample}]"
    params:
        l=config["nanofilt"]["l"],
        q=config["nanofilt"]["q"],
        headcrop=config["nanofilt"]["headcrop"],
        tailcrop=config["nanofilt"]["tailcrop"],
    conda: "../envs/nanoqc.yaml"
    log: OUT_DIR + "/logs/nanofilt/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/nanofilt/{sample}.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        cat {input} | NanoFilt -q {params.q} -l {params.l} \
        --headcrop {params.headcrop} \
        --tailcrop {params.tailcrop} \
        > {output} 2> {log}
        """

use rule qc as post_qc with:
    input: 
        rules.nanofilt.output
    output: 
        directory(OUT_DIR + "/{sample}/qc_clean"),
    message: 
        "NanoQC for clean reads [{wildcards.sample}]"
    log: 
        OUT_DIR + "/logs/qc_clean/{sample}.log"
    benchmark: 
        OUT_DIR + "/benchmarks/qc_clean/{sample}.txt"

rule qc_stat:
    input:
        rules.post_qc.output,
    	fastq = rules.nanofilt.output,
    output: OUT_DIR + "/{sample}/qc_stats.tsv",
    message: "QC statistics [{wildcards.sample}]"
    conda: "../envs/nanoqc.yaml"
    log: OUT_DIR + "/logs/qc_stat/{sample}.log"
    benchmark: OUT_DIR +"/benchmarks/qc_stat/{sample}.txt"
    shell:
    	"NanoStat --fastq {input.fastq} --tsv > {output} 2> {log}"