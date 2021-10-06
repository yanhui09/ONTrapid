rule qc:
    input: OUT_DIR+"/init/{sample}.fastq",
    output: directory(OUT_DIR + "/{sample}/qc_raw"),
    message: "NanoQC for raw reads"
    conda: "../envs/nanoqc.yaml"
    log: OUT_DIR + "/logs/qc_raw/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/qc_raw/{sample}.txt"
    shell:
        "nanoQC -o {output} {input} > {log} 2>&1"

rule porechop:
    input:
        rules.qc.output,
        raw = rules.init.output,
    output: OUT_DIR + "/{sample}/porechopped.fastq",
    message: "Trimming adapters with Porechop"
    conda: "../envs/porechop.yaml"
    log:  OUT_DIR + "/logs/porechop/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/porechop/{sample}.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        porechop -i {input.raw} -o {output} \
        --threads {threads} > {log} 2>&1
        """

def use_porchop(x):
    if not x:
        return rules.init.output
    else:
        return rules.porechop.output

rule nanofilt:
    input: use_porchop(config["use_porchop"])
    output: OUT_DIR + "/{sample}/nanofilted.fastq",
    message: "Filter low-quality reads"
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
        NanoFilt -q {params.q} -l {params.l} \
        --headcrop {params.headcrop} \
        --tailcrop {params.tailcrop} \
        {input} > {output} 2> {log}
        """

use rule qc as post_qc with:
    input: 
        rules.nanofilt.output
    output: 
        directory(OUT_DIR + "/{sample}/qc_clean"),
    message: 
        "NanoQC for clean reads"
    log: 
        OUT_DIR + "/logs/qc_clean/{sample}.log"
    benchmark: 
        OUT_DIR + "/benchmarks/qc_clean/{sample}.txt"

rule qc_stat:
    input:
        rules.post_qc.output,
    	fastq = rules.nanofilt.output,
    output: OUT_DIR + "/{sample}/qc_stats.tsv",
    message: "QC statistics"
    conda: "../envs/nanoqc.yaml"
    log: OUT_DIR + "/logs/qc_stat/{sample}.log"
    benchmark: OUT_DIR +"/benchmarks/qc_stat/{sample}.txt"
    shell:
    	"NanoStat --fastq {input.fastq} --tsv > {output} 2> {log}"