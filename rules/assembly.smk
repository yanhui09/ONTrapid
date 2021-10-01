rule flye:
    input: 
        rules.qc_stat.output,
        fastq = rules.nanofilt.output,
    output:
        fasta = OUT_DIR + "/{sample}/flye/assembly.fasta",
        expand(OUT_DIR + "/{{sample}}/flye/assembly_graph.{s}", s=["gfa","gv"]),
        OUT_DIR + "/{sample}/flye/assembly_info.txt",
    message: "Assembly with Flye"
    params:
        _dir=OUT_DIR + "/{sample}/flye",
        mode=config["flye"]["mode"],
    conda: "../envs/flye.yaml"
    log: OUT_DIR + "/logs/flye/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/flye/{sample}.txt"
    threads: config["threads"]["large"]
    shell:
        """
        flye {params.mode} {input} --out-dir {params._dir} \
        --threads {threads} > {log} 2>&1
        """

# canu
rule canu:
    input:
        rules.qc_stat.output,
        fastq = rules.nanofilt.output,
    output:
        fasta = OUT_DIR + "/{sample}/canu/assembly.fasta"),
    message: "Assembly with canu"
    params:
	p="assembly",
	d=OUT_DIR + "/{sample}/canu",
	size=config["canu"]["size"],
	usegrid=conig["canu"]["usegrid"],
	grid_opts=config["canu"]["grid_opts"],
	ex_args=config["canu"]["args"],
    conda: "../envs/canu.yaml"
    log: OUT_DIR + "/logs/canu/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/canu/{sample}.txt"
    threads: config["threads"]["large"]
    shell:
	"canu -p {params.p} -d {params.d}/ -nanopore-raw {input}"
	" genomeSize={params.size} {params.ex_args}"
	" useGrid={params.usegrid} gridOptions={params.grid_opts}"
	" > {log} 2>&1"

# consensus from multiple contig sets, two-round quickmerge
# To do: include options for hybrid-assembly
rule quickmerge:
    input:
        query = rules.flye.output.fasta,
        ref = rules.canu.output.fasta,
    output: 
        fasta = OUT_DIR + "/{sample}/quickmerge1/assembly.fasta",
    message: "Merge assemblies"
    params: 
        ml=config["quickmerge"]["ml"],
        c=config["quickmerge"]["c"],
        hco=config["quickmerge"]["hco"],
        p=OUT_DIR + "/{sample}/quickmerge1/assembly",
    conda: "../envs/quickmerge.yaml"
    log: OUT_DIR + "/logs/quickmerge1/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/quickmerge1/{sample}.txt"
    shell:
        "merge_wrapper.py {input.query} {input.ref}"
        " -ml {params.ml} -c {params.c}"
        " -hco {params.hco} -p {params.p}"
        " > {log} 2>&1"

use rule quickmerge as quickmerge1 with:
    input:
        query = rules.canu.output.fasta,
        ref = rules.flye.output.fasta,
    output: 
        fasta = OUT_DIR + "/{sample}/quickmerge2/assembly.fasta",
    params: 
        p=OUT_DIR + "/{sample}/quickmerge2/assembly",
    log: OUT_DIR + "/logs/quickmerge2/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/quickmerge2/{sample}.txt"
 
rule quast:
    input: OUT_DIR + "/{sample}/{f}/assembly.fasta"
    output: directory(OUT_DIR + "/{sample}/quast/{f}"),
    message: "Assembly stats with quast"
    conda: "../envs/quast.yaml"
    log: OUT_DIR + "/logs/quast/{f}/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/quast/{f}/{sample}.txt"
    threads: config["threads"]["normal"]
    shell:
        "quast.py {input} -o {output} --threads {threads} > {log} 2>&1"


# polish with racon and medaka
# assuming quickmerge is better (flye > canu)
rule get_polish_input:
    input: 
        expand(OUT_DIR + "/{{sample}}/quast/{f}", 
        f=["flye", "canu", "quickmerge1", "quickmerge2"]),
        fasta = rules.quickmerge.output.fasta,
    output: OUT_DIR + "/{sample}/polish/raw.fasta"
    message: "In preparation for polishing"
    shell: "cp {input.fasta} {output}"

# align merged assemblies with raw reads
# reuse for racon iterations
rule minimap:
    input: 
      ref = OUT_DIR + "/{sample}/polish/{f}.fasta"
      fastq = rules.nanofilt.output,
    output: OUT_DIR + "/{sample}/polish/{f}.paf"
    message: "Alignments against merged assemblies"
    params:
        x=config["minimap"]["x"]
    conda: "../envs/polish.yaml"
    log: OUT_DIR + "/logs/minimap/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/minimap/{sample}.txt"
    threads: config["threads"]["large"]
    shell:
        "minimap2 -t {threads} -x {params.x}"
        " {input.ref} {input.fastq} > {output} 2> {log}"

def get_racon_input(wildcards):
    # adjust input based on racon iteritions
    if int(wildcards.iter) == 1:
        return(rules.minimap.output, rules.get_polish_input.output)
    else:
        prefix = OUT_DIR + "/{sample}/polish/racon_{iter}".format(sample = wildcards.sample, iter = str(int(wildcards.iter) - 1))
        return(prefix + ".paf", prefix + ".fasta")

rule racon:
    input:
        rules.nanofilt.output,
        get_racon_input,
    output: OUT_DIR + "/{sample}/polish/racon_{iter}.fasta"
    message: "Polish with racon, iterations={iter}"
    params:
        m=config["racon"]["m"],
        x=config["racon"]["x"],
        g=config["racon"]["g"],
        w=config["racon"]["w"],
    conda: "../envs/polish.yaml"
    log: OUT_DIR +"/logs/polish/racon/{sample}_round{iter}.log"
    benchmark: OUT_DIR +"/benchmarks//polish/racon/{sample}_round{iter}.txt"
    threads: config["threads"]["large"]
    shell:
        "racon -m {param.m} -x {params.x}"
        " -g {params.g} -w {params.w} -t {threads}"
        " {input} > {output} 2>{log}"

checkpoint medaka_consensus:
    input:
        fasta = expand(OUT_DIR + "/{{sample}}/polish/racon_{iter}.fasta", 
        iter = config["racon"]["iter"]),
        fastq = rules.nanofilt.output,
    output: 
        fasta = OUT_DIR + "/{sample}/polish/medaka/consensus.fasta"
    message: "Generate consensus with medaka"
    params:
        m=config["medaka"]["m"],
        _dir=OUT_DIR + "/{sample}/polish/medaka",
    conda: "../envs/polish.yaml"
    log: OUT_DIR + "/logs/polish/medaka/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/polish/medaka/{sample}.txt"
    threads: config["threads"]["large"]
    shell:
        "medaka_consensus -i {input.fastq}"
        " -d {input.fasta} -o {params._dir}"
        " -t {threads} -m {params.m} > {log} 2>&1"

# circularization with circulator

