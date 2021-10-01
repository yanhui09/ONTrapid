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

rule quast:
    input: rules.flye.output.fasta
    output: directory(OUT_DIR + "/{sample}/quast/flye"),
    message: "Assembly stats with quast"
    conda: "../envs/quast.yaml"
    log: OUT_DIR + "/logs/quast/flye/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/quast/flye/{sample}.txt"
    threads: config["threads"]["normal"]
    shell:
        "quast.py {input} -o {output} --threads {threads} > {log} 2>&1"

# canu
rule canu:
    input:
        rules.qc_stat.output,
        fastq = rules.nanofilt.output,
    output:
        fasta = OUT_DIR + "/{sample}/canu/{sample}.fasta"),
    message: "Assembly with canu"
    params:
	p="{sample}",
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

use rule quast as quast1 with:
    input: rules.canu.output.fasta
    output: directory(OUT_DIR + "/{sample}/quast/canu"),
    log: OUT_DIR + "/logs/quast/canu/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/quast/canu/{sample}.txt"

# consensus from multiple contig sets, two-round quickmerge
# To do: include options for hybrid-assembly
rule quickmerge:
    input:
        query = rules.flye.output.fasta,
        ref = rules.canu.output.fasta,
    output: 
        fasta = OUT_DIR + "/{sample}/quickmerge/1/merged.fasta",
    message: "Merge assemblies"
    params: 
        ml=config["quickmerge"]["ml"],
        c=config["quickmerge"]["c"],
        hco=config["quickmerge"]["hco"],
        p=OUT_DIR + "/{sample}/quickmerge/1/merged",
    conda: "../envs/quickmerge.yaml"
    log: OUT_DIR + "/logs/quickmerge/1/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/quickmerge/1/{sample}.txt"
    shell:
        "merge_wrapper.py {input.query} {input.ref}"
        " -ml {params.ml} -c {params.c}"
        " -hco {params.hco} -p {params.p}"
        " > {log} 2>&1"

use rule quast as quast2 with:
    input: rules.quickmerge.output.fasta
    output: directory(OUT_DIR + "/{sample}/quast/merge_1"),
    log: OUT_DIR + "/logs/quast/merge_1/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/quast/merge_1/{sample}.txt"

use rule quickmerge as quickmerge1 with:
    input:
        query = rules.canu.output.fasta,
        ref = rules.flye.output.fasta,
    output: 
        fasta = OUT_DIR + "/{sample}/quickmerge/2/merged.fasta",
    params: 
        p=OUT_DIR + "/{sample}/quickmerge/2/merged",
    log: OUT_DIR + "/logs/quickmerge/2/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/quickmerge/2/{sample}.txt"
 
use rule quast as quast3 with:
    input: rules.quickmerge1.output.fasta
    output: directory(OUT_DIR + "/{sample}/quast/merge_2"),
    log: OUT_DIR + "/logs/quast/merge_2/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/quast/merge_2/{sample}.txt"

# polish with racon and medaka
# assuming quickmerge is better (flye > canu)
# align merged assemblies with raw reads
# reuse for racon iterations
rule get_polish_input:
    input: 
        rules.quast.output,
        fasta = rules.quickmerge.output.fasta,

rule minimap:
    input: 
      ref = OUT_DIR + "/{sample}/polish/{f}.fas"
      fastq = rules.nanofilt.output,
    output: OUT_DIR + "/{sample}/quickmerge/1/merged.paf"
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
        return(rules.minimap.output, 
        rules..quickmerge.output.fasta)
    else:
        result = "{sample}/2.polish/racon/{sample}_racon_{iter}.fa".format(sample = wildcards.sample, iter = str(int(wildcards.iteration) - 1))
        return(result + '.paf', result)