rule flye:
	input: 
        rules.qc_stat.output,
        fastq = rules.nanofilt.output,
	output:
		fasta = OUT_DIR + "/{sample}/flye/assembly.fasta",
		expand(OUT_DIR + "/{{sample}}/flye/assembly_graph.{s}", s=["gfa","gv"]),
        OUT_DIR + "/{sample}/flye/assembly_info.txt",
	message:
		"Assembly with Flye"
	params:
		_dir=OUT_DIR + "/{sample}/flye",
        mode=config["flye"]["mode"]
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
	message:
		"Assembly stats with quast"
	conda: "../envs/quast.yaml"
    log: OUT_DIR + "/logs/quast/flye/{sample}.log"
	benchmark: OUT_DIR + "/benchmarks/quast/flye/{sample}.txt"
	threads: config["threads"]["normal"]
	shell:
		"""
		quast.py {input} -o {output} \
        --threads {threads} > {log} 2>&1
		"""

# canu
rule canu:
    input:
        rules.qc_stat.output,
        fastq = rules.nanofilt.output,
	output: 
        fasta = OUT_DIR + "/{sample}/canu/"),
	message:
		"Assembly stats with quast"
	conda: "../envs/quast.yaml"
    log: OUT_DIR + "/logs/quast/flye/{sample}.log"
	benchmark: OUT_DIR + "/benchmarks/quast/flye/{sample}.txt"
	threads: config["threads"]["normal"]
	shell:
		"""
		quast.py {input} -o {output} \
        --threads {threads} > {log} 2>&1
		"""

# consensus from multiple contig sets, quickmerge
# To do: include options for other assemblers, hybrid-assembly



# polish with racon and medaka
