rule flye:
    input: 
        rules.qc_stat.output,
        fastq = rules.nanofilt.output,
    output:
        fasta = OUT_DIR + "/{sample}/flye/assembly.fasta",
        _dir = directory(OUT_DIR + "/{sample}/flye"), # clean dir if failed
    message: "Assembly with Flye [{wildcards.sample}]"
    params:
        mode=config["flye"]["mode"],
    conda: "../envs/flye.yaml"
    log: OUT_DIR + "/logs/flye/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/flye/{sample}.txt"
    threads: config["threads"]["large"]
    shell:
        """
        flye {params.mode} {input.fastq} --out-dir {output._dir} \
        --threads {threads} > {log} 2>&1
        """

# canu
rule canu:
    input:
        rules.qc_stat.output,
        fastq = rules.nanofilt.output,
    output:
        fasta = OUT_DIR + "/{sample}/canu/assembly.fasta",
        fastq_cor = OUT_DIR + "/{sample}/canu/assembly.trimmedReads.fasta.gz",
        _dir = directory(OUT_DIR + "/{sample}/canu"), # clean dir if failed
    message: "Assembly with Canu [{wildcards.sample}]"
    params:
        size=config["canu"]["size"],
        usegrid=config["canu"]["usegrid"],
        grid_opts=config["canu"]["grid_opts"],
        stopOnLowCoverage=config["canu"]["stopOnLowCoverage"],
        minInputCoverage=config["canu"]["minInputCoverage"],
        ex_args=config["canu"]["ex_args"],
    conda: "../envs/canu.yaml"
    log: OUT_DIR + "/logs/canu/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/canu/{sample}.txt"
    threads: config["threads"]["large"]
    resources: 
      mem=config["memory"]["large_gb"],
    shell:
       """
       canu -p assembly -d {output._dir}/ -nanopore {input} \
       genomeSize={params.size} {params.ex_args} \
       minThreads={threads} maxThreads={threads} \
       maxMemory={resources.mem} \
       useGrid={params.usegrid} gridOptions={params.grid_opts} \
       stopOnLowCoverage={params.stopOnLowCoverage} \
       minInputCoverage={params.minInputCoverage} \
       > {log} 2>&1
       mv {output._dir}/assembly.contigs.fasta {output.fasta}
       """

# consensus from multiple contig sets, two-round quickmerge
# To do: include options for hybrid-assembly
rule quickmerge:
    input:
        query = rules.flye.output.fasta,
        ref = rules.canu.output.fasta,
    output: 
        fasta = OUT_DIR + "/{sample}/quickmerge1/assembly.fasta",
        _dir = directory(OUT_DIR + "/{sample}/quickmerge1"), 
    message: "Merge assemblies [{wildcards.sample}]"
    params:
        ml=config["quickmerge"]["ml"],
        c=config["quickmerge"]["c"],
        hco=config["quickmerge"]["hco"],
    conda: "../envs/quickmerge.yaml"
    log: OUT_DIR + "/logs/quickmerge1/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/quickmerge1/{sample}.txt"
    shell:
        # can't use absolute path in --prefix 
        """
        mkdir -p {output._dir}; cd {output._dir}
        merge_wrapper.py {input.query} {input.ref} \
        -ml {params.ml} -c {params.c} \
        -hco {params.hco} > {log} 2>&1
        mv merged_out.fasta {output.fasta}; cd - >> {log} 2>&1
        """

use rule quickmerge as quickmerge1 with:
    input:
        query = rules.canu.output.fasta,
        ref = rules.flye.output.fasta,
    output: 
        fasta = OUT_DIR + "/{sample}/quickmerge2/assembly.fasta",
        _dir = directory(OUT_DIR + "/{sample}/quickmerge2"), 
    params: # issues in snakemake, include all if changed
        ml=config["quickmerge"]["ml"],
        c=config["quickmerge"]["c"],
        hco=config["quickmerge"]["hco"],
    log: 
        OUT_DIR + "/logs/quickmerge2/{sample}.log"
    benchmark: 
        OUT_DIR + "/benchmarks/quickmerge2/{sample}.txt"
 
# circularization with circlator
# use polished assemblies and canu corrected|trimmed reads
# circlize genome and use dnaA as start if possible
rule circlator:
    input:
        fasta = rules.quickmerge.output.fasta,
        fastq_cor = rules.canu.output.fastq_cor,
    output:
        fasta = OUT_DIR + "/{sample}/circlator/assembly.fasta",
        _dir = directory(OUT_DIR + "/{sample}/circlator"),
    message: "Circularization with circlator [{wildcards.sample}]"
    params:
        assembler=config["circlator"]["assembler"],
        data_type=config["circlator"]["data_type"],
        bwa_opts=config["circlator"]["bwa_opts"],
        merge_min_id=config["circlator"]["merge_min_id"],
        merge_breaklen=config["circlator"]["merge_breaklen"],
    conda: "../envs/circlator.yaml"
    log: OUT_DIR + "/logs/circlator/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/circlator/{sample}.txt"
    threads: config["threads"]["large"]
    shell:
        # snakemake directory conflicts with os.mkdir() in circlator? existed already?
        # use .tmp to differentiate;
        """
        circlator all \
        --assembler {params.assembler} --data_type {params.data_type} \
        --bwa_opts {params.bwa_opts} --merge_min_id {params.merge_min_id} \
        --merge_breaklen {params.merge_breaklen} --threads {threads} \
        {input.fasta} {input.fastq_cor} {output._dir}.tmp > {log} 2>&1
        cp {output._dir}.tmp/06.fixstart.fasta {output.fasta}
        mv {output._dir}.tmp {output._dir}
        """

rule quast:
    input: OUT_DIR + "/{sample}/{f}/assembly.fasta"
    output: directory(OUT_DIR + "/{sample}/quast/{f}_out"), # avoid ambugious rules with assembly
    message: "{wildcards.f} assembly stats by QUAST [{wildcards.sample}]"
    conda: "../envs/quast.yaml"
    log: OUT_DIR + "/logs/quast/{f}/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/quast/{f}/{sample}.txt"
    threads: config["threads"]["normal"]
    shell:
        "quast.py {input} -o {output} --threads {threads} > {log} 2>&1"

# BUSCO assessment
rule busco:
    input: OUT_DIR + "/{sample}/{f}/assembly.fasta"
    output: directory(OUT_DIR + "/{sample}/busco/{f}_out"), # avoid ambugious rules with assembly
    message: "{wildcards.f} assembly stats by BUSCO [{wildcards.sample}]"
    params:
        o="{f}_out",
        out_path=OUT_DIR + "/{sample}/busco",
        l=config["busco"]["l"],
        m=config["busco"]["m"],
        opts=config["busco"]['opts'],
    conda: "../envs/busco.yaml"
    log: OUT_DIR + "/logs/busco/{f}/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/busco/{f}/{sample}.txt"
    threads: config["threads"]["normal"]
    shell:
        # cd busco_dir for downloads
        """
        cd {params.out_path};
        busco -i {input} -o {params.o} --out_path {params.out_path} \
        -l {params.l} -m {params.m} {params.opts} -c {threads} \
        > {log} 2>&1;
        cd - >> {log} 2>&1
        """

def fasta_to_polish(x):
    if x == '--only-canu':
        return (rules.canu.output.fasta, ["canu_out"])
    elif x == '--only-flye':
        return (rules.flye.output.fasta, ["flye_out"])
    elif x == '--default':
        return (rules.circlator.output.fasta, ["flye_out", "canu_out", "quickmerge1_out", "quickmerge2_out", "circlator_out"])
    else:
        raise Exception('Assembler-opts only allows --only-canu, --only-flye, --default.\n{} is used in the config file'.format(x))

# polish with racon and medaka
# assuming quickmerge is better (flye > canu)
rule get_polish_input:
    input: 
        fasta = fasta_to_polish(config["assembler_opts"])[0],
    output: OUT_DIR + "/{sample}/polish/raw.fasta"
    message: "In preparation for polishing [{wildcards.sample}]"
    shell: "cp {input.fasta} {output}"

# align merged assemblies with raw reads
# reuse for racon iterations
rule minimap:
    input: 
      ref = OUT_DIR + "/{sample}/polish/{f}.fasta",
      fastq = rules.nanofilt.output,
    output: OUT_DIR + "/{sample}/polish/{f}.paf"
    message: "Alignments against {wildcards.f} assembly [{wildcards.sample}]"
    params:
        x=config["minimap"]["x"]
    conda: "../envs/polish.yaml"
    log: OUT_DIR + "/logs/minimap/{sample}/{f}.log"
    benchmark: OUT_DIR + "/benchmarks/minimap/{sample}/{f}.txt"
    threads: config["threads"]["large"]
    shell:
        "minimap2 -t {threads} -x {params.x}"
        " {input.ref} {input.fastq} > {output} 2> {log}"

def get_racon_input(wildcards):
    # adjust input based on racon iteritions
    if int(wildcards.iter) == 1:
        prefix = OUT_DIR + "/{sample}/polish/raw"
        return(prefix + ".paf", prefix + ".fasta")
    else:
        prefix = OUT_DIR + "/{sample}/polish/racon_{iter}".format(sample = wildcards.sample, iter = str(int(wildcards.iter) - 1))
        return(prefix + ".paf", prefix + ".fasta")

rule racon:
    input:
        rules.nanofilt.output,
        get_racon_input,
    output: OUT_DIR + "/{sample}/polish/racon_{iter}.fasta"
    message: "Polish with racon, round={wildcards.iter} [{wildcards.sample}]"
    params:
        m=config["racon"]["m"],
        x=config["racon"]["x"],
        g=config["racon"]["g"],
        w=config["racon"]["w"],
    conda: "../envs/polish.yaml"
    log: OUT_DIR +"/logs/polish/racon/{sample}_round{iter}.log"
    benchmark: OUT_DIR +"/benchmarks/polish/racon/{sample}_round{iter}.txt"
    threads: config["threads"]["large"]
    shell:
        "racon -m {params.m} -x {params.x}"
        " -g {params.g} -w {params.w} -t {threads}"
        " {input} > {output} 2> {log}"

def choose_assembly_qc(x,y):
    if y == 'busco':
        return expand(OUT_DIR + "/{{sample}}/busco/{f}", 
        f=fasta_to_polish(x)[1])
    elif y == 'quast':
        return expand(OUT_DIR + "/{{sample}}/quast/{f}", 
        f=fasta_to_polish(x)[1])
    elif y == 'both':
        return expand(OUT_DIR + "/{{sample}}/{qc}/{f}", 
        f=fasta_to_polish(x)[1], qc=("busco", "quast"))
    else:
        raise Exception('Assembly_qc only allows busco, quast, both.\n{} is used in the config file'.format(y))
 
checkpoint medaka_consensus:
    input:
        choose_assembly_qc(config["assembler_opts"], config["assembly_qc"]),
        fasta = expand(OUT_DIR + "/{{sample}}/polish/racon_{iter}.fasta", 
        iter = config["racon"]["iter"]),
        fastq = rules.nanofilt.output,
    output: 
        fasta = OUT_DIR + "/{sample}/polish/assembly.fasta",
        tag = OUT_DIR + "/{sample}/Assembly.end",
    message: "Generate consensus with medaka [{wildcards.sample}]"
    params:
        m=config["medaka"]["m"],
        _dir=OUT_DIR + "/{sample}/polish",
    conda: "../envs/polish.yaml"
    log: OUT_DIR + "/logs/polish/medaka/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/polish/medaka/{sample}.txt"
    threads: config["threads"]["large"]
    shell:
        """
        medaka_consensus -i {input.fastq} \
        -d {input.fasta} -o {params._dir}/medaka \
        -t {threads} -m {params.m} > {log} 2>&1;
        cp {params._dir}/medaka/consensus.fasta {output.fasta};
        touch {output.tag}
        """
