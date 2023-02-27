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
        ex_args="" if config["flye"]["ex_args"] is None else config["flye"]["ex_args"], 
    conda: "../envs/flye.yaml"
    log: OUT_DIR + "/logs/flye/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/flye/{sample}.txt"
    threads: config["threads"]["large"]
    shell:
        """
        flye {params.mode} {input.fastq} --out-dir {output._dir} \
        --threads {threads} {params.ex_args} > {log} 2>&1
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
        query = OUT_DIR + "/{sample}/%s/assembly.fasta" % config["quickmerge"]["query"],
        ref = OUT_DIR + "/{sample}/%s/assembly.fasta" % config["quickmerge"]["ref"],
    output: 
        fasta = OUT_DIR + "/{sample}/quickmerge/assembly.fasta",
        _dir = directory(OUT_DIR + "/{sample}/quickmerge"), 
    message: "Merge assemblies [{wildcards.sample}]"
    params:
        ml=config["quickmerge"]["ml"],
        c=config["quickmerge"]["c"],
        hco=config["quickmerge"]["hco"],
    conda: "../envs/quickmerge.yaml"
    log: OUT_DIR + "/logs/quickmerge/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/quickmerge/{sample}.txt"
    shell:
        # can't use absolute path in --prefix 
        """
        mkdir -p {output._dir}; cd {output._dir}
        merge_wrapper.py {input.query} {input.ref} \
        -ml {params.ml} -c {params.c} \
        -hco {params.hco} > {log} 2>&1
        mv merged_out.fasta {output.fasta}; cd - >> {log} 2>&1
        """

# circularization with circlator
# use polished assemblies and canu corrected|trimmed reads
# circlize genome and use dnaA as start if possible
rule circlator:
    input:
        fasta = OUT_DIR + "/{sample}/quickmerge2polish/assembly.fasta",
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

# polish with racon and medaka
rule get_polish_input:
    input: OUT_DIR + "/{sample}/{f}/assembly.fasta"
    output: OUT_DIR + "/{sample}/{f}2polish/raw.fasta"
    message: "In preparation for polishing {wildcards.f} assemblies [{wildcards.sample}]"
    shell: "cp {input} {output}"

# align merged assemblies with raw reads
# reused in racon iterations
rule minimap:
    input: 
      ref = OUT_DIR + "/{sample}/{f}2polish/{assembly}.fasta",
      fastq = rules.nanofilt.output,
    output: OUT_DIR + "/{sample}/{f}2polish/{assembly}.paf"
    message: "{wildcards.f}2polish: alignments against {wildcards.assembly} assembly [{wildcards.sample}]"
    params:
        x=config["minimap"]["x"]
    conda: "../envs/polish.yaml"
    log: OUT_DIR + "/logs/minimap/{sample}/{f}2polish/{assembly}.log"
    benchmark: OUT_DIR + "/benchmarks/minimap/{sample}/{f}2polish/{assembly}.txt"
    threads: config["threads"]["large"]
    shell:
        "minimap2 -t {threads} -x {params.x}"
        " {input.ref} {input.fastq} > {output} 2> {log}"

def get_racon_input(wildcards):
    # adjust input based on racon iteritions
    if int(wildcards.iter) == 1:
        prefix = OUT_DIR + "/{sample}/{f}2polish/raw"
        return(prefix + ".paf", prefix + ".fasta")
    else:
        prefix = OUT_DIR + "/{sample}/{f}2polish/racon_{iter}".format(sample=wildcards.sample, f=wildcards.f, iter=str(int(wildcards.iter) - 1))
        return(prefix + ".paf", prefix + ".fasta")

rule racon:
    input:
        rules.nanofilt.output,
        get_racon_input,
    output: OUT_DIR + "/{sample}/{f}2polish/racon_{iter}.fasta"
    message: "Polish {wildcards.f} assemblies with racon, round={wildcards.iter} [{wildcards.sample}]"
    params:
        m=config["racon"]["m"],
        x=config["racon"]["x"],
        g=config["racon"]["g"],
        w=config["racon"]["w"],
    conda: "../envs/polish.yaml"
    log: OUT_DIR +"/logs/polish/racon/{sample}/{f}/round{iter}.log"
    benchmark: OUT_DIR +"/benchmarks/polish/racon/{sample}/{f}/round{iter}.txt"
    threads: config["threads"]["large"]
    shell:
        "racon -m {params.m} -x {params.x}"
        " -g {params.g} -w {params.w} -t {threads}"
        " {input} > {output} 2> {log}"

rule medaka_consensus:
    input:
        fasta = expand(OUT_DIR + "/{{sample}}/{{f}}2polish/racon_{iter}.fasta", 
        iter = config["racon"]["iter"]),
        fastq = rules.nanofilt.output,
    output: 
        fasta = OUT_DIR + "/{sample}/{f}2polish/assembly.fasta",
    message: "Generate consensus for {wildcards.f} assemblies with medaka [{wildcards.sample}]"
    params:
        m=config["medaka"]["m"],
        _dir=OUT_DIR + "/{sample}/{f}2polish",
    conda: "../envs/polish.yaml"
    log: OUT_DIR + "/logs/polish/medaka/{sample}/{f}.log"
    benchmark: OUT_DIR + "/benchmarks/polish/medaka/{sample}/{f}.txt"
    threads: config["threads"]["large"]
    shell:
        """
        medaka_consensus -i {input.fastq} \
        -d {input.fasta} -o {params._dir}/medaka \
        -t {threads} -m {params.m} > {log} 2>&1;
        cp {params._dir}/medaka/consensus.fasta {output.fasta}
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
        mv {params.out_path}/{params.o}/*.txt {params.out_path}/{params.o}/summary.txt;
        cd - >> {log} 2>&1
        """

def choose_assembly(x):
    if x == 'canu' or x == 'flye':
        return x+'_out', x+'2polish_out'
    elif x == 'default':
        out = ("flye", "canu", "quickmerge", "circlator")
        out1 = [x + "_out" for x in out]
        out2 = [x + "2polish_out" for x in out]
        return  out1 + out2
    else:
        raise Exception('Assembler-opts only allows --only-canu, --only-flye, --default.\n{} is used in the config file'.format(x))

# summarize QC output
rule quast_summary:
    input: expand(OUT_DIR + "/{{sample}}/quast/{f}", f=choose_assembly(config["assembler_opts"]))
    output: OUT_DIR + "/{sample}/quast_summary.tsv",
    message: "QUAST summary [{wildcards.sample}]"
    run:
        import pandas as pd
        import os
        quast_stats = pd.DataFrame()
        for dir in input:
            d = pd.read_csv(dir + "/report.tsv", sep="\t", skiprows=1, index_col=0, header=None)
            quast_stats = quast_stats.append(d.T)
        
        quast_stats.index = [os.path.split(x)[-1] for x in input]
        quast_stats.T.to_csv(output[0], sep="\t")
 
rule busco_summary:
    input: expand(OUT_DIR + "/{{sample}}/busco/{f}", f=choose_assembly(config["assembler_opts"]))
    output: OUT_DIR + "/{sample}/busco_summary.tsv",
    message: "BUSCO summary [{wildcards.sample}]"
    run:
        import pandas as pd
        import os
        def busco2tdf(file):
            stats_lines = []
            with open(file) as f:
                for i, line in enumerate(f):
                    if 8 <= i <= 14:
                        stats_lines.append(line.lstrip())
            
            val = [x.split("\t")[0] for x in stats_lines]
            _index = ["Percentage"] + [x.split("\t")[1] for x in stats_lines[1:] if "\t" in x]
            df = pd.DataFrame(val, index=_index)
            return df

        busco_stats = pd.DataFrame()
        for dir in input:
            d = busco2tdf(dir + "/summary.txt")
            busco_stats = busco_stats.append(d.T)
        
        busco_stats.index = [os.path.split(x)[-1] for x in input]
        busco_stats.T.to_csv(output[0], sep="\t")

# collect polished assemblies
checkpoint collect_assembly:
    input: expand(OUT_DIR+ "/{sample}/flye2polish/assembly.fasta", sample = SAMPLES)
    output: directory(OUT_DIR + "/assembly")
    run:
        import os
        import shutil

        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for f in input:
            # extract sample name from path
            sample = os.path.split(f)[0].split("/")[-2]
            shutil.copy(f, output[0] + "/"  + sample + ".fasta")
