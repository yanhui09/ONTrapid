REF=config["variants"]["ref"]

rule medaka_haploid_variants:
    input:
        fasta = REF, 
        fastq = rules.porechop.output,
    output: 
        vcf_anno = OUT_DIR + "/{sample}/variants/medaka/medaka.annotated.vcf",
        vcf = OUT_DIR + "/vcfs/{sample}.vcf",
    message: "Call variants with medaka [{wildcards.sample}]"
    params:
        m=config["variants"]["medaka"]["model"],
        _dir=OUT_DIR + "/{sample}/variants/medaka",
    conda: "../envs/polish.yaml"
    log: OUT_DIR + "/logs/varaints/medaka/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/variants/medaka/{sample}.txt"
    threads: config["threads"]["large"]
    shell:
        """
        medaka_haploid_variant -i {input.fastq} \
        -r {input.fasta} -o {params._dir}  \
        -m {params.m} -s -t {threads} \
        > {log} 2>&1
        cp {output.vcf_anno} {output.vcf}
        rm {input.fasta}.mmi {input.fasta}.fai
        """

# Clair3 https://github.com/HKU-BAL/Clair3
# Pepper https://github.com/kishwarshafin/pepper # Designed for diploid, no conda 
# Magalodon https://github.com/nanoporetech/megalodon # FASTA5 required
# Nanopolish https://github.com/jts/nanopolish # FASTA5 signals required
# Nanocaller https://github.com/WGLab/NanoCaller # Haploid mode not supported yet

rule minimap_bams:
    input:
        fasta = REF, 
        fastq = rules.porechop.output,
    output:
        ref = OUT_DIR + "/{sample}/variants/ref.fasta",
        bam = OUT_DIR + "/{sample}/variants/ref.sorted.bam",
        bai = OUT_DIR + "/{sample}/variants/ref.sorted.bam.bai",
    message: "Generate bam file with minimap2 [{wildcards.sample}]"
    params:
        _dir=OUT_DIR + "/{sample}/variants/clair",
        x=config["minimap"]["x"],
    conda: "../envs/polish.yaml"
    log: OUT_DIR + "/logs/varaints/minimap_bams/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/variants/minimap_bams/{sample}.txt"
    threads: config["threads"]["large"]
    shell:
        """
        [ ! -f {output.ref} ] && cp {input.fasta} {output.ref}
        minimap2 -L -t {threads} -ax {params.x} {output.ref} {input.fastq} 2> {log} | \
        samtools sort -o {output.bam} 2>> {log}
        samtools index {output.bam} 2>> {log}
        """ 

rule fadix:
    input: rules.minimap_bams.output.ref 
    output: OUT_DIR + "/{sample}/variants/ref.fasta.fai", 
    message: "Index the reference with samtools"
    conda: "../envs/polish.yaml"
    log: OUT_DIR + "/logs/varaints/faidx_{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/variants/faidx_{sample}.txt"
    shell: "samtools faidx {input} > {log} 2>& 1"

rule clair3_variants:
    input:
        fai = rules.fadix.output,
        ref = rules.minimap_bams.output.ref, 
        bam = rules.minimap_bams.output.bam,
    output: 
        _dir=directory(OUT_DIR + "/{sample}/variants/clair"),
    message: "Call variants with clair3 [{wildcards.sample}]"
    params:
        platform=config["variants"]["clair"]["platform"],
        model=config["variants"]["clair"]["model"],
        mode=config["variants"]["clair"]["mode"],
    conda: "../envs/clair.yaml"
    log: OUT_DIR + "/logs/varaints/clair/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/variants/clair/{sample}.txt"
    threads: config["threads"]["large"]
    shell:
        """
        run_clair3.sh \
        --bam_fn={input.bam} \
        --ref_fn={input.ref} \
        --threads={threads} \
        --platform={params.platform} \
        --model_path=$CONDA_PREFIX/bin/models/{params.model} \
        --include_all_ctgs \
        {params.mode} \
        --output={output} \
        > {log} 2>&1
        """ 

def get_callers():
     

rule consensus_variants:
    input: rules.minimap_bams.output.ref 
    output: OUT_DIR + "/{sample}/variants/ref.fasta.fai", 
    message: "Take consensus variants from multiple callers."
    conda: "../envs/polish.yaml"
    log: OUT_DIR + "/logs/varaints/faidx_{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/variants/faidx_{sample}.txt"
    shell: "samtools faidx {input} > {log} 2>& 1"

rule vc_finished:
    input: expand(OUT_DIR + "/vcfs/{sample}.vcf", sample=SAMPLES)
    output: OUT_DIR + "/VC_finished"
    message: "Variants calling finished"
    shell: "touch {output}" 