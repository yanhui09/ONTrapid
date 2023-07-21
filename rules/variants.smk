REF=config["variants"]["ref"]

rule medaka_haploid_variants:
    input:
        fastq = rules.porechop.output,
    output: 
        vcf = OUT_DIR + "/variants/{sample}/medaka/medaka.annotated.vcf",
    message: "Call variants with medaka [{wildcards.sample}]"
    params:
        m=config["variants"]["medaka"]["model"],
        _dir=OUT_DIR + "/variants/{sample}/medaka",
        fasta = REF, 
    conda: "../envs/medaka.yaml"
    log: OUT_DIR + "/logs/varaints/medaka/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/variants/medaka/{sample}.txt"
    threads: config["threads"]["large"]
    shell:
        """
        medaka_haploid_variant -i {input.fastq} \
        -r {params.fasta} -o {params._dir}  \
        -m {params.m} -s -t {threads} \
        > {log} 2>&1
        rm {params.fasta}.mmi {params.fasta}.fai
        """

# Clair3 https://github.com/HKU-BAL/Clair3
# Pepper https://github.com/kishwarshafin/pepper # Designed for diploid, no conda 
# Magalodon https://github.com/nanoporetech/megalodon # FASTA5 required
# Nanopolish https://github.com/jts/nanopolish # FASTA5 signals required
# Nanocaller https://github.com/WGLab/NanoCaller # Haploid mode not supported yet
rule faidx:
    input:
        fasta = REF, 
    output:
        ref = OUT_DIR + "/variants/ref.fasta",
        fai = OUT_DIR + "/variants/ref.fasta.fai", 
    message: "Index the reference with samtools"
    conda: "../envs/racon.yaml"
    log: OUT_DIR + "/logs/varaints/faidx.log"
    benchmark: OUT_DIR + "/benchmarks/variants/faidx.txt"
    shell: 
        """
        [ ! -f {output.ref} ] && cp {input.fasta} {output.ref}
        samtools faidx {output.ref} > {log} 2>& 1
        """

rule minimap_bams:
    input:
        ref = rules.faidx.output.ref,
        fastq = rules.porechop.output,
    output:
        bam = OUT_DIR + "/variants/{sample}/ref.sorted.bam",
        bai = OUT_DIR + "/variants/{sample}/ref.sorted.bam.bai",
    message: "Generate bam file with minimap2 [{wildcards.sample}]"
    params:
        _dir=OUT_DIR + "/variants/{sample}/clair",
        x=config["minimap"]["x"],
    conda: "../envs/racon.yaml"
    log: OUT_DIR + "/logs/varaints/minimap_bams/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/variants/minimap_bams/{sample}.txt"
    threads: config["threads"]["large"]
    shell:
        """
        minimap2 -L -t {threads} -ax {params.x} {input.ref} {input.fastq} 2> {log} | \
        samtools sort -o {output.bam} 2>> {log}
        samtools index {output.bam} 2>> {log}
        """ 

rule clair3_variants:
    input:
        bam = rules.minimap_bams.output.bam,
    output: 
        vcf = OUT_DIR + "/variants/{sample}/clair/merge_output.vcf",
    message: "Call variants with clair3 [{wildcards.sample}]"
    params:
        _dir=OUT_DIR + "/variants/{sample}/clair",
        platform=config["variants"]["clair"]["platform"],
        model=config["variants"]["clair"]["model"],
        mode=config["variants"]["clair"]["mode"],
        ref=rules.faidx.output.ref,
    conda: "../envs/clair.yaml"
    log: OUT_DIR + "/logs/varaints/clair/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/variants/clair/{sample}.txt"
    threads: config["threads"]["large"]
    shell:
        """
        run_clair3.sh \
        --bam_fn={input.bam} \
        --ref_fn={params.ref} \
        --threads={threads} \
        --platform={params.platform} \
        --model_path=$CONDA_PREFIX/bin/models/{params.model} \
        --include_all_ctgs --remove_intermediate_dir \
        --no_phasing_for_fa {params.mode} \
        --output={params._dir} \
        > {log} 2>&1
        gunzip {params._dir}/merge_output.vcf.gz
        """ 

rule vcf_sort:
    input: rules.medaka_haploid_variants.output.vcf,
    output: OUT_DIR + "/variants/{sample}/medaka/medaka.annotated.sorted.vcf",
    message: "Sort vcf file [{wildcards.sample}]"
    conda: "../envs/bcftools.yaml"
    log: OUT_DIR + "/logs/varaints/bcftools/sort/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/variants/bcftools/sort/{sample}.txt"
    shell: "bcftools sort {input} > {output} 2> {log}"

rule vcf_comp_index:
    input: rules.vcf_sort.output,
    output:
        vcf_gz=OUT_DIR + "/variants/{sample}/medaka/medaka.annotated.sorted.vcf.gz",
        vcf_csi=OUT_DIR + "/variants/{sample}/medaka/medaka.annotated.sorted.vcf.gz.csi",
    message: "Compress and tbi index vcf file [{wildcards.sample}]"
    conda: "../envs/bcftools.yaml"
    log: OUT_DIR + "/logs/varaints/bcftools/comp_index/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/variants/bcftools/comp_index/{sample}.txt"
    shell:
        """
        bgzip -c {input} > {output.vcf_gz} 2> {log}
        bcftools index {output.vcf_gz} >> {log} 2>&1
        """

use rule vcf_comp_index as vcf_comp_index2 with:
    input:
        rules.clair3_variants.output,
    output:
        vcf_gz=OUT_DIR + "/variants/{sample}/clair/merge_output.vcf.gz",
        vcf_csi=OUT_DIR + "/variants/{sample}/clair/merge_output.vcf.gz.csi",
    log:
        OUT_DIR + "/logs/varaints/bcftools/comp_index2/{sample}.log"
    benchmark:
        OUT_DIR + "/benchmarks/variants/bcftools/comp_index2/{sample}.txt"

rule vcf_isec:
    input:
        rules.vcf_comp_index.output.vcf_gz,
        rules.vcf_comp_index2.output.vcf_gz,
    output:
        isec_medaka=OUT_DIR + "/variants/{sample}/isec/0002.vcf.gz",
        isec_clair=OUT_DIR + "/variants/{sample}/isec/0003.vcf.gz",
    message: "Take intersections of vcf files [{wildcards.sample}]"
    params:
        _dir=OUT_DIR + "/variants/{sample}/isec",
    conda: "../envs/bcftools.yaml"
    log: OUT_DIR + "/logs/varaints/bcftools/isec/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/variants/bcftools/isec/{sample}.txt"
    shell:
        "bcftools isec -Oz -p {params._dir} {input} 2> {log}"

# bcftools reheader
rule vcf_reheader_reindex:
    input: rules.vcf_isec.output.isec_clair
    output:
        vcf=OUT_DIR + "/variants/{sample}/isec/0003_rename.vcf.gz",
        csi=OUT_DIR + "/variants/{sample}/isec/0003_rename.vcf.gz.csi",
    message: "Reheader vcf file [{wildcards.sample}]"
    params:
        s=OUT_DIR + "/variants/{sample}/isec/s",
    conda: "../envs/bcftools.yaml"
    log: OUT_DIR + "/logs/varaints/bcftools/reheader_reindex/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/variants/bcftools/reheader_reindex/{sample}.txt"
    shell:
        """
        echo "SAMPLE {wildcards.sample}\n" > {params.s}
        bcftools reheader -s {params.s} -o {output.vcf} {input} > {log} 2>&1
        bcftools index {output.vcf} >> {log} 2>&1
        rm -f {params.s}
        """

rule vcf_merge:
    input: expand(OUT_DIR + "/variants/{sample}/isec/0003_rename.vcf.gz", sample=SAMPLES)
    output: OUT_DIR + "/variants/vcfs_isc_merged.vcf"
    message: "Merge intersected vcf files"
    conda: "../envs/bcftools.yaml"
    log: OUT_DIR + "/logs/varaints/bcftools/vcf_merge.log"
    benchmark: OUT_DIR + "/benchmarks/variants/bcftools/vcf_merge.txt"
    shell: 
        """
        bcftools merge \
        --merge all {input} -Ov -o {output} > {log} 2>&1
        """

rule variant_annotate:
    input: 
        ref = rules.faidx.output.ref, 
        vcf = rules.vcf_merge.output,
    output:
        vcf2snpEff = OUT_DIR + "/variants/snpEff/vcfs_isc_merged2snpEff.vcf",
        snpEff = OUT_DIR + "/variants/vcfs_isc_merged.ann.vcf",
        _dir=directory(OUT_DIR+"/variants/snpEff"),
    message: "Annotate the merged vcf file with SnpEff"
    conda: "../envs/snpeff.yaml"
    params:
        g_version=config["variants"]["snpEff"]["g_version"],
        chr_exist=config["variants"]["snpEff"]["chr_exist"],
        chr_replace=config["variants"]["snpEff"]["chr_replace"],
    log: OUT_DIR + "/logs/varaints/snpEff/VariantAnnotate.log"
    benchmark: OUT_DIR + "/benchmarks/variants/snpEff/VariantAnnotate.txt"
    shell:
        """
        sed 's/{params.chr_exist}/{params.chr_replace}/' {input.vcf} > {output.vcf2snpEff}
        cd {output._dir}
        snpEff ann -i vcf -o gatk {params.g_version} {output.vcf2snpEff} > {output.snpEff} 2> {log}
        cd - > /dev/null 2>&1
        """

# To do:
# build custom snpEff database
# https://www.biostars.org/p/432180/
# Probably use the consensus assembly to call snp ? a custrom reference?