# setup cogs and kofams for anvio annotation
rule setup_cogs:
    output: directory(COG)
    conda: "../envs/anvio.yaml"
    log: OUT_DIR + "/logs/init/setup_cogs.log"
    benchmark: OUT_DIR + "/benchmarks/init/setup_cogs.txt"
    threads: config["threads"]["large"]
    shell: "anvi-setup-ncbi-cogs --cog-data-dir {output} -T {threads} --just-do-it > {log} 2>&1"

rule setup_kofams:
    output: directory(KEGG)
    conda: "../envs/anvio.yaml"
    log: OUT_DIR + "/logs/init/setup_kofams.log"
    benchmark: OUT_DIR + "/benchmarks/init/setup_kofams.txt"
    shell: "anvi-setup-kegg-kofams --kegg-data-dir {output} --just-do-it -D > {log} 2>&1"
