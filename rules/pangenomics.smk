# check the existence of the requant directory
def dir_check(dir_path, ext):
    if not os.path.exists(dir_path):
        raise ValueError("\n  Directory not found.\n\tMake sure {} is loaded.\n".format(dir_path))
    else:
        fs = [f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))]
        # file with suffix .fasta, fa, or fna
        fas = [f for f in fs if f.endswith(ext)]
        if len(fas) == 0:
            raise ValueError("\n  The {} files \n are not found in {}.\n".format(ext, dir_path))
        else:
            return True

checkpoint update_assembly:
    output: directory(OUT_DIR + '/assembly_updated')
    run:
        dir_check(config["updated_assembly_dir"], ".fasta")
        shutil.copytree(config["updated_assembly_dir"], output[0])

rule external_gff3s:
    output: directory(OUT_DIR + '/external_gff3s')
    run:
        dir_check(config["external_gff3_dir"], ".gff")
        shutil.copytree(config["external_gff3_dir"], output[0]) 

def get_pansamples(wildcards, update):
    if update:
        return glob_wildcards(checkpoints.update_assembly.get(**wildcards).output[0] + "/{sample}.fasta").sample
    else:
        return glob_wildcards(checkpoints.collect_assembly.get(**wildcards).output[0] + "/{sample}.fasta").sample

def get_assembly_dir(update):
    if update:
        return OUT_DIR + '/assembly_updated'
    else:
        return OUT_DIR + '/assembly'

# build anvio contig database
# clean the fasta header if possible, especially for NCBI ref
rule reformat_fasta:
    input: get_assembly_dir(config["update_assembly"]) + "/{sample}.fasta" 
    output: OUT_DIR + "/pangenomics/reformated/{sample}.fasta"
    conda: "../envs/anvio.yaml"
    log: OUT_DIR + "/logs/pangenomics/reformat_fasta/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/pangenomics/reformat_fasta/{sample}.txt"
    shell: "anvi-script-reformat-fasta {input} -o {output} -l 0 --simplify-names > {log} 2>&1"

# generate anvio contig db
rule contig_db:
    input: ancient(OUT_DIR + "/pangenomics/reformated/{sample}.fasta")
    output: OUT_DIR + "/pangenomics/contig_db/{sample}.db"
    conda: "../envs/anvio.yaml"
    log: OUT_DIR + "/logs/pangenomics/contig_db/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/pangenomics/contig_db/{sample}.txt"
    threads: config["threads"]["normal"]
    shell: "anvi-gen-contigs-database -f {input} -o {output} -n 'proj' -T {threads} > {log} 2>&1"

# generate anvio with
# gff_parser
rule gff_parse:
    input: OUT_DIR + "/external_gff3s/{sample}.gff"
    output:
        gene_calls = OUT_DIR + "/pangenomics/gffparse/gene_calls_{sample}.txt",
        gene_annot = OUT_DIR + "/pangenomics/gffparse/gene_annot_{sample}.txt",
    params:
        source = "Prokka",
    conda: "../envs/anvio.yaml"
    log: OUT_DIR + "/logs/pangenomics/gff_parse/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/pangenomics/gff_parse/{sample}.txt"
    shell: 
        "python {workflow.basedir}/scripts/gff_parser.py {input} "
        "{output.gene_calls} {output.gene_annot} --source {params.source} > {log} 2>&1"

rule contig_db_external:
    input: 
        fasta = ancient(OUT_DIR + "/pangenomics/reformated/{sample}.fasta"),
        gene_calls = rules.gff_parse.output.gene_calls,
        gene_annot = rules.gff_parse.output.gene_annot,
    output: OUT_DIR + "/pangenomics/contig_db_external/{sample}.db"
    conda: "../envs/anvio.yaml"
    log: OUT_DIR + "/logs/pangenomics/contig_db_external/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/pangenomics/contig_db_external/{sample}.txt"
    threads: config["threads"]["normal"]
    shell: 
        """
        anvi-gen-contigs-database -f {input.fasta} -o {output} -n 'proj' -T {threads} --external-gene-calls {input.gene_calls} > {log} 2>&1
        anvi-import-functions -c {output} -i {input.gene_annot} >> {log} 2>&1
        """

# scheduler to include external gene calls
def get_contig_db(external_gene_calls=config['external_gene_calls']):
    if not external_gene_calls:
        return rules.contig_db.output
    else:
        return rules.contig_db_external.output

# hmms calculation
rule run_hmms:
    input: ancient(get_contig_db)
    output: OUT_DIR + "/pangenomics/.hmms/.{sample}_DONE"
    conda: "../envs/anvio.yaml"
    log: OUT_DIR + "/logs/pangenomics/run_hmms/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/pangenomics/run_hmms/{sample}.txt"
    threads: config["threads"]["normal"]
    shell: "anvi-run-hmms -c {input} --also-scan-trnas -T {threads} > {log} 2>&1 && touch {output}"

# COG annoatation
rule run_cogs:
    input: ancient(get_contig_db)
    output: OUT_DIR + "/pangenomics/.cogs/.{sample}_DONE"
    params:
        COG = COG,
    conda: "../envs/anvio.yaml"
    log: OUT_DIR + "/logs/pangenomics/run_cogs/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/pangenomics/run_cogs/{sample}.txt"
    threads: config["threads"]["normal"]
    shell: "anvi-run-ncbi-cogs -c {input} --cog-data-dir {params.COG} -T {threads} > {log} 2>&1 && touch {output}"

# KEGG annotation
rule run_kofams:
    input: ancient(get_contig_db)
    output: OUT_DIR + "/pangenomics/.kofams/.{sample}_DONE"
    params:
        KEGG = KEGG,
    conda: "../envs/anvio.yaml"
    log: OUT_DIR + "/logs/pangenomics/run_kofams/{sample}.log"
    benchmark: OUT_DIR + "/benchmarks/pangenomics/run_kofams/{sample}.txt"
    threads: config["threads"]["normal"]
    shell: "anvi-run-kegg-kofams -c {input} --kegg-data-dir {params.KEGG} -T {threads} > {log} 2>&1 && touch {output}"

# create the GENOME storage
rule gen_genomes_list:
    input: ancient(lambda wc: expand(get_contig_db, sample=get_pansamples(wc, config["update_assembly"])))
    output: OUT_DIR + "/pangenomics/genomes.txt"
    run:
        genomes = [x.split("/")[-1].split(".")[0] for x in input]
        with open(output[0], "w") as f:
            f.write("name\tcontigs_db_path\n")
            f.write("\n".join(["\t".join(x) for x in zip(genomes, input)]))  
     
rule gen_genomes_storage:
    input: 
      lambda wc: expand(OUT_DIR + "/pangenomics/.cogs/.{sample}_DONE", sample=get_pansamples(wc, config["update_assembly"])),
      lambda wc: expand(OUT_DIR + "/pangenomics/.kofams/.{sample}_DONE", sample=get_pansamples(wc, config["update_assembly"])),
      lambda wc: expand(OUT_DIR + "/pangenomics/.hmms/.{sample}_DONE", sample=get_pansamples(wc, config["update_assembly"])),
      glist = rules.gen_genomes_list.output,
    output: OUT_DIR + "/pangenomics/GENOMES.db"
    conda: "../envs/anvio.yaml"
    log: OUT_DIR + "/logs/pangenomics/gen_genomes_storage.log"
    benchmark: OUT_DIR + "/benchmarks/pangenomics/gen_genomes_storage.txt"
    shell: "anvi-gen-genomes-storage -e {input.glist} -o {output}"

# pangenomics analysis
rule pan_genome:
    input: rules.gen_genomes_storage.output,
    output: 
        db = OUT_DIR + "/pangenomics/pan_genomes/proj-PAN.db",
        _dir = directory(OUT_DIR + "/pangenomics/pan_genomes"),
    conda: "../envs/anvio.yaml"
    log: OUT_DIR + "/logs/pangenomics/pan_genome.log"
    benchmark: OUT_DIR + "/benchmarks/pangenomics/pan_genome.txt"
    threads: config["threads"]["large"]
    shell:
        """
        anvi-pan-genome -g {input} \
            --project-name 'proj' \
            --output-dir {output._dir} \
            --num-threads {threads} \
            --minbit 0.5 \
            --mcl-inflation 10 \
            > {log} 2>&1
        """

# incoporated with clade data
#anvi-import-misc-data layer-additional-data.txt \
#                      -p meta-pan/metapredict-PAN.db \
#                      --target-data-table layers

# calculate the average genome-similarity
# Install coverage python module
# AttributeError: type object 'DataFrame' has no attribute 'from_csv'
#pip install pandas==0.23.1

rule calc_genome_similarity:
    input:
       glist = rules.gen_genomes_list.output, 
       db = rules.pan_genome.output.db,
    output: directory(OUT_DIR + "/pangenomics/ANI")
    conda: "../envs/anvio.yaml"
    log: OUT_DIR + "/logs/pangenomics/calc_genome_similarity.log"
    benchmark: OUT_DIR + "/benchmarks/pangenomics/calc_genome_similarity.txt"
    threads: config["threads"]["large"]
    shell:
        """
        anvi-compute-genome-similarity --external-genomes {input.glist} \
            --program pyANI \
            --output-dir {output} \
            --num-threads {threads} \
            --pan-db {input.db} \
            > {log} 2>&1
        """

# display
#anvi-display-pan -g metapredict-GENOMES.db -p meta-pan/metapredict-PAN.db --server-only -P 8080

# genearate a summary table 
# add a default collection containing every gene cluster
rule add_default_collection:
    input: rules.pan_genome.output.db,
    output: OUT_DIR + "/pangenomics/.add_default_collection_DONE",
    conda: "../envs/anvio.yaml"
    log: OUT_DIR + "/logs/pangenomics/add_default_collection.log"
    benchmark: OUT_DIR + "/benchmarks/pangenomics/add_default_collection.txt"
    shell: "anvi-script-add-default-collection -p {input} > {log} 2>&1 && touch {output}"

rule anvi_sum:
    input: 
      rules.add_default_collection.output,
      gdb = rules.gen_genomes_storage.output,
      pandb = OUT_DIR + "/pangenomics/pan_genomes/proj-PAN.db",
    output: directory(OUT_DIR + "/pangenomics/pan_genomes/summary_gcs"),
    conda: "../envs/anvio.yaml"
    log: OUT_DIR + "/logs/pangenomics/anvi_sum.log"
    benchmark: OUT_DIR + "/benchmarks/pangenomics/anvi_sum.txt"
    shell:
        """
        anvi-summarize -p {input.pandb} \
            -g {input.gdb} \
            -C DEFAULT \
            -o {output} \
            > {log} 2>&1           
        """