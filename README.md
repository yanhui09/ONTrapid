# Ready-to-use workflow for long read assembly
## Reconciling with multiple assemblers, e.g., [Canu](https://github.com/marbl/canu) and [Flye](https://github.com/fenderglass/Flye).
![rulegraph](/dag.png)
## Rapid assembly soly with a trusted assembler i.e., [Flye](https://github.com/fenderglass/Flye).
![rulegraph-flye](/dag_flye.png)

# Installation
Make sure conda is installed. [Miniconda](https://docs.conda.io/en/latest/miniconda.html) is enough for the whole pipeline.

Install latest version of mamba and snakemake, and activate the conda environment to run snakemake.
```
conda install -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
```

# Usage
Manually edit the working directory path and params in the `config.yaml`.
## Sample Initialization
```
mkdir /path/to/working/directory
python init.py -p /path/to/raw/data -o /path/to/working/directory
```
`/path/to/raw/data` shall be the main directory containing respective runs. E.g. `/mnt/md0/metapredict/raw`,
```
/mnt/md0/metapredict/raw
├── RAPID19
│   └── basecalled_fq
│       ├── barcode01
│       ├── barcode02
│       ├── barcode03
│       ├── barcode04
│       ├── barcode05
│       ├── barcode06
│       ├── barcode07
│       └── barcode08
├── RAPID20
│   └── basecalled_fq
│       ├── barcode06
│       ├── barcode07
│       ├── barcode08
│       ├── barcode09
│       ├── barcode10
│       ├── barcode11
│       └── barcode12
```

## Assembly
```
snakemake -j 24 --use-conda --config work_dir=/path/to/working/directory
```

## Variants calling
```
snakemake -j 6 --use-conda /mnt/md0/nano_snp/out_dh5a/results/variants/vcfs_isc_merged.ann.vcf --config work_dir="/mnt/md0/nano_snp/out_dh5a"
```