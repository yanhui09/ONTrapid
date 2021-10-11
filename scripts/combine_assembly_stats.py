import pandas as pd
import os 

def busco2tdf(file):
    stats_lines = []
    with open(file) as f:
        for i, line in enumerate(f):
            if "***** Results: *****" in line:
                num = i
        for i, line in enumerate(f):    
            if num + 2 <= i <= num +8:
                stats_lines.append(line.lstrip())

    val = [x.split("\t")[0] for x in stats_lines]
    _index = ["Percentage"] + [x.split("\t")[1] for x in stats_lines if "\t" in x]
    df = pd.DataFrame(val, index=_index)
    return df

def main(input, quast_out, busco_out):
    quast_stats = pd.DataFrame()
    busco_stats = pd.DataFrame()
    for file in input:
        # quast summary
        if "quast" in file:
            d = pd.read_csv(file + "/report.tsv", sep="\t", skiprows=1, index_col=0, header=None)
            quast_stats = quast_stats.append(d.T)
        elif "busco" in file:
            assembly_out=os.path.split(file)[-1]
            d = busco2tdf(file + "/summary.txt")
            busco_stats = busco_stats.append(d.T)
        else:
            raise Exception("QUAST or BUSCO output is not found.")
    
    if not quast_stats.empty:
        quast_stats.index = [os.path.split(x)[-1] for x in input if "quast" in x]
        quast_stats.T.to_csv(quast_out, sep="\t")
    if not busco_stats.empty:
        busco_stats.index = [os.path.split(x)[-1] for x in input if "busco" in x]
        busco_stats.T.to_csv(busco_out, sep="\t")

if __name__ == "__main__":
    main(
        input=snakemake.input,
        quast_out=snakemake.output.quast,
        busco_out=snakemake.output.busco,
    )