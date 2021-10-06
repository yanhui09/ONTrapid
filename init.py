import os
from collections import defaultdict
import pandas as pd
import argparse

def add_sample(sample_dict, sample_id, header, path):
    sample_dict[sample_id][header]=path

def get_samples(path):
    """
    create a table containing the paths to basecalled dir of each barcodes
    """
    samples = defaultdict(dict)
    for root, dirs, fqs in os.walk(os.path.abspath(path)):
        for dirname in dirs:

            # only check fq files
            if "barcode" in dirname:
                dir_path = os.path.join(root, dirname)
                root_split = root.split(os.sep)
                sample_barcode =  root_split[-2] + "_" + dirname # specify by run
                add_sample(samples, sample_barcode, "fq_dir", dir_path)
    samples_dt = pd.DataFrame(samples).T
    return samples_dt

def parse_arguments():
    """Read arguments from the console"""
    parser = argparse.ArgumentParser(description="Note: generate sample.tsv")
    parser.add_argument("-p", "--path", help='path to raw data')
    parser.add_argument("-o", "--out", help='path to working directory')

    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    sample_dt = get_samples(args.path)
    sample_dt.to_csv(args.out + "/samples.tsv", sep="\t")

if __name__ == "__main__":
    main()

for root, dirs, files in os.walk("raw", topdown=False):
    for name in dirs:
       print(os.path.join(root, name))