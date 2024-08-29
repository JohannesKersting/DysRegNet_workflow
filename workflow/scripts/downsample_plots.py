import os
import pandas as pd
import pyarrow
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', nargs='+', required=True,
                        help='.csv files with network stats')
    parser.add_argument('--output', type=str, required=True,
                        help='Output path to the result .png (will also be saved as .pdf)')
    args = parser.parse_args()

    print(f"Python version:\n{sys.version}")
    print(f"Arguments:\n{args}\n")

    input_paths = args.inputs
    output_path = args.output


    # read and format data
    df_list = []
    for path in input_paths:
        n_samples = int(os.path.basename(path).split("-")[0])
        df = pd.read_csv(path, index_col=0)
        df["n_samples"] = n_samples
        df_list.append(df)
    df = pd.concat(df_list, ignore_index=True)
    df_list = None


    # plot
    sns.set(style="whitegrid")
    ax = sns.violinplot(df, x="n_samples", y="pval_birth_days_to", cut=0, color="#127611")
    ax.axhline(0.05)
    plt.gca().invert_xaxis()
    ax.set_xlabel("Number of control samples")
    ax.set_ylabel('P-value of the "age" coefficient')
    plt.gcf().tight_layout()

    plt.savefig(output_path, dpi=300)
    plt.savefig(os.path.splitext(output_path)[0] + '.pdf', dpi=300)

if __name__ == "__main__":
    main()