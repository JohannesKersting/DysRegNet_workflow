import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import argparse
import os
import sys

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--input',  nargs='+', required=True,
                        help='Snakemake benchmark files, tile name will be used to identify the method.')
    parser.add_argument('--output',  type=str, required=True,
                        help='Output path to the result .png (will also be saved as .pdf)')

    args = parser.parse_args()

    print(f"Python version:\n{sys.version}")
    print(f"Arguments:\n{args}\n")

    inputs = args.input
    output_path = args.output


    # Read data
    df_list = []
    method_lut = {"dysregnet": "DysRegNet",
                  "ssn": "SSN",
                  }

    for input_path in inputs:
        df = pd.read_csv(input_path, sep="\t")
        # file name is assumed to be the method name
        df["method"] = os.path.splitext(os.path.basename(input_path))[0]
        df_list.append(df)

    df = pd.concat(df_list, ignore_index=True)
    # map method name
    df["method"] = df["method"].map(method_lut)
    print(df)

    print("\nMean run-time [s]:")
    print(df.groupby("method")["s"].mean().to_string())


    # Plot
    sns.set(style="whitegrid")

    sns.boxplot(data=df, x="method", y="s")
    plt.ylabel("Time [s]")
    plt.xlabel("Method")
    plt.tight_layout()

    plt.savefig(output_path, dpi=300)
    plt.savefig(os.path.splitext(output_path)[0] + '.pdf', dpi=300)

if __name__ == "__main__":
    main()