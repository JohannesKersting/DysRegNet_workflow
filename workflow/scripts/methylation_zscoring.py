import argparse
import sys
import pandas as pd
from scipy.stats import zscore

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--methylation',  type=str, required=True,
                        help='.npy objects with the local methylation test results. Order has to match inputs_global.')
    parser.add_argument('--meta',  type=str, required=True,
                        help='.npy objects with the local methylation test results. Order has to match inputs_global.')
    parser.add_argument('--output',  type=str, required=True,
                        help='Output path to the result .png (will also be saved as .pdf)')

    args = parser.parse_args()

    print(f"Python version:\n{sys.version}")
    print(f"Arguments:\n{args}\n")


    methylation_path = args.methylation
    meta_path = args.meta
    output_path = args.output


    print("\nReading methylation data...")
    methylation = pd.read_csv(methylation_path, index_col=0)
    print(methylation)


    print("\nReading meta data...")
    meta = pd.read_csv(meta_path, sep="\t", index_col=0)
    print(meta)


    # keep only meta samples which have methylation data
    shared_samples = set(meta.index.values).intersection(set(methylation.columns.values))
    meta = meta.loc[list(shared_samples),]
    methylation = methylation[meta.index.values]


    print("\nAvailable samples per cancer type:")
    print(meta["_primary_disease"].value_counts())

    print("\nAvailable control samples per cancer type:")
    print(meta[meta["sample_type"] == "Solid Tissue Normal"]["_primary_disease"].value_counts())

    #  calcutate row-wise z-scores per cacner type
    df_list = []
    for cancer_type, group in meta.groupby("_primary_disease"):
        methylation_group = methylation[group.index.values]
        df_list.append(zscore(methylation_group, axis=1))

    # cbind the pandas for different cancer types
    methylation_zscores = pd.concat(df_list, axis=1)

    # get absolute values
    methylation_zscores = methylation_zscores.abs()

    print("Absolute methylation z-scores:")
    print(methylation_zscores)

    # save output
    methylation_zscores.to_csv(output_path)


if __name__ == "__main__":
    main()