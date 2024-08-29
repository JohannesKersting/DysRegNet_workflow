import os
import pandas as pd
import argparse
import sys

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, required=True,
                        help='meta csv file to downsample')
    parser.add_argument('--outputs', nargs='+', required=True,
                        help='Output paths to downsampled meta files. File name is used to determine the downsample size.')
    args = parser.parse_args()

    print(f"Python version:\n{sys.version}")
    print(f"Arguments:\n{args}\n")

    input_path = args.input
    output_paths = args.outputs

    meta = pd.read_csv(input_path)
    meta_case = meta[meta["condition"] == 1]
    meta_control = meta[meta["condition"] == 0]

    for path in output_paths:
        # get sample size
        size = int(os.path.basename(path).split("-")[0])

        # randomly select samples
        sample = meta_control.sample(n=size, random_state=42)

        # concat with case and save
        pd.concat([sample, meta_case], ignore_index=True).to_csv(path, index=False)


if __name__ == "__main__":
    main()