import pandas as pd
import glob
import os
import sys
import argparse


def get_input_info(path):
    '''
    Returns a tuple for each input path to extract information from it.
    Order is path, cancer
    '''

    # split path by separator
    path = os.path.normpath(path)
    splitted_path = path.split(os.sep)

    # get cancer and filename
    file_name = splitted_path[-1]
    cancer = splitted_path[-2]

    return path, cancer


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', nargs='+', required=True,
                        help='.tsv meta files')
    parser.add_argument('--output', type=str, required=True,
                        help='Output path to the .tsv summary file (also gets saved as .tex)')

    args = parser.parse_args()

    print(f"Python version:\n{sys.version}")
    print(f"Arguments:\n{args}\n")

    input_paths = args.inputs
    output_path = args.output

    # get input info
    input_infos = [get_input_info(path) for path in input_paths]
    input_df = pd.DataFrame(input_infos, columns=["path", "cancer"])
    print(input_df.to_string())

    # read data and generate a summary
    summary_list = []
    for index, row in input_df.iterrows():
        meta = pd.read_csv(row["path"])
        summary_list.append({
            "Cancer": row["cancer"],
            "Control samples": (meta["condition"] == 0).sum(),
            "Female control samples": meta.query("condition == 0 and gender == 'FEMALE'").shape[0],
            "Male control samples": meta.query("condition == 0 and gender == 'MALE'").shape[0],
            "Case samples": (meta["condition"] == 1).sum(),
        })

    summary_df = pd.DataFrame.from_dict(summary_list)
    summary_list = None

    # save summary
    summary_df.to_csv(output_path, sep="\t", index=False)
    summary_df.style.hide(axis=0).to_latex(buf=os.path.splitext(output_path)[0] + '.tex')


if __name__ == "__main__":
    main()
