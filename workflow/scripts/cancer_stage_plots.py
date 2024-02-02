import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os
import sys

from commons import method_lut, network_lut


def get_input_info(path):
    '''
    Returns a tuple for each input path to extract information from it.
    Order is path, cancer, method, norm_method, ref_network_start
    '''

    # split path by separator
    path = os.path.normpath(path)
    splitted_path = path.split(os.sep)

    # the cancer type is ecpected to be in the third last postion

    # get cancer, method and filename
    file_name = splitted_path[-1]
    method = splitted_path[-2]
    cancer = splitted_path[-3]

    # get overlap_type, norm method and the reference network
    stem_name, file_extension = os.path.splitext(file_name)
    test_type, norm_method, ref_network = stem_name.split("-")

    # shorten refernece network name (for genie3 shared and individual)
    ref_network_split = ref_network.split("_")
    if (len(ref_network_split) > 1):
        ref_network_start = "_".join(ref_network_split[0:2])
    else:
        ref_network_start = ref_network_split[0]

    return path, cancer, method, norm_method, ref_network_start


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', nargs='+', required=True,
                        help='.csv files with test results')
    parser.add_argument('--output', type=str, required=True,
                        help='Output path to the result .png (will also be saved as .pdf)')

    args = parser.parse_args()

    print(f"Python version:\n{sys.version}")
    print(f"Arguments:\n{args}\n")

    input_paths = args.inputs
    output_path = args.output

    # get input info
    input_infos = [get_input_info(path) for path in input_paths]
    input_df = pd.DataFrame(input_infos, columns=["path", "cancer", "method", "norm_method", "network"])
    input_df["method"] = input_df["method"].map(lambda x: method_lut[x] if x in method_lut.keys() else x)
    input_df["network"] = input_df["network"].map(lambda x: network_lut[x] if x in network_lut.keys() else x)

    print(f"input_df: \n{input_df.to_string()}\n")

    # read input files
    summary_list = []
    for index, row in input_df.iterrows():

        # ignore empty files
        if os.path.getsize(row["path"]) == 0:
            continue

        test_data = pd.read_csv(row["path"], index_col=0)
        summary_list.append({
            "% of significant associations": 100 * sum(test_data["fdr_bh"] < 0.05) / test_data.shape[0],
            "cancer": row["cancer"],
            "method": row["method"],
            "network": row["network"]
        })
    summary_df = pd.DataFrame(summary_list)
    print(f"summary_df: \n{summary_df.to_string()}\n")

    # plot
    sns.set(style="whitegrid")

    y = "cancer"
    col = "network"
    hue = "method"

    g = sns.catplot(
        summary_df, kind="bar",
        x="% of significant associations", y=y, col=col, hue=hue, palette="Set1",
        height=6, aspect=0.5, legend=None
    )
    g.set_titles(col_template="{col_name}", row_template="{row_name}")
    g.add_legend(title="Method")

    plt.savefig(output_path, dpi=300)
    plt.savefig(os.path.splitext(output_path)[0] + '.pdf', dpi=300)


if __name__ == "__main__":
    main()
