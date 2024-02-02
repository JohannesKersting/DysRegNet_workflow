import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pyarrow
import argparse
import os
import sys

from commons import method_lut, network_lut



def get_input_info(path):
    '''
    Returns a tuple for each input path to extract information from it.
    Order is path, method, norm_method, ref_network_start
    '''

    # split path by separator
    path = os.path.normpath(path)
    splitted_path = path.split(os.sep)

    # get method and filename
    file_name = splitted_path[-1]
    method = splitted_path[-2]

    # get overlap_type, norm method and the reference network
    stem_name, file_extension = os.path.splitext(file_name)
    test_type, norm_method, ref_network = stem_name.split("-")

    # shorten refernece network name (for genie3 shared and individual)
    ref_network_split = ref_network.split("_")
    if (len(ref_network_split) > 1):
        ref_network_start = "_".join(ref_network_split[0:2])
    else:
        ref_network_start = ref_network_split[0]

    return path, method, norm_method, ref_network_start


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs_local',  nargs='+', required=True,
                        help='.npy objects with the local methylation test results. Order has to match inputs_global.')
    parser.add_argument('--inputs_global',  nargs='+', required=True,
                        help='.npy objects with the global methylation test results. Order has to match inputs_local.')
    parser.add_argument('--output',  type=str, required=True,
                        help='Output path to the result .png (will also be saved as .pdf)')

    args = parser.parse_args()

    print(f"Python version:\n{sys.version}")
    print(f"Arguments:\n{args}\n")

    input_paths_local = args.inputs_local
    input_paths_global = args.inputs_global
    output_path = args.output

    input_infos = [get_input_info(path) for path in input_paths_local]

    input_df = pd.DataFrame(input_infos, columns=["local_path", "method", "norm_method", "network"])
    input_df["method"] = input_df["method"].map(lambda x: method_lut[x] if x in method_lut.keys() else x)
    input_df["network"] = input_df["network"].map(lambda x: network_lut[x] if x in network_lut.keys() else x)
    input_df.insert(1, "global_path", input_paths_global)

    print(f"input_df: \n{input_df.to_string()}\n")


    # read and format input data
    methylation_test_data = []
    for index, row in input_df.iterrows():
        local_tests = np.load(row["local_path"], allow_pickle='TRUE').item()
        global_tests = np.load(row["global_path"], allow_pickle='TRUE').item()

        for cancer_type in sorted(global_tests.keys()):
            global_pval, random_pval, n_tests = global_tests[cancer_type]

            # local_tests[cancer_type] is list with tuples (tg, p-val corrected, p-val)
            pvals_corrected = [pval_corrected for tg, pval_corrected, p_val in local_tests[cancer_type]]
            n_signif = len([pval for pval in pvals_corrected if pval < 0.05])

            result = (
                row["network"],
                row["method"],
                cancer_type,
                n_signif,
                n_tests,
                global_pval
            )
            methylation_test_data.append(result)

    methylation_test_df = pd.DataFrame(methylation_test_data,
                                       columns=["network", "method", "cancer", "number significant", "number of tests",
                                                "global p-value"]
                                       )
    methylation_test_data = None
    methylation_test_df["% of significant associations"] = 100 * methylation_test_df["number significant"] / \
                                                           methylation_test_df["number of tests"]

    print(f"methylation_test_df: \n{methylation_test_df.to_string()}\n")


    # plot
    sns.set(style="whitegrid")
    y = "cancer"
    col = "network"
    hue = "method"

    g = sns.catplot(
        methylation_test_df, kind="bar",
        x="% of significant associations", y=y, col=col, hue=hue, palette="Set1",
        height=6, aspect=0.5, legend=None
    )
    g.set_titles(col_template="{col_name}", row_template="{row_name}")
    g.add_legend(title="Method")

    # iterate through subplots
    for i, ax in enumerate(g.axes.ravel()):
        # iterate through containers (one per hue)
        for j, c in enumerate(ax.containers):
            sub_df = methylation_test_df.query(f'{hue}=="{c.get_label()}" and {col}=="{ax.get_title()}"')

            labels = [u'\u2217 ' if row["global p-value"] < 0.05 else "" for index, row in sub_df.iterrows()]
            ax.bar_label(c, labels=labels)

    plt.savefig(output_path, dpi=300)
    plt.savefig(os.path.splitext(output_path)[0] + '.pdf', dpi=300)



if __name__ == "__main__":
    main()