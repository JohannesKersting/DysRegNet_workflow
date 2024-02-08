import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
import sys
import colorcet as cc
from commons import method_lut, network_lut


def get_input_info(path):
    """
    Returns a tuple for each input path to extract information from it.
    Order is path, method, norm_method, ref_network_start
    """

    # split path by separator
    path = os.path.normpath(path)
    splitted_path = path.split(os.sep)

    # get method and filename
    file_name = splitted_path[-1]
    method = splitted_path[-2]

    # get overlap_type, norm method and the reference network
    stem_name, file_extension = os.path.splitext(file_name)
    test_type, norm_method, ref_network, min_patients = stem_name.split("-")

    # shorten refernece network name (for genie3 shared and individual)
    ref_network_split = ref_network.split("_")
    if (len(ref_network_split) > 1):
        ref_network_start = "_".join(ref_network_split[0:2])
    else:
        ref_network_start = ref_network_split[0]

    return path, method, norm_method, ref_network_start


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs_local',  nargs='+', required=False,
                        help='.csv objects with the local methylation test results. Order has to match inputs_global. Currently unused.')
    parser.add_argument('--inputs_global',  nargs='+', required=True,
                        help='.csv objects the global methylation test results. Order has to match inputs_local.')
    parser.add_argument('--output',  type=str, required=True,
                        help='Output path to the result .png (will also be saved as .pdf)')

    parser.add_argument('--output_global',  type=str, required=True,
                        help='Output path to the result .png (will also be saved as .pdf)')

    args = parser.parse_args()

    print(f"Python version:\n{sys.version}")
    print(f"Arguments:\n{args}\n")

    input_paths_local = args.inputs_local
    input_paths_global = args.inputs_global
    output_path = args.output
    output_global_path = args.output_global


    input_infos = [get_input_info(path) for path in input_paths_local]

    input_df = pd.DataFrame(input_infos, columns=["local_path", "method", "norm_method", "network"])
    input_df["method"] = input_df["method"].map(lambda x: method_lut[x] if x in method_lut.keys() else x)
    input_df["network"] = input_df["network"].map(lambda x: network_lut[x] if x in network_lut.keys() else x)
    input_df.insert(1, "global_path", input_paths_global)

    print(f"input_df: \n{input_df.to_string()}\n")


    # parse input data
    mutation_test_data = []
    for index, row in input_df.iterrows():
        global_tests_df = pd.read_csv(row["global_path"])
        global_tests_df.sort_values(by=['cancer'], inplace=True)
        global_tests_df["network"] = row["network"]
        global_tests_df["method"] = row["method"]
        mutation_test_data.append(global_tests_df)

    mutation_test_df = pd.concat(mutation_test_data, ignore_index=True)
    mutation_test_data = None

    mutation_test_df["% of significant associations"] = 100 * mutation_test_df["n_sign_fdr"] / mutation_test_df["n_tfs"]
    mutation_test_df["Globally sig."] = np.where(mutation_test_df["global_pval"] < 0.05, 'yes', 'no')
    mutation_test_df = mutation_test_df.rename(columns={"cancer": "Cancer"})

    # drop rows with less than 30 tested tfs
    mutation_test_df = mutation_test_df.query(f'n_tfs >= 30')

    print(f"mutation_test_df: \n{mutation_test_df.to_string()}\n")


    # plot
    sns.set(style="whitegrid")
    y = "% of significant associations"
    col = "network"
    hue = "method"
    x = "Cancer"

    g = sns.catplot(
        mutation_test_df, kind="bar", col_wrap=2,
        x=x, y=y, col=col, hue=hue, palette="Set1",
        height=6, aspect=1.5, legend=None, sharex=False,
    )
    g.set_titles(col_template="{col_name}", row_template="{row_name}")
    g.add_legend(title="Method")

    # iterate through subplots
    for i, ax in enumerate(g.axes.ravel()):
        # iterate through containers (one per hue)
        for j, c in enumerate(ax.containers):
            sub_df = mutation_test_df.query(f'{hue}=="{c.get_label()}" and {col}=="{ax.get_title()}"')

            labels = [u'\u2217 ' if row["global_pval"] < 0.05 else "" for index, row in sub_df.iterrows()]
            print(f"{ax.get_title()} {c.get_label()} {len(labels)}")
            ax.bar_label(c, labels=labels)

    plt.savefig(output_path, dpi=300)
    plt.savefig(os.path.splitext(output_path)[0] + '.pdf', dpi=300)


    # global p-values plot
    sns.set(style="whitegrid")

    mutation_test_df["transformed_pval"] = -np.log10(mutation_test_df["global_pval"])

    g = sns.catplot(
        mutation_test_df, kind="bar",
        x="Cancer", y="transformed_pval", col="network", hue="method", palette="Set1", col_wrap=2,
        height=6, aspect=1.5, legend=None, sharex=False, sharey=False
    )
    g.set_titles(col_template="{col_name}", row_template="{row_name}")
    g.set(ylabel=r'$-\log_{10}$(p-value)')
    g.add_legend(title="Method")
    g.refline(y=-np.log10(0.05))

    plt.savefig(output_global_path, dpi=300)
    plt.savefig(os.path.splitext(output_global_path)[0] + '.pdf', dpi=300)


if __name__ == "__main__":
    main()