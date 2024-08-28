import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
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

    # get method and filename
    file_name = splitted_path[-1]
    method = splitted_path[-2]
    cancer = splitted_path[-3]

    # get  norm method and the reference network
    stem_name, file_extension = os.path.splitext(file_name)
    _, norm_method, ref_network = stem_name.split("-")

    # shorten reference network name (for genie3 shared and individual)
    ref_network_split = ref_network.split("_")
    if (len(ref_network_split) > 1):
        ref_network_start = "_".join(ref_network_split[0:2])
    else:
        ref_network_start = ref_network_split[0]

    return path, cancer, method, norm_method, ref_network_start



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', nargs='+', required=True,
                        help='.csv files with network stats')
    parser.add_argument('--output_edges', type=str, required=True,
                        help='Output path to the result .png (will also be saved as .pdf)')
    parser.add_argument('--output_nodes', type=str, required=True,
                        help='Output path to the result .png (will also be saved as .pdf)')
    parser.add_argument('--output_components', type=str, required=True,
                        help='Output path to the result .png (will also be saved as .pdf)')
    parser.add_argument('--output_nodes_per_component', type=str, required=True,
                        help='Output path to the result .png (will also be saved as .pdf)')

    args = parser.parse_args()

    print(f"Python version:\n{sys.version}")
    print(f"Arguments:\n{args}\n")

    input_paths = args.inputs
    output_edges_path = args.output_edges
    output_nodes_path = args.output_nodes
    output_components_path = args.output_components
    output_nodes_per_component_path = args.output_nodes_per_component


    # get input info
    input_infos = [get_input_info(path) for path in input_paths]
    input_df = pd.DataFrame(input_infos, columns=["path", "cancer", "method", "norm_method", "network"])

    input_df["method"] = input_df["method"].map(lambda x: method_lut[x] if x in method_lut.keys() else x)
    input_df["network"] = input_df["network"].map(
        lambda x: "genie3_individual" if (("genie3" in x) and (not "shared" in x)) else x)
    input_df["network"] = input_df["network"].map(lambda x: network_lut[x] if x in network_lut.keys() else x)
    print(input_df.to_string())


    # read files
    df_list = []
    for index, row in input_df.iterrows():
        df = pd.read_csv(row["path"])
        df["cancer"] = row["cancer"]
        df["method"] = row["method"]
        df["norm_method"] = row["norm_method"]
        df["network"] = row["network"]
        df_list.append(df)
    df = pd.concat(df_list, ignore_index=True)
    df_list = None
    df["mean_nodes_per_component"] = df["nodes"] / df["connected_components"]

    print(df)


    # plot edges
    sns.set(style="whitegrid")
    g = sns.catplot(
        df, kind="box",
        y="edges", x="cancer", col="network", col_wrap=2, hue="method", palette="Set1",
        height=3, aspect=2, legend=None, sharey=False
    )
    g.set_titles(col_template="{col_name}", row_template="{row_name}")
    [plt.setp(ax.get_xticklabels(), rotation=45) for ax in g.axes.flat]
    g.add_legend(title="Method")
    g.set_xlabels("Cancer")
    g.set_ylabels("Dysregulated edges")
    g.tight_layout()

    plt.savefig(output_edges_path, dpi=300)
    plt.savefig(os.path.splitext(output_edges_path)[0] + '.pdf', dpi=300)

    # plot nodes
    sns.set(style="whitegrid")
    g = sns.catplot(
        df, kind="box",
        y="nodes", x="cancer", col="network", col_wrap=2, hue="method", palette="Set1",
        height=3, aspect=2, legend=None, sharey=False
    )
    g.set_titles(col_template="{col_name}", row_template="{row_name}")
    [plt.setp(ax.get_xticklabels(), rotation=45) for ax in g.axes.flat]
    g.add_legend(title="Method")
    g.set_xlabels("Cancer")
    g.set_ylabels("Number of nodes")
    g.tight_layout()

    plt.savefig(output_nodes_path, dpi=300)
    plt.savefig(os.path.splitext(output_nodes_path)[0] + '.pdf', dpi=300)

    # plot number of connected components
    sns.set(style="whitegrid")
    g = sns.catplot(
        df, kind="box",
        y="connected_components", x="cancer", col="network", col_wrap=2, hue="method", palette="Set1",
        height=3, aspect=2, legend=None, sharey=False
    )
    g.set_titles(col_template="{col_name}", row_template="{row_name}")
    [plt.setp(ax.get_xticklabels(), rotation=45) for ax in g.axes.flat]
    g.add_legend(title="Method")
    g.set_xlabels("Cancer")
    g.set_ylabels("Connected components")
    g.tight_layout()

    plt.savefig(output_components_path, dpi=300)
    plt.savefig(os.path.splitext(output_components_path)[0] + '.pdf', dpi=300)

    # plot mean nodes per component
    sns.set(style="whitegrid")
    g = sns.catplot(
        df, kind="box",
        y="mean_nodes_per_component", x="cancer", col="network", col_wrap=2, hue="method", palette="Set1",
        height=3, aspect=2, legend=None, sharey=False
    )
    g.set_titles(col_template="{col_name}", row_template="{row_name}")
    [plt.setp(ax.get_xticklabels(), rotation=45) for ax in g.axes.flat]
    g.add_legend(title="Method")
    g.set_xlabels("Cancer")
    g.set_ylabels("Mean nodes per component")
    g.tight_layout()

    plt.savefig(output_nodes_per_component_path, dpi=300)
    plt.savefig(os.path.splitext(output_nodes_per_component_path)[0] + '.pdf', dpi=300)


if __name__ == "__main__":
    main()
