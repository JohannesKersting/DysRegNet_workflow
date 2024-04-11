import os
import pandas as pd
import pyarrow
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

    # get overlap_type, norm method and the reference network
    stem_name, file_extension = os.path.splitext(file_name)
    norm_method, ref_network = stem_name.split("-")

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
                        help='.fea files with dysregulated networks')
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
    input_df["network"] = input_df["network"].map(
        lambda x: "genie3_individual" if (("genie3" in x) and (not "shared" in x)) else x)
    input_df["network"] = input_df["network"].map(lambda x: network_lut[x] if x in network_lut.keys() else x)
    print(input_df.to_string())


    # read files
    summary_list = []
    for index, row in input_df.iterrows():
        net = pd.read_feather(row["path"])
        net = net.set_index("patient id")
        for patient_index, patient_row in net.iterrows():
            summary_list.append({
                "cancer": row["cancer"],
                "method": row["method"],
                "patient": patient_index,
                "norm_method": row["norm_method"],
                "network": row["network"],
                "dysregulated_edges": np.sum(patient_row != 0),
            })
    summary_df = pd.DataFrame.from_dict(summary_list)
    summary_list = None

    print(summary_df)


    # plot
    sns.set(style="whitegrid")
    g = sns.catplot(
        summary_df, kind="box",
        y="dysregulated_edges", x="cancer", col="network", col_wrap=2, hue="method", palette="Set1",
        height=3, aspect=2, legend=None, sharey=False
    )
    g.set_titles(col_template="{col_name}", row_template="{row_name}")
    [plt.setp(ax.get_xticklabels(), rotation=45) for ax in g.axes.flat]
    g.add_legend(title="Method")
    g.set_xlabels("Cancer")
    g.set_ylabels("Dysregulated edges")
    g.tight_layout()

    plt.savefig(output_path, dpi=300)
    plt.savefig(os.path.splitext(output_path)[0] + '.pdf', dpi=300)


if __name__ == "__main__":
    main()
