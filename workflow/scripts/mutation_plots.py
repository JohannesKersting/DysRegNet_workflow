import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
import sys
import colorcet as cc


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

    args = parser.parse_args()

    print(f"Python version:\n{sys.version}")
    print(f"Arguments:\n{args}\n")

    input_paths_local = args.inputs_local
    input_paths_global = args.inputs_global
    output_path = args.output

    # map input files
    method_lut = {"dysregnet": "DysRegNet",
                  "ssn": "SSN",
                  }
    network_lut = {"exp": "experimental",
                   "string": "string",
                   "genie3_shared": "genie3 shared",
                   "genie3_individual": "genie3 individual",
                   }

    input_infos = [get_input_info(path) for path in input_paths_local]

    input_df = pd.DataFrame(input_infos, columns=["local_path", "method", "norm_method", "network"])
    input_df["method"] = input_df["method"].map(method_lut)
    input_df["network"] = input_df["network"].map(network_lut)
    input_df.insert(1, "global_path", input_paths_global)

    print(f"input_df: \n{input_df.to_string()}\n")


    # parse input data
    mutation_test_data = []
    for index, row in input_df.iterrows():
        global_tests_df = pd.read_csv(row["global_path"])
        global_tests_df.sort_values(by=['cancer'], inplace=True)
        global_tests_df["network"] = row["network"]
        mutation_test_data.append(global_tests_df)

    mutation_test_df = pd.concat(mutation_test_data, ignore_index=True)
    mutation_test_data = None

    mutation_test_df["% of significant associations"] = 100 * mutation_test_df["n_sign_fdr"] / mutation_test_df["n_tfs"]
    mutation_test_df["Globally sig."] = np.where(mutation_test_df["global_pval"] < 0.05, 'yes', 'no')
    mutation_test_df = mutation_test_df.rename(columns={"cancer": "Cancer"})

    print(f"mutation_test_df: \n{mutation_test_df.to_string()}\n")


    # plot
    sns.set(style="whitegrid")
    palette = sns.color_palette(cc.glasbey_light, n_colors=len(set(mutation_test_df["Cancer"])) + 1)[1:]

    g = sns.relplot(data=mutation_test_df, x="n_tfs", y="% of significant associations", hue="Cancer",
                    style="Globally sig.", col="network",
                    s=90, style_order=["no", "yes"], height=3, aspect=1.6, col_wrap=2, palette=palette)
    g.set_axis_labels("Number of tests")
    g.set_titles(col_template="{col_name}")
    g.tight_layout()

    plt.savefig(output_path, dpi=300)
    plt.savefig(os.path.splitext(output_path)[0] + '.pdf', dpi=300)


if __name__ == "__main__":
    main()