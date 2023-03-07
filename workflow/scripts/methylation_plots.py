import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pyarrow
import argparse
import os
import sys



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


    # map input paths
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
                cancer_type,
                n_signif,
                n_tests,
                global_pval
            )
            methylation_test_data.append(result)

    methylation_test_df = pd.DataFrame(methylation_test_data,
                                       columns=["network", "cancer", "number significant", "number of tests",
                                                "global p-value"]
                                       )
    methylation_test_data = None
    methylation_test_df["% of significant associations"] = 100 * methylation_test_df["number significant"] / \
                                                           methylation_test_df["number of tests"]

    print(f"methylation_test_df: \n{methylation_test_df.to_string()}\n")


    # plot
    sns.set(style="whitegrid")
    plt.figure(figsize=(4, 8))

    ax = sns.barplot(x="% of significant associations", hue="network", y="cancer", data=methylation_test_df,
                     palette="Set2")

    ax.legend(title="Reference network", bbox_to_anchor=(0, 1.03), loc="lower left", borderaxespad=0)
    plt.tight_layout()
    plt.ylabel("Cancer")

    for index, row in methylation_test_df.iterrows():
        p = ax.patches[index]
        if row["global p-value"] <= 0.5:
            ax.annotate("*",
                        (p.get_x() + p.get_width() + 0.3, p.get_y() + 0.25),
                        xytext=(0, 0),
                        textcoords='offset points'
                        )

    plt.savefig(output_path, dpi=300)
    plt.savefig(os.path.splitext(output_path)[0] + '.pdf', dpi=300)



if __name__ == "__main__":
    main()