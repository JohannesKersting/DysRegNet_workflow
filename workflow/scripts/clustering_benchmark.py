import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from munkres import Munkres
import pyarrow
from sklearn.cluster import SpectralClustering
from sklearn.metrics.cluster import contingency_matrix
from sklearn.preprocessing import normalize
import argparse
import os
import sys

from commons import method_lut, network_lut


def get_input_info(path):
    '''
    Returns a tuple for each input path to extract information from it.
    The order is path, method, overlap_type, norm_method, ref_network_start.
    '''

    # split path by separator
    path = os.path.normpath(path)
    splitted_path = path.split(os.sep)

    # get method and filename
    file_name = splitted_path[-1]
    method = splitted_path[-2]

    # get overlap_type, norm method and the reference network
    stem_name, file_extension = os.path.splitext(file_name)
    overlap_type, norm_method, ref_network = stem_name.split("-")

    # shorten refernece network name (for genie3)
    ref_network_start = ref_network.split("_")[0]

    return path, method, overlap_type, norm_method, ref_network_start


def f_matrix(labels_pred, labels_true):
    # Calculate F1 matrix
    cont_mat = contingency_matrix(labels_pred=labels_pred, labels_true=labels_true)
    precision = normalize(cont_mat, norm='l1', axis=0)
    recall = normalize(cont_mat, norm='l1', axis=1)
    som = precision + recall
    f1 = np.round(np.divide((2 * recall * precision), som, out=np.zeros_like(som), where=som != 0), 3)
    return f1


def f1_hungarian(f1):
    m = Munkres()
    inverse = 1 - f1
    indices = m.compute(inverse.tolist())
    fscore = sum([f1[i] for i in indices]) / len(indices)
    return fscore


def clustering_benchmark(overlaps_path, meta):
    """

    :param overlaps_path: A path to an overlap file in feather format
    :param meta: The meta sheet as a pandas df
    :return: The clustering F1 score
    """
    # read overlaps file
    overlaps = pd.read_feather(overlaps_path)
    overlaps = overlaps.set_index("patient id")

    # select relevant meta data
    patient_list = list(overlaps.index)
    meta_temp = meta.loc[patient_list]

    n_clusters = len(set(meta_temp._primary_disease))
    clustering = SpectralClustering(affinity="precomputed", n_clusters=n_clusters, random_state=0).fit(overlaps)

    # get class names
    class_names = meta_temp._primary_disease.values
    cs = list(set(class_names))
    cs.sort()

    # indices of true lables
    true = [cs.index(x) for x in class_names]

    f_mx = f_matrix(true, clustering.labels_)
    hungarian = f1_hungarian(f_mx)
    return hungarian


def main():

    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', nargs='+', required=True,
                        help='Multiple overlap files in feather format')
    parser.add_argument('--meta',  type=str, required=True,
                        help='Paths to multiple patient specific networks')
    parser.add_argument('--output', type=str, required=True,
                        help='Output for the clustering benchmark plot')

    args = parser.parse_args()

    print(f"Python version:\n{sys.version}")
    print(f"Arguments:\n{args}\n")

    input_paths = args.input
    meta_path = args.meta
    output_path = args.output


    # read meta sheet
    meta = pd.read_csv(meta_path, index_col=0, sep = "\t")
    print(f"meta: \n{meta}\n")


    # get infos from input paths

    input_infos = [get_input_info(path) for path in input_paths]

    input_df = pd.DataFrame(input_infos, columns=["path", "method", "overlap_type", "norm_method", "network"])
    input_df["method"] = input_df["method"].map(lambda x: method_lut[x] if x in method_lut.keys() else x)
    input_df["network"] = input_df["network"].map(lambda x: network_lut[x] if x in network_lut.keys() else x)
    print(f"input_df: \n{input_df.to_string()}\n")


    # cluster benchmarking
    input_df["F1"] = [clustering_benchmark(path, meta) for path in input_df["path"].values]
    print(f"input_df with F1: \n{input_df.to_string()}\n")


    # plotting
    plt.figure(figsize=(2*len(set(input_df["method"])),3))
    sns.set(font_scale=0.8, style="whitegrid")

    g = sns.barplot(x="method", hue="network", y="F1", data=input_df, palette="Set2")
    g.set(ylim=(0, 1), xlabel='Method')
    g.legend(title="Reference network", bbox_to_anchor=(1, 0), loc="lower right", borderaxespad=1)
    plt.tight_layout() # else saving will cut of the axis labels

    plt.savefig(output_path, dpi=600)
    plt.savefig(os.path.splitext(output_path)[0] + '.pdf', dpi=600)
    plt.savefig(os.path.splitext(output_path)[0] + '.eps', dpi=600)



if __name__ == "__main__":
    main()