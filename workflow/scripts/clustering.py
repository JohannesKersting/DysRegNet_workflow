import seaborn as sns
import colorcet as cc
import pandas as pd
import pyarrow
import sys
import argparse
import os
from matplotlib.patches import Patch
import matplotlib.pyplot as plt

def get_clustermap(overlaps, meta):

    # styling
    cmap = sns.cubehelix_palette(start=.5, rot=-.5, as_cmap=True)
    sns.set(font_scale=0.8)

    # color mapping for legend
    cancer_types = list(set(meta.loc[overlaps.index]._primary_disease))  # list of unique classes
    cancer_types.sort()  # sort for reproducible color mapping
    colors = sns.color_palette(cc.glasbey_light, n_colors=len(cancer_types) + 1)[1:]  # list of colors
    lut = {cancer_type: color for cancer_type, color in zip(cancer_types, colors)}  # map classes to colors

    col_colors = meta.loc[overlaps.index]._primary_disease.map(lut)
    col_colors = col_colors.rename("Cancer")

    g = sns.clustermap(overlaps, row_colors=col_colors, row_cluster=True, col_cluster=True,
                       col_colors=col_colors, figsize=(8, 8), cmap=cmap, yticklabels=False, xticklabels=False)
    ax = g.ax_heatmap
    ax.set_ylabel("Patients")
    ax.set_xlabel("Patients")
    g.ax_col_dendrogram.set_visible(False)

    return g, lut


def add_legend(g, lut):
    sns.set(font_scale=0.8)

    handles = [Patch(facecolor=lut[name]) for name in lut]
    ax = g.ax_heatmap
    ax.legend(handles, lut, title="Cancer", ncol=2, bbox_to_anchor=(0.5, 1.32), loc="upper center", borderaxespad=0)

    return g


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--overlaps_edges',  type=str, required=True,
                        help='Paths to edge overlaps')
    parser.add_argument('--overlaps_nodes',  type=str, required=True,
                        help='Paths to node overlaps')
    parser.add_argument('--meta',  type=str, required=True,
                        help='Paths to multiple patient specific networks')
    parser.add_argument('--output_edges', type=str, required=True,
                        help='Path to the edge cluster heatmap')
    parser.add_argument('--output_nodes', type=str, required=True,
                        help='Path to the node cluster heatmap')
    args = parser.parse_args()

    print(f"Python version:\n{sys.version}")
    print(f"Arguments:\n{args}\n")


    # read meta sheet
    meta = pd.read_csv(args.meta, index_col=0, sep = "\t")
    print(f"meta: \n{meta}\n")


    # create edge heatmap
    overlaps = pd.read_feather(args.overlaps_edges)
    overlaps = overlaps.set_index("patient id")
    print(f"overlaps edges: \n{overlaps}\n")

    g, lut = get_clustermap(overlaps, meta)
    g.savefig(args.output_edges, dpi=350)

    g = add_legend(g, lut)
    g.savefig(os.path.splitext(args.output_edges)[0] + '_legend.png', dpi=350, bbox_inches="tight")


    # create nodes heatmap
    overlaps = pd.read_feather(args.overlaps_nodes)
    overlaps = overlaps.set_index("patient id")
    print(f"overlaps nodes: \n{overlaps}\n")

    g, lut = get_clustermap(overlaps, meta)
    g.savefig(args.output_nodes, dpi=350)

    g = add_legend(g, lut)
    g.savefig(os.path.splitext(args.output_nodes)[0] + '_legend.png', dpi=350, bbox_inches="tight")


if __name__ == "__main__":
    main()
