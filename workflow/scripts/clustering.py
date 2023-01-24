import seaborn as sns
import colorcet as cc
import pandas as pd
import pyarrow
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--overlaps',  type=str, required=True,
                    help='Paths to multiple patient specific networks')
parser.add_argument('--meta',  type=str, required=True,
                    help='Paths to multiple patient specific networks')
parser.add_argument('--output', type=str, required=True,
                    help='Path to the output file (will be stored in feather format)')

args = parser.parse_args()

print("Python version:")
print(sys.version)
print("Arguments:")
print(args)
print("\n")

meta = pd.read_csv(args.meta, index_col=0, sep = "\t")
print(meta)

overlaps = pd.read_feather(args.overlaps)
overlaps = overlaps.set_index("patient id")
print(overlaps)


colours = sns.color_palette(cc.glasbey_light, n_colors=12)[1:]

classes = list(set(meta.loc[overlaps.index]._primary_disease))



lut = {classes[i]: colours[i] for i in range(len(classes))}
print(lut)

col_colors = meta.loc[overlaps.index]._primary_disease.map(lut)
col_colors = col_colors.rename("cancer")

print(col_colors)

sns.set(font_scale=0.8)
cmap = sns.cubehelix_palette(start=.5, rot=-.5, as_cmap=True)

g = sns.clustermap(overlaps, row_colors=col_colors, row_cluster=True, col_cluster=True,
                   col_colors=col_colors, figsize=(5, 5), cmap= cmap, yticklabels = False, xticklabels = False)
ax = g.ax_heatmap
ax.set_xlabel("Patients")
ax.set_ylabel("Patients")
g.ax_col_dendrogram.set_visible(False)
g.savefig(args.output, dpi = 350)
# g.savefig("figures/heatmap_edges.pdf")

