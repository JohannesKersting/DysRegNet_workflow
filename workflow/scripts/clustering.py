import seaborn as sns
import colorcet as cc
import pandas as pd
import pyarrow
import sys
import argparse
import os
from matplotlib.patches import Patch
import matplotlib.pyplot as plt

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

print("Python version:")
print(sys.version)
print("Arguments:")
print(args)
print("\n")

# read meta sheet
meta = pd.read_csv(args.meta, index_col=0, sep = "\t")
print(meta)

# styling
sns.set(font_scale=0.8)
cmap = sns.cubehelix_palette(start=.5, rot=-.5, as_cmap=True)


# create edge heatmap
overlaps = pd.read_feather(args.overlaps_edges)
overlaps = overlaps.set_index("patient id")
print(overlaps)

colours = sns.color_palette(cc.glasbey_light, n_colors=12)[1:] # list of colors
classes = list(set(meta.loc[overlaps.index]._primary_disease)) # list of classes
classes.sort() # sort for reproducible color mapping

lut = {classes[i]: colours[i] for i in range(len(classes))} # map classes to colors

col_colors = meta.loc[overlaps.index]._primary_disease.map(lut)
col_colors = col_colors.rename("cancer")

g = sns.clustermap(overlaps, row_colors=col_colors, row_cluster=True, col_cluster=True,
                   col_colors=col_colors, figsize=(5, 5), cmap= cmap, yticklabels = False, xticklabels = False)
ax = g.ax_heatmap
ax.set_xlabel("Patients")
ax.set_ylabel("Patients")
g.ax_col_dendrogram.set_visible(False)
g.savefig(args.output_edges, dpi = 350)


handles = [Patch(facecolor=lut[name]) for name in lut]
plt.legend(handles, lut, title='Cancer',
           bbox_to_anchor=(0.2, 1), ncol = 2, bbox_transform=plt.gcf().transFigure, loc='upper center')
plt.subplots_adjust(top=0.3)
g.savefig(os.path.splitext(args.output_edges)[0] + '_legend.png', dpi = 350, bbox_inches="tight")



# create node heatmap
overlaps = pd.read_feather(args.overlaps_nodes)
overlaps = overlaps.set_index("patient id")
print(overlaps)

# prepare colors
colours = sns.color_palette(cc.glasbey_light, n_colors=12)[1:] # list of colors
classes = list(set(meta.loc[overlaps.index]._primary_disease)) # list of classes
classes.sort() # sort for reproducible color mapping

lut = {classes[i]: colours[i] for i in range(len(classes))} # map classes to colors

col_colors = meta.loc[overlaps.index]._primary_disease.map(lut)
col_colors = col_colors.rename("cancer")

g = sns.clustermap(overlaps, row_colors=col_colors, row_cluster=True, col_cluster=True,
                   col_colors=col_colors, figsize=(5, 5), cmap= cmap, yticklabels = False, xticklabels = False)
ax = g.ax_heatmap
ax.set_xlabel("Patients")
ax.set_ylabel("Patients")
g.ax_col_dendrogram.set_visible(False)
g.savefig(args.output_nodes, dpi = 350)

handles = [Patch(facecolor=lut[name]) for name in lut]
plt.legend(handles, lut, title='Cancer',
           bbox_to_anchor=(1, 1), ncol = 2, bbox_transform=plt.gcf().transFigure, loc='upper center')
#plt.subplots_adjust(top=0.75)

g.savefig(os.path.splitext(args.output_nodes)[0] + '_legend.png', dpi = 350, bbox_inches="tight")