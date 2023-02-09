import pandas as pd
import numpy as np

from tqdm import tqdm
from functools import partialmethod
tqdm.__init__ = partialmethod(tqdm.__init__, disable=True)

#import pickle
#import json
import pyarrow
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--input', nargs='+', required=True,
                    help='Paths to multiple patient specific networks')
parser.add_argument('--output_edges', type=str, required=True,
                    help='Path to the edge output file (will be stored in feather format)')
parser.add_argument('--output_nodes', type=str, required=True,
                    help='Path to the node output file (will be stored in feather format)')

args = parser.parse_args()

print("Python version:")
print(sys.version)
print("Arguments:")
print(args)
print("\n")


nodes = []
edges = []
patient_list = []

index_col = "patient id"

# loop over cancer types
for input_path in args.input:
    print(f"Reading {input_path}")
    net = pd.read_feather(input_path)
    net = net.set_index(index_col)

    crc = net.values
    all_edges = list(net.columns)
    all_patients = list(net.index)

    # loop over patients
    for i in tqdm(range(crc.shape[0])):
        patient_list.append(all_patients[i])
        count = 0
        node_set = set()
        edge_list = []
        # loop over edges
        for edge in all_edges:

            # get nodes from edge names
            v1, v2 = edge.split(", ")
            v1 = v1[2:-1]
            v2 = v2[1:-2]

            if abs(crc[i, count]) > 0:
                node_set.add(v1)
                node_set.add(v2)
                edge_list.append((v1, v2))

            count += 1
        nodes.append(node_set)
        edges.append(edge_list)

#with open(f'benchmark/edges_{out}.pkl', 'wb') as f:
#    pickle.dump(edges, f)

def get_overlap(x, y):
    return (len(x.intersection(y)) / min(len(x), len(y)))

def overlap_mat(patient_list, networks):

    overlap_mat = np.zeros((len(patient_list), len(patient_list)))

    idx = [(i, j) for i in range(len(patient_list)) for j in range(i, len(patient_list))]
    for idx in tqdm(idx):
        i, j = idx
        set_1 = networks[i]
        set_2 = networks[j]
        try:
            overlap = get_overlap(set(set_1), set(set_2))
        except ZeroDivisionError:
            overlap = 0

        overlap_mat[i, j] = overlap
        overlap_mat[j, i] = overlap_mat[i, j]

    return pd.DataFrame(overlap_mat, index=patient_list, columns=patient_list)


overlap_edges = overlap_mat(patient_list, edges).reset_index(names=index_col)
print("Edge overlaps:")
print(overlap_edges)
overlap_edges.to_feather(args.output_edges)

overlap_nodes = overlap_mat(patient_list, nodes).reset_index(names=index_col)
print("Node overlaps:")
print(overlap_nodes)
overlap_nodes.to_feather(args.output_nodes)