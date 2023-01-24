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
parser.add_argument('--output', type=str, required=True,
                    help='Path to the output file (will be stored in feather format)')

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

for input_path in args.input:
    print(f"Reading {input_path}")
    net = pd.read_feather(input_path)
    net = net.set_index(index_col)

    crc = net.values
    all_edges = list(net.columns)
    all_patients = list(net.index)

    for i in tqdm(range(crc.shape[0])):
        patient_list.append(all_patients[i])
        count = 0
        edge_list_pos = []
        edge_list_neg = []

        nodes_small = []
        for edge in all_edges:
            v1, v2 = edge.split(", ")
            v1 = v1[2:-1]
            v2 = v2[1:-2]
            if crc[i, count] > 0:
                edge_list_pos.append((v1, v2))
                nodes_small.append(v1)
                nodes_small.append(v2)
            elif crc[i, count] < 0:
                edge_list_neg.append((v1, v2))
                nodes_small.append(v1)
                nodes_small.append(v2)
            count += 1
        nodes.append(nodes_small)
        edges.append((edge_list_pos, edge_list_neg))

#with open(f'benchmark/edges_{out}.pkl', 'wb') as f:
#    pickle.dump(edges, f)


def overlap(x, y):
    return (len(x.intersection(y)) / min(len(x), len(y)))


def overlap_mat(patient_list, networks):

    overlap_mat = np.zeros((len(patient_list), len(patient_list)))
    idx = [(i, j) for i in range(len(patient_list)) for j in range(i, len(patient_list))]
    for idx in tqdm(idx):
        i, j = idx
        s1 = networks[i]
        s2 = networks[j]
        try:
            overlap_mat[i, j] = overlap(s1, s2)
        except ZeroDivisionError:
            print(i, j)
        overlap_mat[j, i] = overlap_mat[i, j]
    return (pd.DataFrame(overlap_mat, index=patient_list, columns=patient_list))


def overlap_edges(patient_list, networks):
    overlap_mat_pos = np.zeros((len(patient_list), len(patient_list)))
    idx = [(i, j) for i in range(len(patient_list)) for j in range(i, len(patient_list))]
    for idx in tqdm(idx):
        i, j = idx
        edges_pos1, edges_neg1 = networks[i]
        edges_pos2, edges_neg2 = networks[j]
        try:
            o_pos = overlap(set(edges_pos1), set(edges_pos2))
        except ZeroDivisionError:
            o_pos = 0

        overlap_mat_pos[i, j] = o_pos
        overlap_mat_pos[j, i] = overlap_mat_pos[i, j]
    overlap_mat_pos = pd.DataFrame(overlap_mat_pos, index=patient_list, columns=patient_list)

    return (overlap_mat_pos)


overlap_pos = overlap_edges(patient_list, edges).reset_index(names=index_col)
print(overlap_pos)
overlap_pos.to_feather(args.output)