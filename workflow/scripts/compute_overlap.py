import pandas as pd
import numpy as np

from tqdm import tqdm
from functools import partialmethod

tqdm.__init__ = partialmethod(tqdm.__init__, disable=True)

# import pickle
# import json
import pyarrow
import sys
import argparse


def get_sets(input_list, index_col, only_shared_edges, sign):
    """

    :param input_list: list with paths to patient-specific networks
    :param index_col: name of colum, which should be used as index
    :param only_shared_edges: filter all edges and nodes, which do not appear in all input networks
    :param sign: 0 -> look at all edges wich are not 0, sign > 0 -> look only at positive edges, sigh < 0 -> look only at negativ edges
    :return: (patient_list, edges, nodes)
    """
    # list of node sets (one set for every patient)
    nodes = []

    # list of edge lists (one list for every patient)
    edges = []

    # stores edges present in all networks
    shared_edges = None

    # list of all patients in all networks
    patient_list = []

    # loop over cancer types
    for input_path in input_list:

        print(f"Reading {input_path}")
        net = pd.read_feather(input_path)
        net = net.set_index(index_col)

        # values as matrix
        crc = net.values

        # parse index to list of node name tuples
        all_edges = [edge.split(", ") for edge in list(net.columns)]
        all_edges = [(v1[2:-1], v2[1:-2]) for v1, v2 in all_edges]
        all_nodes = {node for edge in all_edges for node in edge}

        all_patients = list(net.index)

        print(f"Patients: {len(all_patients)}")
        print(f"Edges: {len(all_edges)}")
        print(f"Nodes: {len(all_nodes)}")

        # update shared_edges
        if only_shared_edges:
            if shared_edges is None:
                shared_edges = set(all_edges)
            else:
                shared_edges = shared_edges.intersection(set(all_edges))

        # loop over patients
        for i in tqdm(range(crc.shape[0])):
            patient_list.append(all_patients[i])
            count = 0
            node_set = set()
            edge_list = []

            # loop over edges
            for v1, v2 in all_edges:
                if sign == 0 and abs(crc[i, count]) > 0:
                    node_set.add(v1)
                    node_set.add(v2)
                    edge_list.append((v1, v2))
                elif sign > 0 and crc[i, count] > 0:
                    node_set.add(v1)
                    node_set.add(v2)
                    edge_list.append((v1, v2))
                elif sign < 0 and crc[i, count] < 0:
                    node_set.add(v1)
                    node_set.add(v2)
                    edge_list.append((v1, v2))
                count += 1

            nodes.append(node_set)
            edges.append(edge_list)

    print(f"\nTotal number of dysregulated edges in all networks and patients:", end=' ')
    print(sum([len(patient_edges) for patient_edges in edges]))
    print(f"Total number of dysregulated nodes in all networks and patients:", end=' ')
    print(sum([len(patient_nodes) for patient_nodes in nodes]))
    print(f"Total number of patients: {len(patient_list)}")

    if only_shared_edges:
        print("\nFiltering for shared edges and nodes...")

        shared_nodes = {node for edge in shared_edges for node in edge}

        print(f"Shared edges: {len(shared_edges)}")
        print(f"Shared nodes: {len(shared_nodes)}")

        # filter edges
        filtered_edges = []
        for patient_edges in edges:
            filtered_patient_edges = [edge for edge in patient_edges if edge in shared_edges]
            filtered_edges.append(filtered_patient_edges)
        edges = filtered_edges

        # filter nodes
        filtered_nodes = []
        for patient_nodes in nodes:
            filtered_patient_nodes = {node for node in patient_nodes if node in shared_nodes}
            filtered_nodes.append(filtered_patient_nodes)
        nodes = filtered_nodes

        print(f"Filtered number of dysregulated edges in all networks and patients:", end=' ')
        print(sum([len(patient_edges) for patient_edges in edges]))
        print(f"Filtered number of dysregulated nodes in all networks and patients:", end=' ')
        print(sum([len(patient_nodes) for patient_nodes in nodes]))

    return patient_list, edges, nodes

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

def combined_overlaps(patient_list, index_col, pos, neg):
    overlap_pos = overlap_mat(patient_list, pos)
    overlap_neg = overlap_mat(patient_list, neg)
    overlap = (overlap_pos + overlap_neg) / 2
    overlap = overlap.reset_index(names=index_col)
    return overlap


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--input', nargs='+', required=True,
                        help='Paths to multiple patient specific networks')
    parser.add_argument("--shared", action="store_true",
                        help="Consider only edges which appear in all input networks")
    parser.add_argument("--signed", action="store_true",
                        help="The overlap calculation should distinguish between positive and negative values in the input network")
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

    input_list = args.input
    output_edges = args.output_edges
    output_nodes = args.output_nodes
    only_shared_edges = args.shared
    signed = args.signed

    index_col = "patient id"

    if not signed:

        patient_list, edges, nodes = get_sets(input_list, index_col=index_col, only_shared_edges=only_shared_edges, sign=0)

        overlap_edges = overlap_mat(patient_list, edges).reset_index(names=index_col)
        print("Edge overlaps:")
        print(overlap_edges)
        overlap_edges.to_feather(output_edges)
        overlap_edges = None
        edges = None

        overlap_nodes = overlap_mat(patient_list, nodes).reset_index(names=index_col)
        print("Node overlaps:")
        print(overlap_nodes)
        overlap_nodes.to_feather(output_nodes)

        overlap_nodes = None
        nodes = None

    else:
        print("\nReading positive values:")
        patient_list_pos, edges_pos, nodes_pos = get_sets(input_list, index_col=index_col, only_shared_edges=only_shared_edges, sign=1)

        print("\nReading negative values:")
        patient_list_neg, edges_neg, nodes_neg = get_sets(input_list, index_col=index_col, only_shared_edges=only_shared_edges, sign=-1)

        assert patient_list_pos == patient_list_neg

        overlap_edges = combined_overlaps(patient_list_pos, index_col, edges_pos, edges_neg)
        edges_pos = 0
        edges_neg = 0

        print("Edge overlaps:")
        print(overlap_edges)
        overlap_edges.to_feather(output_edges)
        overlap_edges = None


        overlap_nodes = combined_overlaps(patient_list_pos, index_col, nodes_pos, nodes_neg)
        nodes_pos = None
        nodes_neg = None

        print("Node overlaps:")
        print(overlap_nodes)
        overlap_nodes.to_feather(output_nodes)
        overlap_nodes= None


if __name__ == "__main__":
    main()
