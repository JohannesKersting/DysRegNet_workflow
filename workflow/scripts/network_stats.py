import glob
import os
import sys
import pandas as pd
import pyarrow
import numpy as np
import networkx as nx
from ast import literal_eval as make_tuple
import argparse

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, required=True,
                        help='Paths to a patient-specific networks file (fea)')
    parser.add_argument('--output', type=str, required=True,
                        help='Path to output file for network statistics (csv)')
    args = parser.parse_args()

    print(f"Python version:\n{sys.version}")
    print(f"Arguments:\n{args}\n")

    input_path = args.input
    output_path = args.output

    # read and format data
    net = pd.read_feather(input_path)
    net = net.set_index("patient id")
    net.columns = [make_tuple(col) for col in net.columns]

    # iterate over patients
    summary_list = []
    for patient_index, patient_row in net.iterrows():

        # select only dysregulated edges
        edge_list = list(patient_row[patient_row != 0].index)

        # convert to networkx
        G = nx.Graph(edge_list)

        # collect stats
        summary_list.append({
            "patient": patient_index,
            "nodes": G.number_of_nodes(),
            "edges": len(edge_list),
            "connected_components": nx.number_connected_components(G),
        })

    # combine into one df and save
    summary_df = pd.DataFrame.from_dict(summary_list)
    summary_list = None
    summary_df.to_csv(output_path, index=False)


if __name__ == "__main__":
    main()