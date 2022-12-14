
import argparse
import pandas as pd
import sys


parser = argparse.ArgumentParser()

parser.add_argument('--input', nargs='+', help='All input files (separated with whitespaces)', required=True)
parser.add_argument('--output', type=str, help="Output file for the summed up network", required=True)

args = parser.parse_args()

print("Python version:")
print(sys.version)
print("Arguments:")
print(args)
print("\n")

input_files = args.input
output_file = args.output

# put all edges from genie3 in one dictionary and add their weights
edges = dict()
for path in input_files:
    print(f"\n{path}")
    normal_net = pd.read_csv(path)
    print(f"Number of edges: {normal_net.shape[0]}")
    normal_net = normal_net[normal_net.weight > 0]
    for tup in normal_net.itertuples():
        edge = (tup[1], tup[2])
        if edge in edges:
            edges[edge].append(tup[3])
        else:
            edges[edge] = [tup[3]]

# keep only edges present in all networks and some them up
edges_small = dict()
for edge in edges:
    if len(edges[edge]) == len(input_files):
        edges_small[edge] = sum(edges[edge])

df = pd.DataFrame([(k[0], k[1], v) for k, v in edges_small.items()], columns =["regulatoryGene", "targetGene", "weight"])
df = df.sort_values(by=['weight'], ascending=False)
print(f"Saving a shared network with {df.shape[0]} edges (sorted by weight)...")
df.to_csv(output_file, index=False)