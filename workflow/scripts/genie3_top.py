
import argparse
import pandas as pd
import sys


parser = argparse.ArgumentParser()

parser.add_argument('--input', type=str, help='Input genie3 network', required=True)
parser.add_argument('--top_k', type=int, help="Number of edges (in thousands) that should be in the output")
parser.add_argument('--output', type=str, help="Output file for top edges", required=True)

args = parser.parse_args()

print("Python version:")
print(sys.version)
print("Arguments:")
print(args)
print("\n")

input_file = args.input
output_file = args.output
top = args.top_k * 1000

normal_net = pd.read_csv(input_file)
normal_net = normal_net[normal_net.weight > 0]
normal_net = normal_net.sort_values(by=['weight'], ascending=False)
normal_net = normal_net.head(top)

normal_net.to_csv(output_file, index=False)

