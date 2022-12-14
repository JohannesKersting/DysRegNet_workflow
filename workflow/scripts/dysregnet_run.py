import pandas as pd
import dysregnet
import argparse
import sys



parser = argparse.ArgumentParser()

parser.add_argument('--expr', type=str, required=True)
parser.add_argument('--meta', type=str, required=True)
parser.add_argument('--grn', type=str, required=True)

args = parser.parse_args()

print("Python version:")
print(sys.version)
print("Arguments:")
print(args)
print("\n")

print("test")

