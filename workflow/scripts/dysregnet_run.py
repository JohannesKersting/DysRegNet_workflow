import pandas as pd
import dysregnet
import argparse
import pyarrow
import sys



parser = argparse.ArgumentParser()

parser.add_argument('--expr', type=str, required=True)
parser.add_argument('--meta', type=str, required=True)
parser.add_argument('--grn', type=str, required=True)
parser.add_argument('--output', type=str, required=True)

args = parser.parse_args()

print("Python version:")
print(sys.version)
print("Arguments:")
print(args)
print("\n")

df = pd.DataFrame({'Temperature': ['Hot', 'Cold', 'Warm', 'Cold'],
                   })
df.to_feather(args.output)

