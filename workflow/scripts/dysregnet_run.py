# Silence tqdm to keep log file clean
from tqdm import tqdm
from functools import partialmethod

tqdm.__init__ = partialmethod(tqdm.__init__, disable=True)


import pandas as pd
import numpy as np
import dysregnet
import argparse
import pyarrow
import sys




parser = argparse.ArgumentParser()

parser.add_argument('--expr', type=str, required=True,
                    help='Path to the expression CSV matrix (genes as rows, samples as columns)')
parser.add_argument('--meta', type=str, required=True,
                    help='Path to the CSV meta file. Sample ids must be in "sample" column'
                    )
parser.add_argument('--grn', type=str, required=True,
                    help='Path to a reference network in CSV format. Must have exactly two columns and a header.'
                    )
parser.add_argument('--output', type=str, required=True,
                    help='The output file path and name.'
                    )

args = parser.parse_args()


print("Python version:")
print(sys.version)
print("Arguments:")
print(args)
print("\n")

conCov = ["birth_days_to"]
catCov = ["gender", "race"]
conCol = "condition"


# Read data
meta = pd.read_csv(args.meta)
expr = pd.read_csv(args.expr)
grn = pd.read_csv(args.grn)

# Pre-process data

meta = meta[ ["sample",conCol] + conCov + catCov ]
meta["birth_days_to"] = meta["birth_days_to"].fillna(meta["birth_days_to"].mean())
meta['race'] = meta['race'].fillna('not reported')
meta["race"] = meta["race"].replace({"[Unknown]": 'not reported', "[Not Evaluated]":'not reported'})

expr = expr.set_index(expr.columns[0])
if not all(expr.columns == meta.iloc[:,0]):
    raise ValueError("Column names of expression differ from 'sample' column in meta")
expr = expr.T
expr.insert(0, 'sample_ids', expr.index)

grn = grn.iloc[:, 0:2]


# Print data summary
print("\nData as submitted to dysregnet:")
print("\nexpr:")
print(expr)

print("\nmeta:")
print(meta)

print("\ngrn:")
print(grn)

print("\nMeta composition:\n")
print(meta.groupby(conCol).size())
print("\n")

for covar in catCov:
    print(pd.crosstab(meta[covar],meta[conCol]))
    print("")

# Run dysregnet
data=dysregnet.run(expression_data=expr,
                   meta=meta,
                   CatCov=catCov,
                   ConCov=conCov,
                   GRN=grn,
                   conCol=conCol)

result = data.get_results()

# feather needs str column names
result.columns = [str(x) for x in result.columns]

# Print some result statistics
print("Result:")
print(result)

# This works only with applied condition criterion for dysregnet
col_sums = result[result.columns[1:]].apply(np.sum, axis=0)

print(f"Number of edges with at least one dysregulation: {np.sum(col_sums!=0)}")
print(f"Number of positive slopes: {np.sum(col_sums>0)}")
print(f"Number of negative slopes: {np.sum(col_sums<0)}")


# Write output
result.to_feather(args.output)
