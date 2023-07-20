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
                    help='Path to the CSV meta file.'
                    )
parser.add_argument('--grn', type=str, required=True,
                    help='Path to a reference network in CSV format. The first two columns will be selected and TF and TG cols.'
                    )
parser.add_argument('--output', type=str, required=True,
                    help='The output file path (feather).'
                    )
parser.add_argument('--output_stats', type=str, required=True,
                    help='The stats output file path (csv).'
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
idCol = "sample"


# Read data
meta = pd.read_csv(args.meta)
expr = pd.read_csv(args.expr)
grn = pd.read_csv(args.grn)

# Pre-process data

meta = meta[ [idCol, conCol] + conCov + catCov ]
meta["birth_days_to"] = meta["birth_days_to"].fillna(meta["birth_days_to"].mean())
meta['race'] = meta['race'].fillna('not reported')
meta["race"] = meta["race"].replace({"[Unknown]": 'not reported', "[Not Evaluated]":'not reported'})

expr = expr.set_index(expr.columns[0])
if not all(expr.columns == meta[idCol]):
    raise ValueError(f"Column names of expression differ from '{idCol}' column in meta")
expr = expr.T
expr.insert(0, idCol, expr.index) # dysregnet wants ids in first column

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
stats = data.get_model_stats()

# feather needs str column names
result.columns = [str(x) for x in result.columns]

# Print some result statistics
print("Result:")
print(result)

print("Stats:")
print(stats)



# This works only with applied condition criterion for dysregnet
col_sums = result.apply(np.sum, axis=0)
not_zero = np.sum(np.sum(result!=0))

print(f"Number of total dysregulations: {not_zero}")
print(f"Number of edges with at least one dysregulation: {np.sum(col_sums!=0)}")
print(f"Number of positive slopes: {np.sum(col_sums>0)}")
print(f"Number of negative slopes: {np.sum(col_sums<0)}")


# Write output
result.reset_index(names="patient id", inplace=True)
result.to_feather(args.output)

stats.to_csv(args.output_stats)
