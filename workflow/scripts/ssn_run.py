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
from scipy import stats

import warnings
warnings.filterwarnings("ignore")


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
                    help='The output file path and name.'
                    )
parser.add_argument("--pvalue", help="Define a p-value threshold for edges", default=0.005,
                    required=False)

args = parser.parse_args()


print("Python version:")
print(sys.version)
print("Arguments:")
print(args)
print("\n")


conCol = "condition"



def run_ssn(expression_data, meta, GRN, conCol, pvalue=0.005):


    # Filter are genes that aren't included in both, network and expression data
    GRN_genes = list(set(GRN.iloc[:, 0].values.tolist() + GRN.iloc[:, 1].values.tolist()))
    GRN_genes = [g for g in GRN_genes if g in expression_data.index]

    if not GRN_genes:
        raise ValueError('Gene id or name in GRN DataFrame do not match the ones in expression_data DataFrame')

    expression_data = expression_data.loc[GRN_genes]
    GRN = GRN[GRN.iloc[:, 0].isin(GRN_genes)]
    GRN = GRN[GRN.iloc[:, 1].isin(GRN_genes)].drop_duplicates()

    case = expression_data[meta[meta[conCol]==1].index]
    control = expression_data[meta[meta[conCol]==0].index]
    normal_net = GRN

    genes = list(case.index)
    genes_dict = {genes[i]: i for i in range(len(genes))}
    patients = list(case.columns)
    n = len(control.columns)
    genes = set(genes)

    print(f"Number of unique genes in expr AND grn: {len(genes)}")

    cor = np.corrcoef(control)
    crc = np.zeros((len(patients), normal_net.shape[0]))
    pat_count = 0
    for pat in tqdm(patients):
        data_p = control.copy()
        data_p[pat] = case[pat]
        corr_p = np.corrcoef(data_p)
        dcor = corr_p - cor
        zscore = dcor / ((1 - cor ** 2) / (n - 1))
        pv = stats.norm.sf(zscore)
        pv = pv < pvalue / len(patients)
        edge_count = 0
        edges = []
        for tup in normal_net.itertuples():
            if tup[1] in genes and tup[2] in genes:

                edge = (genes_dict[tup[1]], genes_dict[tup[2]])
                if pv[edge]:
                    crc[pat_count, edge_count] = 1

            edge_count += 1
            edges.append((tup[1], tup[2]))

        pat_count += 1

    edges = [str(x) for x in edges]
    crc = pd.DataFrame(crc, columns=edges)
    crc['patient id'] = patients

    return crc




# Read data
meta = pd.read_csv(args.meta, index_col = 0)
expr = pd.read_csv(args.expr, index_col = 0)
grn = pd.read_csv(args.grn)

# Pre-process data

if not all(expr.columns == meta.index):
    raise ValueError(f"Column names of expression differ from row names in meta")

grn = grn.iloc[:, 0:2]
meta = meta[[conCol]]


# Print data summary
print("\nData as submitted to ssn_run:")
print("\nexpr:")
print(expr)

print("\nmeta:")
print(meta)

print("\ngrn:")
print(grn)

print("\nMeta composition:\n")
print(meta.groupby(conCol).size())
print("\n")


result = run_ssn(expression_data=expr, meta=meta, GRN=grn, conCol=conCol, pvalue=args.pvalue)

# feather needs str column names
result.columns = [str(x) for x in result.columns]

# Print some result statistics
print("Result:")
print(result)


col_sums = result[result.columns[1:]].apply(np.sum, axis=0)
not_zero = np.sum(np.sum(result[result.columns[1:]]!=0))

print(f"Number of total dysregulations: {not_zero}")
print(f"Number of edges with at least one dysregulation: {np.sum(col_sums!=0)}")


# Write output
result.to_feather(args.output)
