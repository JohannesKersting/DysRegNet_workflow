import pandas as pd
import numpy as np
import math
import pyarrow

import scipy.stats as ss

import statsmodels.stats.multitest as mt
import statsmodels.formula.api as smf
from statsmodels.formula.api import ols


import sys
import os

import argparse


def flatten(t):
    return [item for sublist in t for item in sublist]


def pval(model):
    df = model.df_resid
    return(ss.t.sf(model.tvalues[1], df))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--methylation',  type=str, required=True,
                        help='Methylation data in .csv format')
    parser.add_argument('--networks', nargs='+', required=True,
                        help='Paths to multiple patient specific networks')
    parser.add_argument('--output_local', type=str, required=True,
                        help='Path to output of local tests')
    parser.add_argument('--output_global', type=str, required=True,
                        help='Path to output of global test')
    args = parser.parse_args()

    print(f"Python version:\n{sys.version}")
    print(f"Arguments:\n{args}\n")

    methylation_path = args.methylation
    input_paths = args.networks
    output_local = args.output_local
    output_global = args.output_global


    # cancer type is expected to be in the third last position of the splitted path
    input_map = {os.path.normpath(path).split(os.sep)[-3]: path for path in input_paths}
    print("\nNetwork inputs:")
    for key, value in input_map.items():
        print(f"{key}: {value}")


    print("\nReading and transposing methylation data...")
    met = pd.read_csv(methylation_path, index_col=0)
    met = met.T
    print(met)

    # output dicts
    met_test_local = dict()
    met_test_global = dict()


    for c, path in input_map.items():

        print("\n" + c)
        net = pd.read_feather(input_map[c])
        net = net.set_index("patient id")

        # recovere the reference network based on the edge names
        normal_net = []
        for edge in net.columns:
            source, target = edge[1:-1].split(", ")
            source = source[1:-1]
            target = target[1:-1]
            normal_net.append((source, target))
        normal_net = pd.DataFrame(normal_net, columns=["regulatoryGene", "targetGene"])
        normal_net = normal_net.set_index("regulatoryGene")
        normal_net_edges = [(tup[0], tup[1]) for tup in normal_net.itertuples()]

        # lists with unique transcription factors and target genes
        tfs = list(set(normal_net.index))
        tgs = list(set(normal_net.targetGene))

        # list patients with available mutation and methylation data respectively
        patients_met = list(set(met.index).intersection(set(net.index)))

        tfs_met = list(set(met.columns).intersection(set(tfs)))
        tgs_met = list(set(met.columns).intersection(set(tgs)))

        print(f"{c} Patients with methylation data: {len(patients_met)}/{len(net.index)}")
        print(f"{c} TFs with methylation data: {len(tfs_met)}/{len(tfs)}")
        print(f"{c} TGs with methylation data: {len(tgs_met)}/{len(tgs)}")

        # get methylation data and patient-specfic networks for patients with methylation data
        met_local = met.copy()
        met_local = met_local.loc[patients_met]

        net_met = net
        net_met = net_met.loc[patients_met]
        net_met = np.abs(net_met) > 0


        # local tests
        met_test_local[c] = []  # test results (tg, p-val corrected, p-val)
        dat_global = []  # data for global test
        pv_nc = []  # not corrected p-values
        genes_nc = []  # tested target genes
        fail_count = 0

        for ntf in tgs_met:

            # number/fraction of dysregulated tfs of the target gene for each patient
            in_degree = net_met[net_met.columns[normal_net.targetGene == ntf]].sum(axis=1).values
            in_degree_norm = in_degree / len(net_met.columns[normal_net.targetGene == ntf])

            # methylation values of the target gene for each patient
            methylation = met_local[ntf].values

            # append dysregulation score, methylation value, and target gene for global test
            dat_global.append([(in_degree_norm[i], methylation[i], ntf) for i in range(len(in_degree_norm))])

            # data for local test
            dat = [(in_degree_norm[i], methylation[i]) for i in range(len(in_degree_norm))]
            dat = pd.DataFrame(dat, columns=["degree", "methylation"])

            # train model, append p-value and target gene to lists
            model_full = ols("methylation ~ degree", data=dat).fit()
            if not math.isnan(model_full.pvalues[1]):
                pv_nc.append(pval(model_full))
                genes_nc.append(ntf)
            else:
                fail_count += 1

        # adjustments for individual tests
        _, pvs, _, _ = mt.multipletests(pv_nc, method="fdr_bh")
        print(
            f'{c} local methylation tests: {sum(pvs < 0.05)}/{len(pvs)} significant ({fail_count} tests produced NANs)')

        for i in range(len(pvs)):
            met_test_local[c].append((genes_nc[i], pvs[i], pv_nc[i]))


        # global test
        dat_global = flatten(dat_global)
        dat_global = pd.DataFrame(dat_global, columns=["degree", "methylation", "gene"])

        model_full = smf.mixedlm("methylation ~ degree", data=dat_global, groups=dat_global["gene"]).fit()
        pv_correct = pval(model_full)

        # randomly permutated
        dat_global["degree"] = np.random.permutation(dat_global["degree"].values)
        model_full = smf.mixedlm("methylation ~ degree", data=dat_global, groups=dat_global["gene"]).fit()
        pv_perm = pval(model_full)

        print(f'{c} global methylation test: p-value full: {pv_correct}; p-value permuted: {pv_perm}')
        met_test_global[c] = (pv_correct, pv_perm, len(tgs_met))


    # save output
    np.save(output_local, met_test_local)
    np.save(output_global, met_test_global)




if __name__ == "__main__":
    main()