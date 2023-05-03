import os
import sys

import pandas as pd
import numpy as np
import math
import pyarrow
import statsmodels.formula.api as smf
import statsmodels.stats.multitest as mt
from ast import literal_eval as make_tuple
import scipy.stats as ss
import warnings
from pymer4.models import Lmer

import argparse


def get_mutation_dict(mutations_path):
    """
    Reads the mutations file and returns a dictionary,
    which maps a sample/patient id to a set containing all its mutated genes.
    """

    fields = ["sample", "gene", "effect"]
    drop_effects = ["Silent", "3'UTR", "5'UTR", "3'Flank", "5'Flank"]

    print("\nReading mutation data...")
    mutations_df = pd.read_csv(mutations_path, usecols=fields)
    print(mutations_df)
    print(f"Patient number: {len(set(mutations_df['sample']))}")
    print(f"\nMutation effects:\n{mutations_df.effect.value_counts()}")

    print(f"\nRemoving mutations with effect in {drop_effects}...")
    mutations_df = mutations_df.query("effect not in @drop_effects")
    print(f"Filtered mutation effects:\n{mutations_df.effect.value_counts()}")

    print("\nTransforming to dicts...")
    sample2genes = mutations_df.groupby("sample").apply(lambda x: set(x["gene"]))
    print(f"Patients with at least one mutated gene: {len(sample2genes)}")

    return sample2genes


def get_gene2samples(sample2genes, patients_mut):
    """
    Returns a dictionary, which maps the genes in sample2genes to all patient, it is mutated in.
    Only patients in 'patients_mut' are considered
    """
    gene2samples = {}
    for sample, gene_set in sample2genes.items():
        if sample in patients_mut:
            for gene in gene_set:
                if gene not in gene2samples:
                    gene2samples[gene] = set()
                gene2samples[gene].add(sample)
    return gene2samples


def get_dys_score_df(net, tfs):
    print("Calculating tf dysregulation scores...")
    dys_score_list = []

    for tf in tfs:
        # subset net
        tf_net = net[[edge for edge in net.columns if edge[0] == tf]]

        # for each patient, divide number of dysregulated edges by total number
        dys_scores = tf_net.apply(lambda x: sum(abs(x) > 0) / len(x), axis=1)
        dys_score_list.append(dys_scores)

    tf_dys_scores = pd.DataFrame(dys_score_list, index=tfs)
    print(f"TFs with std(dys_scores)==0: {sum(tf_dys_scores.std(axis=1) == 0)}/{len(tf_dys_scores.index)}")
    return tf_dys_scores


def get_mutation_df(gene2samples, tfs, patients):
    print("Creating binary mutation matrix...")
    df = pd.Series(gene2samples).str.join('|').str.get_dummies()
    df = df.loc[tfs, patients]
    return df


def pval(model):
    df = model.df_resid
    return (ss.t.sf(model.tvalues[1], df))


def local_test_pvals(mutations_df, dys_score_df):
    print("Performing local tests...")
    warnings.simplefilter('ignore')
    pvals = []
    n_failed = 0

    # iterate over tfs in both data frames
    for (tf, mut), (_, dys) in zip(mutations_df.iterrows(), dys_score_df.iterrows()):
        try:
            data = pd.DataFrame({"mut": mut.values, "dys": dys.values})
            log_reg = smf.logit("mut ~ dys", data).fit(disp=0)
            # pvals.append(log_reg.pvalues[1])
            pvals.append(pval(log_reg))
        except:
            pvals.append(1)
            n_failed += 1

    warnings.resetwarnings()
    print(f"Failed tests: {n_failed}/{len(pvals)}")
    return (pvals, n_failed)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mutations',  type=str, required=True,
                        help='Mutation data in .csv format')
    parser.add_argument('--networks', nargs='+', required=True,
                        help='Paths to multiple patient specific networks')
    parser.add_argument('--output_local', type=str, required=True,
                        help='Path to output of local tests')
    parser.add_argument('--output_global', type=str, required=True,
                        help='Path to output of global test')
    args = parser.parse_args()

    print(f"Python version:\n{sys.version}")
    print(f"Arguments:\n{args}\n")

    mutations_path = args.mutations
    input_paths = args.networks
    output_local = args.output_local
    output_global = args.output_global
    n_min_patients = 2

    # cancer type is expected to be in the third last position of the splitted path
    input_map = {os.path.normpath(path).split(os.sep)[-3]: path for path in input_paths}
    print("\nNetwork inputs:")
    for key, value in input_map.items():
        print(f"{key}: {value}")

    sample2genes = get_mutation_dict(mutations_path)

    local_test_results = []
    global_test_results = []
    for c, path in input_map.items():
        print(f"\n\n{c}")

        net = pd.read_feather(path)
        net = net.set_index("patient id")

        # parse string col names to tuples
        net.columns = [make_tuple(col) for col in net.columns]

        # select patients, with at least one mutated gene
        patients_mut = {patient for patient in net.index if patient in sample2genes.keys()}

        # get a gene2samples map, mapping mutated genes to samples of the cancer type
        gene2samples = get_gene2samples(sample2genes, patients_mut)

        # get tfs, mutated tfs, and tfs with enough mutated patients
        tfs = {tf for tf, tg in net.columns}
        tfs_mut = tfs.intersection(gene2samples.keys())
        tfs_enough_mut = [tf for tf in tfs_mut if len(gene2samples[tf]) >= n_min_patients]

        print(f"Patients with mutation data: {len(patients_mut)}/{len(net.index)}")
        print(f"TFs with > 0 mutated patients: {len(tfs_mut)}/{len(tfs)}")
        print(f"TFs with >= {n_min_patients} mutated patients: {len(tfs_enough_mut)}/{len(tfs)}")

        if len(tfs_enough_mut) == 0:
            continue

        # subset network for patients with mutation data
        net = net.query("index in @patients_mut")

        # get dysregulation score matrix (tfs as rows and patients as columns)
        dys_score_df = get_dys_score_df(net, tfs=tfs_enough_mut)

        # get mutation matrix (binary, tfs as rows and patients as columns)
        mutations_df = get_mutation_df(gene2samples, tfs=tfs_enough_mut, patients=dys_score_df.columns.values)

        # some checks
        assert dys_score_df.columns.equals(mutations_df.columns)
        assert dys_score_df.index.equals(mutations_df.index)
        print(f"{mutations_df.shape = }")
        net = None
        gene2samples = None

        # local tests

        pvals, n_failed = local_test_pvals(mutations_df=mutations_df, dys_score_df=dys_score_df)

        # multiple testing correction
        _, fdrs, _, _ = mt.multipletests(pvals, method="fdr_bh")

        # summary
        n_sign_pval = sum([pv < 0.05 for pv in pvals])
        n_sign_fdr = sum([fdr < 0.05 for fdr in fdrs])
        print(f"Significant p-values: {n_sign_pval}/{len(pvals)}")
        print(f"Significant corrected p-values: {n_sign_fdr}/{len(fdrs)}")

        local_test_results.append(pd.DataFrame({"pval": pvals, "fdr": fdrs, "cancer": c, "tf": tfs_enough_mut}))

        # global test

        if len(tfs_enough_mut) <= 1:
            continue

        # get long data
        mutations_df = pd.melt(mutations_df.T, ignore_index=False, var_name="tf", value_name="mut").reset_index(
            names="sample")
        dys_score_df = pd.melt(dys_score_df.T, ignore_index=False, var_name="tf", value_name="dys").reset_index(
            names="sample")
        assert dys_score_df["tf"].equals(mutations_df["tf"])
        assert dys_score_df["sample"].equals(mutations_df["sample"])

        data = dys_score_df
        data["mut"] = mutations_df["mut"]

        # perform test (logistic regression with random intercept for tf)
        print("Performing global test...")
        model = Lmer("mut  ~ dys  + (1|tf)", data=data, family='binomial')
        model.fit()
        print(model.summary())
        global_pval = model.coefs.loc["dys", "P-val"]

        # perform same test with randomly permuteted dysregulation scores
        data["dys"] = np.random.permutation(data["dys"].values)
        model = Lmer("mut  ~ dys  + (1|tf)", data=data, family='binomial')
        model.fit()
        permuted_global_pval = model.coefs.loc["dys", "P-val"]

        # summary
        print(f"Global p-value: {global_pval}")
        print(f"Permuted global p-value: {permuted_global_pval}")

        global_test_results.append(pd.DataFrame({"cancer": c,
                                                 "global_pval": global_pval,
                                                 "permuted_global_pval": permuted_global_pval,
                                                 "n_patients": len(patients_mut),
                                                 "n_tfs": len(tfs_enough_mut),
                                                 "n_failed": n_failed,
                                                 "n_sign_pval": n_sign_pval,
                                                 "n_sign_fdr": n_sign_fdr,
                                                 }, index=[0]))


    if local_test_results:
        local_test_results = pd.concat(local_test_results).reset_index(drop=True)
        print(f"\nLocal test results:\n{local_test_results}")
        local_test_results.to_csv(output_local, index=False)

    if global_test_results:
        global_test_results = pd.concat(global_test_results).reset_index(drop=True)
        print(f"\nGlobal test results:\n{global_test_results.to_string()}")
        global_test_results.to_csv(output_global, index=False)


if __name__ == "__main__":
    main()

