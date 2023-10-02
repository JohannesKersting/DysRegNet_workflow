import pandas as pd
import numpy as np

import argparse
import sys

import statsmodels.api as sm


from scipy.stats import mannwhitneyu

from ast import literal_eval as make_tuple


# take patients with stage 1 and (3 or 4)
def parse_stage(x):
    # if x in ['Stage IIA',"Stage IIB",'Stage I']:
    if x in ['Stage I']:
        return 'early'

    elif x in ['Stage IIIA', "Stage IIIB", 'Stage IIIC', 'Stage IV', 'Stage IVA', 'Stage IVB', 'Stage IVC', 'Stage III']:
        return 'late'

    else:
        return 'others'


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


def get_pvals(early_stage_dys, late_stage_dys, alternative, tfs):
    """
    Returns a pandas data frame with raw and bh corrected p-values for cancer stages based on a mann whitney u test.
    """

    assert alternative in ["less", "greater"]

    p_values = []
    # test for every tf
    for i, g in enumerate(tfs):
        _, p = mannwhitneyu(early_stage_dys[:, i], late_stage_dys[:, i], alternative=alternative)
        p_values.append(p)

    p_values = np.array(p_values)
    fdr_bh = sm.stats.multipletests(p_values, method='fdr_bh', alpha=0.05)[1]

    print(f"\n{alternative=}")
    print(f"{sum(p_values<0.05)=}")
    print(f"{sum(fdr_bh<0.05)=}")

    result = pd.DataFrame({"p_values": p_values, "fdr_bh": fdr_bh})
    result.index = tfs
    return result

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--network', type=str, required=True,
                        help='Paths to a patient-specific networks file')
    parser.add_argument('--meta', type=str, required=True,
                        help='Paths to file with meta information')
    parser.add_argument('--output_pvals', type=str, required=True,
                        help='Path to output file for p-values')
    args = parser.parse_args()

    print(f"Python version:\n{sys.version}")
    print(f"Arguments:\n{args}\n")

    network_path = args.network
    meta_path = args.meta
    pvals_out = args.output_pvals
    min_patients = 30

    # read and format data
    net = pd.read_feather(network_path)
    net = net.set_index("patient id")
    net.columns = [make_tuple(col) for col in net.columns]

    meta = pd.read_csv(meta_path)
    meta = meta.set_index("sample")

    # select only meta data for samples in network
    meta = meta.loc[net.index.tolist()]

    print("\nStages:")
    print(meta.ajcc_pathologic_tumor_stage.value_counts(dropna=False))

    print("\nParsed stages:")
    meta['stage'] = meta['ajcc_pathologic_tumor_stage'].apply(lambda x: parse_stage(x))
    print(meta.stage.value_counts(dropna=False))

    patients = list(meta[meta.stage.isin(['early', 'late'])].index)  # select patients with early and late stages
    unique_tfs = list(set([edge[0] for edge in net.columns.values]))  # select unique transcription factors

    dys_scores = get_dys_score_df(net, unique_tfs).T  # get dysregulation scores
    print(f"{dys_scores.shape=}")

    # separate into early and late stage
    late_stage_dys = np.array(dys_scores.loc[meta[meta.stage == 'late'].index])
    early_stage_dys = np.array(dys_scores.loc[meta[meta.stage == 'early'].index])

    print(f"{early_stage_dys.shape=}")
    print(f"{late_stage_dys.shape=}")

    # check sufficient sample number for stages
    if early_stage_dys.shape[0] < min_patients or late_stage_dys.shape[0] < min_patients:
        print(f"Insufficient sample number (min required: {min_patients})")
        sys.exit(0)

    results_less = get_pvals(early_stage_dys=early_stage_dys, late_stage_dys=late_stage_dys, alternative="less",
                             tfs=dys_scores.columns)
    results_greater = get_pvals(early_stage_dys=early_stage_dys, late_stage_dys=late_stage_dys, alternative="greater",
                                tfs=dys_scores.columns)

    # format and save results
    results_less.to_csv(pvals_out)


if __name__ == "__main__":
    main()