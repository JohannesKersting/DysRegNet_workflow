import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from munkres import Munkres
import argparse
import os
import sys
import math
import statsmodels.stats.multitest as mt
import colorcet as cc
from commons import method_lut, network_lut


def get_input_info(path):
    '''
    Returns a tuple for each input path to extract information from it.
    Order is path, cancer, method, norm_method, ref_network_start
    '''

    # split path by separator
    path = os.path.normpath(path)
    splitted_path = path.split(os.sep)

    # get method and filename
    file_name = splitted_path[-1]
    method = splitted_path[-2]
    cancer = splitted_path[-3]

    # get overlap_type, norm method and the reference network
    stem_name, file_extension = os.path.splitext(file_name)
    norm_method, ref_network, _ = stem_name.split("-")

    # shorten refernece network name (for genie3 shared and individual)
    ref_network_split = ref_network.split("_")
    if (len(ref_network_split) > 1):
        ref_network_start = "_".join(ref_network_split[0:2])
    else:
        ref_network_start = ref_network_split[0]

    return path, cancer, method, norm_method, ref_network_start


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs',  nargs='+', required=True,
                        help='.csv stats files produced by dysregnet')
    parser.add_argument('--norm_method',  type=str, required=True,
                        help='Normalization method, used to name the output files')
    parser.add_argument('--output_dir',  type=str, required=True,
                        help='Output directory for the plots')
    parser.add_argument('--output_stats_summary',  type=str, required=True,
                        help='Output file for saving a stat summary in .tsv')

    args = parser.parse_args()

    print(f"Python version:\n{sys.version}")
    print(f"Arguments:\n{args}\n")

    input_paths = args.inputs
    norm_method = args.norm_method
    output_dir = args.output_dir
    output_stats_summary = args.output_stats_summary


    # get input path infos
    input_infos = [get_input_info(path) for path in input_paths]

    input_df = pd.DataFrame(input_infos, columns=["path", "cancer", "method", "norm_method", "network"])
    input_df["method"] = input_df["method"].map(lambda x: method_lut[x] if x in method_lut.keys() else x)
    input_df["network"] = input_df["network"].map(
        lambda x: "genie3_individual" if (("genie3" in x) and (not "shared" in x)) else x)
    input_df["network"] = input_df["network"].map(lambda x: network_lut[x] if x in network_lut.keys() else x)

    print(input_df.to_string())


    # Read input files
    stats_data = []
    for index, row in input_df.iterrows():
        df = pd.read_csv(row["path"], index_col=0)
        df["edge"] = df.index
        df["Cancer"] = row["cancer"]
        df["network"] = row["network"]
        stats_data.append(df)

    stats_df = pd.concat(stats_data, ignore_index=True)
    stats_data = None
    print(stats_df)


    # R2 distribution plot
    sns.set(style="whitegrid")
    cancer_types = sorted(list(set(input_df["cancer"])))
    palette = sns.color_palette(cc.glasbey_light, n_colors=len(cancer_types) + 1)[1:]
    palette_dict = dict(zip(cancer_types, palette))

    g = sns.FacetGrid(data=stats_df, col="network", height=3, aspect=1.6, col_wrap=2)
    g.set_xticklabels(rotation=45)
    g.map_dataframe(sns.boxplot, x="Cancer", y="R2", hue="Cancer", hue_order=cancer_types, order=cancer_types,
                    palette=palette_dict, dodge=False)
    g.set_ylabels(r"$R^2$")
    g.set_titles(col_template="{col_name}")
    g.add_legend(title="Cancer")
    g.tight_layout()

    plt.savefig(os.path.join(output_dir, f"model_stats-{norm_method}-R2.png"), dpi=300)
    plt.savefig(os.path.join(output_dir, f"model_stats-{norm_method}-R2.pdf"), dpi=300)


    # TF confounder distribution plot
    sns.set(style="whitegrid")
    cancer_types = sorted(list(set(input_df["cancer"])))
    palette = sns.color_palette(cc.glasbey_light, n_colors=len(cancer_types) + 1)[1:]
    palette_dict = dict(zip(cancer_types, palette))

    tf_stats_df = stats_df[["coef_TF", "pval_TF", "network", "Cancer"]].copy()
    tf_stats_df = pd.melt(tf_stats_df, id_vars=["Cancer", "network"])

    g = sns.FacetGrid(tf_stats_df, col="network", row="variable", height=3, aspect=1.2, sharey="row",
                      margin_titles=True)
    g.set_xticklabels(rotation=45)
    g.map_dataframe(sns.violinplot, x="Cancer", y="value", hue="Cancer", hue_order=cancer_types, order=cancer_types,
                    palette=palette_dict, dodge=False, cut=0)

    g.set_titles(col_template="{col_name}", row_template="")
    g.set_xlabels("")
    g.add_legend(title="Cancer")

    for (row_val, col_val), ax in g.axes_dict.items():
        if row_val == "coef_TF":
            ax.set_ylabel("TF coefficient")
            ax.axhline(0)
        elif row_val == "pval_TF":
            ax.set_ylabel("TF p-value")
            ax.axhline(0.05)

    g.tight_layout()

    plt.savefig(os.path.join(output_dir, f"model_stats-{norm_method}-TF.png"), dpi=300)
    plt.savefig(os.path.join(output_dir, f"model_stats-{norm_method}-TF.pdf"), dpi=300)
    tf_stats_df = None


    # all confounders distribution plots
    networks = list(set(input_df["network"]))

    count_sign_list = []
    pval_stats_list = []
    for network in networks:
        coef_stats_df = stats_df[stats_df["network"] == network].copy()

        # get column name mappings
        coef_stats_slim_df = coef_stats_df.drop(columns=["network", "R2", "Cancer", "edge"])
        type_lut = {col: col.split("_")[0] for col in coef_stats_slim_df.columns.values}
        confounder_lut = {col: col.split("_")[1] for col in coef_stats_slim_df.columns.values}
        coef_stats_slim = None

        # cast to long data and add categories
        coef_stats_df = coef_stats_df.drop(columns=["R2", "edge"])
        coef_stats_df = pd.melt(coef_stats_df, id_vars=["Cancer", "network"])
        coef_stats_df["type"] = coef_stats_df["variable"].map(type_lut)
        coef_stats_df["confounder"] = coef_stats_df["variable"].map(confounder_lut)
        coef_stats_df.drop(columns="variable", inplace=True)

        # rename birth category
        coef_stats_df.loc[coef_stats_df["confounder"] == "birth", "confounder"] = "age"

        # rename ethnicity category
        coef_stats_df.loc[coef_stats_df["confounder"] == "race", "confounder"] = "ethnicity"

        # rename gender category
        coef_stats_df.loc[coef_stats_df["confounder"] == "gender", "confounder"] = "sex"

        # drop NAs
        coef_stats_df = coef_stats_df.dropna()

        # add p-values to pval_stats_list
        pval_stats_list.append(coef_stats_df.query("type=='pval' and confounder!='intercept'").drop("type", axis=1))

        sns.set(style="whitegrid")
        cancer_types = sorted(list(set(input_df["cancer"])))
        palette = sns.color_palette(cc.glasbey_light, n_colors=len(cancer_types) + 1)[1:]
        palette_dict = dict(zip(cancer_types, palette))

        g = sns.FacetGrid(coef_stats_df[coef_stats_df["confounder"]!="intercept"], row="type", height=3, aspect=5, sharey="row")
        g.map_dataframe(sns.violinplot, x="confounder", y="value", hue="Cancer", hue_order=cancer_types,
                        palette=palette_dict, cut=0)

        g.set_titles(col_template="", row_template="")
        g.set_xlabels("")
        g.add_legend(title="Cancer")

        for row_val, ax in g.axes_dict.items():
            if row_val == "coef":
                ax.set_ylabel("Coefficient")
                ax.axhline(0)
            elif row_val == "pval":
                ax.set_ylabel("P-value")
                ax.axhline(0.05)

        g.tight_layout()

        plt.savefig(os.path.join(output_dir, f"model_stats-{norm_method}-confounders-{network}.png"), dpi=300)
        plt.savefig(os.path.join(output_dir, f"model_stats-{norm_method}-confounders-{network}.pdf"), dpi=300)


        # plot without coeficcients
        g = sns.FacetGrid(coef_stats_df.query("confounder!='intercept' and type!='coef'"), row="type", height=3,
                          aspect=5, sharey="row")
        g.map_dataframe(sns.violinplot, x="confounder", y="value", hue="Cancer", hue_order=cancer_types,
                        palette=palette_dict, cut=0)

        g.set_titles(col_template="", row_template="")
        g.set_xlabels("")
        g.add_legend(title="Cancer")

        for row_val, ax in g.axes_dict.items():
            if row_val == "coef":
                ax.set_ylabel("Coefficient")
                ax.axhline(0)
            elif row_val == "pval":
                ax.set_ylabel("P-value")
                ax.axhline(0.05)

        g.tight_layout()

        plt.savefig(os.path.join(output_dir, f"model_stats-{norm_method}-pvals-{network}.png"), dpi=300)
        plt.savefig(os.path.join(output_dir, f"model_stats-{norm_method}-pvals-{network}.pdf"), dpi=300)
        plt.show()

        # count significant p-values
        count_sign_df = pd.DataFrame(
            columns=["network", "cancer", "confounder", "count", "sign_pvals", "sign_pvals_corrected"]
        )

        for (cancer, confounder), group in coef_stats_df[coef_stats_df["type"] == "pval"].groupby(
                ["Cancer", "confounder"], sort=False):
            pvals = group["value"]
            _, pvals_corrected, _, _ = mt.multipletests(pvals, method="fdr_bh")

            count = len(pvals)
            sign_pvals = sum(pvals < 0.05)
            sign_pvals_corrected = sum(pvals_corrected < 0.05)

            row = pd.DataFrame({
                "network": network,
                "cancer": cancer,
                "confounder": confounder,
                "count": count,
                "sign_pvals": sign_pvals,
                "sign_pvals_corrected": sign_pvals_corrected}, index=[0])

            count_sign_df = pd.concat([count_sign_df, row], ignore_index=True)

        count_sign_df["frac_sign_pvals"] = count_sign_df["sign_pvals"] / count_sign_df["count"]
        count_sign_df["frac_sign_pvals_corrected"] = count_sign_df["sign_pvals_corrected"] / count_sign_df["count"]

        count_sign_list.append(count_sign_df)

    count_sign_df = pd.concat(count_sign_list, ignore_index=True)
    pval_stats_df = pd.concat(pval_stats_list, ignore_index=True)

    # save model stats summary
    count_sign_df.to_csv(output_stats_summary, index=False, sep="\t")

   # plot p-value distributions across all networks
    sns.set(style="whitegrid")
    g = sns.FacetGrid(pval_stats_df, row="network", height=3, aspect=5, sharey="row")
    g.map_dataframe(sns.violinplot, x="confounder", y="value", hue="Cancer", hue_order=cancer_types, palette=palette_dict, cut=0)
    g.set_titles(col_template="{col_name}", row_template="{row_name}")
    g.set_xlabels("")
    g.set_ylabels("P-value")
    g.add_legend(title="Cancer")

    for row_val, ax in g.axes_dict.items():
        ax.axhline(0.05)

    g.tight_layout()

    plt.savefig(os.path.join(output_dir, f"model_stats-{norm_method}-pvals-all.png"), dpi=300)
    plt.savefig(os.path.join(output_dir, f"model_stats-{norm_method}-pvals-all.pdf"), dpi=300)



if __name__ == "__main__":
    main()