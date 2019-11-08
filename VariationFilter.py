import os

import pandas as pd


def variation_filter(drug, phenotype):
    op_phenotype = "res"
    if phenotype == "res":
        op_phenotype = "sus"

    one_float = float(1.000)

    # filter each dataframe based on the values of the lcs.
    df = pd.read_csv(os.path.join(os.getcwd(), "sorted_data", drug, f"{phenotype}_Variations.csv"), index_col=0)
    boring_df = df.loc[(df['res_lcs'] == one_float) & (df['sus_lcs'] == one_float)]
    one_boring_df = df.loc[(df['res_lcs'] == one_float) & (df['sus_lcs'] != one_float)]
    opp_one_boring_df = df.loc[(df['sus_lcs'] == one_float) & (df['res_lcs'] != one_float)]
    both_varied_df = df.loc[(df['res_lcs'] != df['sus_lcs']) & (df['res_lcs'] != one_float) & (df['sus_lcs'] != one_float)]
    both_varied_same_df = df.loc[(df['res_lcs'] == df['sus_lcs']) & (df['res_lcs'] != one_float)]

    # output each subset as a csv for further analysis.
    boring_df.to_csv(os.path.join(os.getcwd(), "sorted_data", drug, "variations", f"{phenotype}_boring.csv"),
                     index=False)
    one_boring_df.to_csv(os.path.join(os.getcwd(), "sorted_data", drug, "variations", f"{phenotype}_one_boring.csv"),
                         index=False)
    opp_one_boring_df.to_csv(os.path.join(os.getcwd(), "sorted_data", drug, "variations", f"{op_phenotype}_one_boring.csv"),
                             index=False)
    both_varied_df.to_csv(os.path.join(os.getcwd(), "sorted_data", drug, "variations", f"{phenotype}_both_var.csv"),
                          index=False)
    both_varied_same_df.to_csv(os.path.join(os.getcwd(), "sorted_data", drug, "variations", f"{phenotype}_both_var_same.csv"),
                               index=False)


phenotypes = ["res"]
drugs = ["CIPRO"]

for drug in drugs:
    for phenotype in phenotypes:
        variation_filter(drug, phenotype)
