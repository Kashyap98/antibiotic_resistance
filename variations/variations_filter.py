import os

import pandas as pd
from utils.dir_utils import DrugDirs
import utils.gen_utils as gen_utils


def variation_filter(drug, phenotype):
    op_phenotype = gen_utils.get_op_phenotype(phenotype)
    drug_dirs = DrugDirs(drug, phenotype)
    one_float = float(1.000)

    # filter each dataframe based on the values of the lcs.
    df = pd.read_csv(os.path.join(drug_dirs.drug_dir, f"{phenotype}_Variations.csv"), index_col=0)
    boring_df = df.loc[(df['res_lcs'] == one_float) & (df['sus_lcs'] == one_float)]
    one_boring_df = df.loc[(df['res_lcs'] == one_float) & (df['sus_lcs'] != one_float)]
    opp_one_boring_df = df.loc[(df['sus_lcs'] == one_float) & (df['res_lcs'] != one_float)]
    both_varied_df = df.loc[(df['res_lcs'] != df['sus_lcs']) & (df['res_lcs'] != one_float) &
                            (df['sus_lcs'] != one_float)]
    both_varied_same_df = df.loc[(df['res_lcs'] == df['sus_lcs']) & (df['res_lcs'] != one_float)]

    # output each subset as a csv for further analysis.
    boring_df.to_csv(os.path.join(drug_dirs.variations_dir, f"{phenotype}_boring.csv"),
                     index=False)
    one_boring_df.to_csv(os.path.join(drug_dirs.variations_dir, f"{phenotype}_one_boring.csv"),
                         index=False)
    opp_one_boring_df.to_csv(os.path.join(drug_dirs.variations_dir, f"{op_phenotype}_one_boring.csv"),
                             index=False)
    both_varied_df.to_csv(os.path.join(drug_dirs.variations_dir, f"{phenotype}_both_var.csv"),
                          index=False)
    both_varied_same_df.to_csv(os.path.join(drug_dirs.variations_dir, f"{phenotype}_both_var_same.csv"),
                               index=False)


PHENOTYPES = ["res"]
DRUGS = ["CIPRO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        variation_filter(DRUG, PHENOTYPE)
