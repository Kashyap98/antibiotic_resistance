import os

import pandas as pd
from utils.dir_utils import DrugDirs


# This is used to only get the interesting mutations
def get_mutations(drug, phenotype):
    drug_dirs = DrugDirs(drug, phenotype)
    # Main Organism
    unique_df = pd.read_csv(os.path.join(drug_dirs.drug_dir, f"{phenotype}_InfoMutations.csv"))
    unique_df = unique_df.loc[(unique_df["polar_nonpolar"] > 0) | (unique_df["polar_acidic"] > 0)
                              | (unique_df["polar_basic"] > 0) | (unique_df["nonpolar_polar"] > 0)
                              | (unique_df["nonpolar_acidic"] > 0) | (unique_df["nonpolar_basic"] > 0)
                              | (unique_df["basic_acidic"] > 0) | (unique_df["basic_polar"] > 0)
                              | (unique_df["basic_nonpolar"] > 0) | (unique_df["acidic_polar"] > 0)
                              | (unique_df["acidic_nonpolar"] > 0) | (unique_df["acidic_basic"] > 0)
                              | (unique_df["Gap Reference"] > 0) | (unique_df["Gap Matched"] > 0)]

    unique_df.to_csv(os.path.join(drug_dirs.drug_dir, f"{phenotype}_FilteredMutations.csv"), index=False)


PHENOTYPES = ["res", "sus"]
DRUGS = ["SULF"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_mutations(DRUG, PHENOTYPE)
