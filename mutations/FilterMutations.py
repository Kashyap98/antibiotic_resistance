import os

import pandas as pd


## This is used to find the genes that are unique to each phenotype after being recip blasted to the same phenotype
def get_mutations(drug, phenotype):

    # Main Organism
    unique_df = pd.read_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_InfoMutations.csv"))
    unique_df = unique_df.loc[(unique_df["polar_nonpolar"] > 0) | (unique_df["polar_acidic"] > 0)
                              | (unique_df["polar_basic"] > 0) | (unique_df["nonpolar_polar"] > 0)
                              | (unique_df["nonpolar_acidic"] > 0) | (unique_df["nonpolar_basic"] > 0)
                              | (unique_df["basic_acidic"] > 0) | (unique_df["basic_polar"] > 0)
                              | (unique_df["basic_nonpolar"] > 0) | (unique_df["acidic_polar"] > 0)
                              | (unique_df["acidic_nonpolar"] > 0) | (unique_df["acidic_basic"] > 0)
                              | (unique_df["Gap Reference"] > 0) | (unique_df["Gap Matched"] > 0)]

    unique_df.to_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_FilteredMutations.csv"), index=False)


phenotypes = ["res", "sus"]
drugs = ["SULF"]

for drug in drugs:
    for phenotype in phenotypes:
        get_mutations(drug, phenotype)
