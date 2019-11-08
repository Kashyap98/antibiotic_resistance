import os

import pandas as pd


# Used to remove values that are not as important for our analysis
def get_mutations(drug, phenotype):

    bad_proteins = []
    # Main Organism
    unique_df = pd.read_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_InfoMutations.csv"))
    bad_proteins = unique_df["Original Gene"].loc[(unique_df["polar_nonpolar"] == 0) & (unique_df["polar_acidic"] == 0)
                              & (unique_df["polar_basic"] == 0) & (unique_df["nonpolar_polar"] == 0)
                              & (unique_df["nonpolar_acidic"] == 0) & (unique_df["nonpolar_basic"] == 0)
                              & (unique_df["basic_acidic"] == 0) & (unique_df["basic_polar"] == 0)
                              & (unique_df["basic_nonpolar"] == 0) & (unique_df["acidic_polar"] == 0)
                              & (unique_df["acidic_nonpolar"] == 0) & (unique_df["acidic_basic"] == 0)
                              & (unique_df["Gap Reference"] == 0) & (unique_df["Gap Matched"] == 0)]

    print(len(set(bad_proteins)))

    for protein in bad_proteins:
        unique_df = unique_df.loc[unique_df["Original Gene"] != protein]

    gap_bad_proteins = []
    unique_gap_df = pd.read_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_InfoGaps.csv"))
    gap_bad_proteins = unique_gap_df["Original Gene"].loc[(unique_gap_df["polar_reference"] == 0)
                              & (unique_gap_df["nonpolar_reference"] == 0) & (unique_gap_df["basic_refernce"] == 0)
                              & (unique_gap_df["acidic_reference"] == 0) & (unique_gap_df["polar_matched"] == 0)
                              & (unique_gap_df["nonpolar_matched"] == 0) & (unique_gap_df["baic_matched"] == 0)
                              & (unique_gap_df["acidic_matched"] == 0)]

    for protein in gap_bad_proteins:
        unique_gap_df = unique_gap_df.loc[unique_gap_df["Original Gene"] != protein]
    print(len(set(gap_bad_proteins)))

    unique_df.to_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_CleansedMutations.csv"), index=False)
    unique_gap_df.to_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_CleansedGaps.csv"), index=False)


phenotypes = ["res", "sus"]
drugs = ["AMOXO", "CIPRO", "FOSF", "SULF"]

for drug in drugs:
    for phenotype in phenotypes:
        get_mutations(drug, phenotype)
