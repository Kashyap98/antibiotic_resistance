import os

import pandas as pd
from utils.dir_utils import DrugDirs


# Used to remove values that are not as important for our analysis
def get_mutations(drug, phenotype):
    drug_dirs = DrugDirs(drug, phenotype)
    # Remove all the proteins with no changes
    unique_df = pd.read_csv(os.path.join(drug_dirs.drug_dir, f"{phenotype}_InfoMutations.csv"))
    bad_proteins = unique_df["Original Gene"].loc[(unique_df["polar_nonpolar"] == 0) & (unique_df["polar_acidic"] == 0)
                                                  & (unique_df["polar_basic"] == 0) & (unique_df["nonpolar_polar"] == 0)
                                                  & (unique_df["nonpolar_acidic"] == 0) & (
                                                              unique_df["nonpolar_basic"] == 0)
                                                  & (unique_df["basic_acidic"] == 0) & (unique_df["basic_polar"] == 0)
                                                  & (unique_df["basic_nonpolar"] == 0) & (
                                                              unique_df["acidic_polar"] == 0)
                                                  & (unique_df["acidic_nonpolar"] == 0) & (
                                                              unique_df["acidic_basic"] == 0)
                                                  & (unique_df["Gap Reference"] == 0) & (unique_df["Gap Matched"] == 0)]

    print(len(set(bad_proteins)))

    for protein in bad_proteins:
        unique_df = unique_df.loc[unique_df["Original Gene"] != protein]

    # Do the same as above but with gaps
    unique_gap_df = pd.read_csv(os.path.join(drug_dirs.drug_dir, f"{phenotype}_InfoGaps.csv"))
    gap_bad_proteins = unique_gap_df["Original Gene"].loc[(unique_gap_df["polar_reference"] == 0)
                                                          & (unique_gap_df["nonpolar_reference"] == 0) & (
                                                                      unique_gap_df["basic_refernce"] == 0)
                                                          & (unique_gap_df["acidic_reference"] == 0) & (
                                                                      unique_gap_df["polar_matched"] == 0)
                                                          & (unique_gap_df["nonpolar_matched"] == 0) & (
                                                                      unique_gap_df["baic_matched"] == 0)
                                                          & (unique_gap_df["acidic_matched"] == 0)]

    for protein in gap_bad_proteins:
        unique_gap_df = unique_gap_df.loc[unique_gap_df["Original Gene"] != protein]
    print(len(set(gap_bad_proteins)))

    unique_df.to_csv(os.path.join(drug_dirs.drug_dir, f"{phenotype}__CleansedMutations.csv"), index=False)
    unique_gap_df.to_csv(os.path.join(drug_dirs.drug_dir, f"{phenotype}_CleansedGaps.csv"), index=False)


PHENOTYPES = ["res", "sus"]
DRUGS = ["AMOXO", "CIPRO", "FOSF", "SULF"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_mutations(DRUG, PHENOTYPE)
