import os

import pandas as pd

# quick and dirty way used to filter mutations in earlier forms of the experiment
def get_mutations(drug, phenotype):

    # Main Organism
    mutations = pd.read_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_Mutations.csv"))
    gaps = pd.read_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_Gaps.csv"))
    data = {}

    mutations = mutations.drop_duplicates(subset="Original Gene")
    gaps = gaps.drop_duplicates(subset="Original Gene")

    # Mutations = 6, 19, Gaps = 6, 13
    for x in range(6, 20):
        col = mutations.columns[x]
        data[col] = sum(mutations[col])

    for y in range(6, 14):
        col = gaps.columns[y]
        data[col] = sum(gaps[col])

    df = pd.DataFrame(data, index=[0], columns=["nonpolar_polar", "nonpolar_acidic"
        , "nonpolar_basic", "polar_nonpolar", "polar_acidic", "polar_basic"
        , "acidic_nonpolar", "acidic_polar", "acidic_basic", "basic_nonpolar"
        , "basic_polar", "basic_acidic", "Gap Reference", "Gap Matched", "nonpolar_reference"
        , "polar_reference", "acidic_reference", "basic_refernce", "nonpolar_matched"
        , "polar_matched", "acidic_matched", "baic_matched"])

    print(df["baic_matched"])
    df.to_csv(os.path.join(os.getcwd(), drug + "_" + phenotype + "_Data.csv"), index=False)


phenotypes = ["res", "sus"]
drugs = ["AMOXO", "CIPRO", "FOSF", "SULF"]

for drug in drugs:
    for phenotype in phenotypes:
        get_mutations(drug, phenotype)
