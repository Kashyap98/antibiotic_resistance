import os

import pandas as pd

from models.Gene import Gene
from utils.dir_utils import DrugDirs


def get_mutations(drug, phenotype):
    drug_dirs = DrugDirs(drug, phenotype)
    # Main Organism
    data_file = pd.read_csv(os.path.join(drug_dirs.drug_dir, f"{phenotype}_Mutations.csv"))
    print(data_file["Original Organism"][0])
    data_file.insert(loc=len(data_file.columns), column="Original Gene Function", value=None)

    for index, row in data_file.iterrows():
        original_gene = Gene(row["Original Organism"], row["Original Gene"].replace(".fasta", ""), get_info=True)
        data_file.set_value(index, "Original Gene Function", original_gene.function_data)

    data_file.to_csv(os.path.join(drug_dirs.drug_dir, f"{phenotype}_InfoMutations.csv"), index=False)


PHENOTYPES = ["res"]
DRUGS = ["CIPRO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_mutations(DRUG, PHENOTYPE)
