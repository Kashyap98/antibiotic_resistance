import csv
import os
import threading
import time
from queue import Queue

import pandas as pd

from models import blast as b
from models.gene import Gene
from utils.dir_utils import OrganismDirs, DrugDirs
from utils import dir_utils, gen_utils

print_lock = threading.Lock()


def get_recips(drug, phenotype):
    # Create output file
    drug_dirs = DrugDirs(drug, phenotype)
    complete_comparisons_dir = os.path.join(drug_dirs.drug_dir, "complete_comparison")

    # Get all the organisms for the phenotype
    all_organisms_raw = list(csv.reader(open(drug_dirs.target_phenotype_file), delimiter=','))
    res_organisms = []

    for org in all_organisms_raw:
        res_organisms.append(org[0])

    # get all the genes for the target organism
    print(f"{phenotype} Organisms: ", res_organisms)

    complete_df = pd.DataFrame()
    filtered_df = pd.DataFrame()

    for organism_file in os.listdir(complete_comparisons_dir):

        organism_file_path = os.path.join(complete_comparisons_dir, organism_file)
        organism_data = pd.read_csv(organism_file_path, header=0)

        complete_df = pd.concat([complete_df, organism_data])
        temp_filtered_df = pd.concat([organism_data[organism_data['gene_description'].str.contains("topoisomerase")],
                                      organism_data[organism_data['gene_description'].str.contains("gyrase subunit")]])
        filtered_df = pd.concat([filtered_df, temp_filtered_df])

    complete_df.to_csv(os.path.join(drug_dirs.drug_dir, "complete_organism_comparison.csv"), sep=',', index=False)
    filtered_df.to_csv(os.path.join(drug_dirs.drug_dir, "filtered_organism_comparison.csv"), sep=',', index=False)


PHENOTYPES = ["res"]
DRUGS = ["CIPRO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_recips(DRUG, PHENOTYPE)
