import os

import pandas as pd

from models.Gene import Gene
from utils.dir_utils import OrganismDirs, DrugDirs


def get_unique_gene_info(drug, phenotype):
    drug_dirs = DrugDirs(drug, phenotype)
    # Main Organism
    organism_file = pd.read_csv(drug_dirs.target_phenotype_file, header=None)
    all_orgs = list(organism_file.iloc[:, 0])
    target_organism_dirs = OrganismDirs(all_orgs[0])

    gene_file = pd.read_csv(os.path.join(drug_dirs.drug_dir, f"{phenotype}_real_full_genes.csv"), header=None)
    all_genes = list(gene_file.iloc[:, 0])
    with open(os.path.join(drug_dirs.drug_dir, f"{phenotype}_real_full_genes_detailed.csv"), "w") as unique_gene_info:
        unique_gene_info.write("Organism,Gene,Gene Length,Function,Evidence Codes\n")
        for gene in all_genes:
            gene_data = Gene(target_organism_dirs.organism, gene, get_info=True)
            unique_gene_info.write(f"{gene_data.organism},{gene_data.name},{len(gene_data.aa_sequence)},"
                                   f"{gene_data.function_data},{gene_data.evidence_codes}\n")


PHENOTYPES = ["res"]
DRUGS = ["AMOXO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_unique_gene_info(DRUG, PHENOTYPE)
