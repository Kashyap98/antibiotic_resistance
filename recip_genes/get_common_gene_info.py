import os

import pandas as pd

from models.gene import Gene
from utils.dir_utils import OrganismDirs, DrugDirs


def get_common_gene_info(drug, phenotype):
    drug_dirs = DrugDirs(drug, phenotype)
    # Main Organism
    organism_file = pd.read_csv(drug_dirs.target_phenotype_file, header=None)
    all_orgs = list(organism_file.iloc[:, 0])
    target_organism_dirs = OrganismDirs(all_orgs[0])

    gene_file = pd.read_csv(os.path.join(drug_dirs.drug_dir, f"{phenotype}_aa_recip_detailed_info.csv"), header=0)
    all_genes = list(gene_file.iloc[:, 0])
    with open(os.path.join(drug_dirs.drug_dir, f"{phenotype}_aa_recip_detailed_info_with_codes.csv"), "w") as unique_gene_info:
        unique_gene_info.write("gene,gene_function,hit_count,qcov_average,pident_average,perfect_matches\n")
        for i in range(len(all_genes)):
            gene = all_genes[i]
            gene_data = Gene(target_organism_dirs.organism, gene, get_info=True)
            unique_gene_info.write(f"{gene},{gene_file['gene_function'][i]},{gene_file['hit_count'][i]},"
                                   f"{gene_file['qcov_average'][i]},{gene_file['pident_average'][i]},"
                                   f"{gene_file['perfect_matches'][i]},{gene_data.evidence_codes}\n")


PHENOTYPES = ["res"]
DRUGS = ["SULF"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_common_gene_info(DRUG, PHENOTYPE)
