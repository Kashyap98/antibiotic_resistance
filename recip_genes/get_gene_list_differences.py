import os

import pandas as pd

from models.Gene import Gene
from utils.dir_utils import OrganismDirs, DrugDirs


def get_gene_list_difference(drug, phenotype):
    drug_dirs = DrugDirs(drug, phenotype)

    # Main Organism
    organism_file = pd.read_csv(drug_dirs.target_phenotype_file, header=None)
    all_orgs = list(organism_file.iloc[:, 0])

    gene_file = pd.read_csv(os.path.join(drug_dirs.drug_dir, f"{phenotype}_aa_recip_detailed_info.csv"), header=0)
    sus_gene_file = pd.read_csv(os.path.join(drug_dirs.drug_dir, f"sus_aa_recip_detailed_info.csv"), header=0)
    all_res_genes = list(gene_file.iloc[:, 1])
    all_sus_genes = list(sus_gene_file.iloc[:, 1])

    all_res_genes_names = list(gene_file.iloc[:, 0])
    all_sus_genes_names = list(sus_gene_file.iloc[:, 0])

    for x in range(len(all_res_genes)):
        all_res_genes[x] = all_res_genes[x].replace(f' {all_res_genes_names[x].replace(".fasta", "")} MULTISPECIES: ', "").replace(f' {all_res_genes_names[x].replace(".fasta", "")} ', "")
    for x in range(len(all_sus_genes)):
        all_sus_genes[x] = all_sus_genes[x].replace(f' {all_sus_genes_names[x].replace(".fasta", "")} MULTISPECIES: ', "").replace(f' {all_sus_genes_names[x].replace(".fasta", "")} ', "")
    with open(os.path.join(drug_dirs.drug_dir, f"{phenotype}_aa_difference_genes.csv"), "w") as unique_gene_info:
        unique_gene_info.write("gene,gene_description,hit_count,qcov_average,pident_average,perfect_matches\n")
        for i in range(len(all_res_genes)):
            gene = all_res_genes[i]
            if gene not in all_sus_genes:
                unique_gene_info.write(f"{gene_file['gene_description'][i]},{gene},{gene_file['hit_count'][i]},"
                                       f"{gene_file['qcov_average'][i]},{gene_file['pident_average'][i]},"
                                       f"{gene_file['perfect_matches'][i]}\n")


PHENOTYPES = ["res"]
DRUGS = ["AMOXO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_gene_list_difference(DRUG, PHENOTYPE)
