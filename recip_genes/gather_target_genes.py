import csv
import os

import pandas as pd

from models import blast as b
from models.gene import Gene
from utils import gen_utils
from utils.dir_utils import OrganismDirs, DrugDirs


def gather_target_genes(drug, phenotype):
    drug_dirs = DrugDirs(drug, phenotype)
    op_phenotype = gen_utils.get_op_phenotype(phenotype)
    drug_dirs.set_opposite_phenotype_file(op_phenotype)

    organism_file = pd.read_csv(drug_dirs.target_phenotype_file, header=None)
    sus_organism_file = pd.read_csv(drug_dirs.op_phenotype_file, header=None)
    all_orgs = list(organism_file.iloc[:, 0])
    sus_orgs = list(sus_organism_file.iloc[:, 0])
    target_org = all_orgs.pop(0)
    target_organism_dirs = OrganismDirs(target_org)
    target_genes_file = pd.read_csv(os.path.join(drug_dirs.drug_dir, f"{phenotype}_real_full_genes_detailed.csv"))
    all_rows = [["organism", "genes", "function_data", "nuc sequence"]]
    for gene in target_genes_file["Gene"]:
        current_gene = Gene(target_org, gene, get_info=True)
        all_rows.append([target_org, gene, current_gene.function_data, current_gene.aa_sequence])
        print(current_gene.name)
        for organism in all_orgs:
            recip_organism_dirs = OrganismDirs(organism)
            database_path = recip_organism_dirs.database_dir
            blast_data = b.blast(current_gene, database_path, organism)
            if blast_data:
                all_rows.append([organism, blast_data.blast_gene.name, blast_data.blast_gene.function_data,
                                 blast_data.blast_gene.aa_sequence])
        all_rows.append([])

    with open(os.path.join(drug_dirs.drug_dir, f"{phenotype}_aa_target_gene_info.csv"), "w", newline='') as t_file:
        wr = csv.writer(t_file, dialect='excel')
        wr.writerows(all_rows)


PHENOTYPES = ["res"]
DRUGS = ["AMOXO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        gather_target_genes(DRUG, PHENOTYPE)
