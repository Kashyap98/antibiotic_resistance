import csv
import os

import pandas as pd

from models import blast as b
from models.gene import Gene
from utils import gen_utils
from utils.dir_utils import OrganismDirs, DrugDirs, generate_dir


def compare_target_res_to_sus(drug, phenotype):
    drug_dirs = DrugDirs(drug, phenotype)
    op_phenotype = gen_utils.get_op_phenotype(phenotype)
    drug_dirs.set_opposite_phenotype_file(op_phenotype)
    target_gene_folder = os.path.join(drug_dirs.drug_dir, "target_genes")
    generate_dir(target_gene_folder)
    organism_file = pd.read_csv(drug_dirs.target_phenotype_file, header=None)
    sus_organism_file = pd.read_csv(drug_dirs.op_phenotype_file, header=None)
    all_orgs = list(organism_file.iloc[:, 0])
    sus_orgs = list(sus_organism_file.iloc[:, 0])
    all_orgs.extend(sus_orgs)
    all_orgs_copy = [""]
    all_orgs_copy.extend(list(organism_file.iloc[:, 0]))
    all_orgs_copy.extend(list(sus_organism_file.iloc[:, 0]))
    target_genes_file = pd.read_csv(os.path.join(drug_dirs.drug_dir, f"{phenotype}_real_full_genes_detailed.csv"))
    target_organism = all_orgs[0]
    all_orgs_mapped = {}
    all_genes_mapped = {}
    for gene in target_genes_file["Gene"]:
        target_gene = Gene(target_organism, gene, get_info=True)
        gene_mapped = {}
        print(gene)
        for org in all_orgs:
            recip_organism_dirs = OrganismDirs(org)
            database_path = recip_organism_dirs.database_dir
            blast_data = b.blast(target_gene, database_path, org)
            if blast_data:
                gene_mapped[org] = blast_data.blast_gene
            else:
                gene_mapped[org] = None

        all_genes_mapped[gene] = gene_mapped

    for gene in target_genes_file["Gene"]:
        all_orgs_header = [all_orgs_copy]
        all_rows = list(all_orgs_header)
        mapped_gene = all_genes_mapped[gene]
        print(gene)
        for org, org_gene in mapped_gene.items():
            current_row = [org]
            if org_gene is not None:
                for organism in all_orgs:
                    recip_organism_dirs = OrganismDirs(organism)
                    database_path = recip_organism_dirs.database_dir
                    blast_data = b.blast(org_gene, database_path, organism)
                    if blast_data:
                        current_row.append(blast_data.pident)
                    else:
                        current_row.append(0)
            else:
                current_row.append(0)
            all_rows.append(current_row)

        with open(os.path.join(target_gene_folder, f"{gene.replace('.fasta', '')}.csv"), "w", newline='') as t_file:
            wr = csv.writer(t_file, dialect='excel')
            wr.writerows(all_rows)


PHENOTYPES = ["res"]
DRUGS = ["AMOXO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        compare_target_res_to_sus(DRUG, PHENOTYPE)
