import os

import pandas as pd

from models import Blast as b
from models.Gene import Gene
from utils import dir_utils
from utils.dir_utils import OrganismDirs, DrugDirs
import utils.gen_utils as gen_utils


def get_uniques(drug, phenotype):
    # create output file and write the headers
    drug_dirs = DrugDirs(drug, phenotype)
    op_phenotype = gen_utils.get_op_phenotype(phenotype)
    drug_dirs.set_opposite_phenotype_file(op_phenotype)

    unique_genes_folder = os.path.join(drug_dirs.drug_dir, "unique_genes")
    dir_utils.generate_dir(unique_genes_folder)

    # Main Organism
    organism_file = pd.read_csv(drug_dirs.target_phenotype_file, header=None)
    all_orgs = list(organism_file.iloc[:, 0])
    target_organism_dirs = OrganismDirs(all_orgs[0])
    op_organism_file = pd.read_csv(drug_dirs.op_phenotype_file, header=None)
    op_orgs = list(op_organism_file.iloc[:, 0])

    count = 0
    # for org in all_orgs[4:]:
    org = all_orgs[0]
    res_organism_paths = OrganismDirs(org)
    res_genes = os.listdir(res_organism_paths.gene_folder)
    res_genes_copy = {}

    with open(os.path.join(unique_genes_folder, f"{org}_unique.csv"), "w") as organism_unique_genes:
        organism_unique_genes.write(
            "gene,hit_count,qcov_average,pident_average,perfect_matches,gene_length,function\n")

    for res_gene_name in res_genes:
        res_gene = Gene(org, res_gene_name, get_info=True)
        print(res_gene_name)
        for op_organism in op_orgs:
            op_organism = str(op_organism)
            op_organism_dirs = OrganismDirs(op_organism)

            # First blast the first organism gene against the database of the gene in the list.
            blast_data = b.blast(res_gene, op_organism_dirs.database_dir, op_organism)

            # if there is a hit
            if not blast_data:
                if res_gene_name not in res_genes_copy:
                    res_genes_copy[res_gene_name] = [0, 0, 0, 0, res_gene.length, res_gene.function_data]
                continue

            if res_gene_name in res_genes_copy:
                old_data = res_genes_copy[res_gene_name]
                new_count = old_data[0] + 1
                new_qcov = (old_data[1] + blast_data.qcov) / 2
                new_pident = (old_data[2] + blast_data.pident) / 2
                if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                    new_perfect_matches = old_data[3] + 1
                else:
                    new_perfect_matches = old_data[3]
                res_genes_copy[res_gene_name] = [new_count, new_qcov, new_pident, new_perfect_matches,
                                                 res_gene.length, res_gene.function_data]
            else:
                hit_count = 0
                if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                    perfect_matches = 1
                else:
                    perfect_matches = 0
                res_genes_copy[res_gene_name] = [hit_count, blast_data.qcov, blast_data.pident, perfect_matches,
                                                 res_gene.length, res_gene.function_data]

    with open(os.path.join(unique_genes_folder, f"{org}_unique.csv"), "a") as organism_unique_genes:
        for gene, info in res_genes_copy.items():
            if info[3] == 0:
                organism_unique_genes.write(f"{gene},{str(info[0])},{str(info[1])},{str(info[2])},{str(info[3])},"
                                            f"{str(info[4])},{str(info[5])}\n")

    count += 1
    print(f"{count} / {len(all_orgs)}")


PHENOTYPES = ["res"]
DRUGS = ["AMOXO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_uniques(DRUG, PHENOTYPE)
