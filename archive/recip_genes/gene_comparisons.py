import os

import pandas as pd

from models import blast as b
from models.gene import Gene
from utils.dir_utils import OrganismDirs, DrugDirs
from utils import dir_utils
import utils.gen_utils as gen_utils


def get_uniques(drug, phenotype):
    # create output file and write the headers
    drug_dirs = DrugDirs(drug, phenotype)
    gene_finder_dir = os.path.join(drug_dirs.drug_dir, "gene_finder")
    op_phenotype = gen_utils.get_op_phenotype(phenotype)

    dir_utils.generate_dir(gene_finder_dir)

    with open(os.path.join(gene_finder_dir, f"{phenotype}_compare_topo_iv_a.csv"), "w") as total_matched_unique_genes:
        total_matched_unique_genes.write(
            "organism,resistance,gene,qcov_average,pident_average,blast_organism,blast_resistance,blast_gene\n")

    # Main Organism
    organism_file = pd.read_csv(drug_dirs.target_phenotype_file, header=None)
    all_orgs = list(organism_file.iloc[:, 0])
    target_organism_dirs = OrganismDirs(all_orgs[0], converted_ncbi_data=True)
    op_organism_file = pd.read_csv(drug_dirs.op_phenotype_file, header=None)
    op_orgs = list(op_organism_file.iloc[:, 0])

    all_orgs.extend(op_orgs)
    genes_found = []
    count = 0
    for org in all_orgs:
        org_dir = OrganismDirs(org, converted_ncbi_data=True)
        for gene in os.listdir(org_dir.gene_folder):
            temp_gene = Gene(org, gene)
            if "lactamase" in temp_gene.description:
                genes_found.append(temp_gene)
                count += 1
                print(f"{count} / {len(all_orgs)}")

    for res_gene in genes_found:

        for compare_org in all_orgs:

            compare_org_dirs = OrganismDirs(compare_org, converted_ncbi_data=True)

            if compare_org in op_orgs:
                compare_org_resistance = "sus"
            else:
                compare_org_resistance = "res"

            # First blast the first organism gene against the database of the gene in the list.
            blast_data = b.blast(res_gene, compare_org_dirs.database_dir, compare_org)

            # if there is not a hit
            if not blast_data:
                continue

            if res_gene.organism in op_orgs:
                org_resistance = "sus"
            else:
                org_resistance = "res"

            with open(os.path.join(gene_finder_dir, f"{phenotype}_compare_lactamase.csv"), "a") as total_matched_unique_genes:
                total_matched_unique_genes.write(
                    f"{res_gene.organism},{org_resistance},{res_gene.description},{blast_data.qcov},"
                    f"{blast_data.pident},{compare_org},{compare_org_resistance},{blast_data.blast_gene.description}\n")


PHENOTYPES = ["res"]
DRUGS = ["AMOXO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_uniques(DRUG, PHENOTYPE)
