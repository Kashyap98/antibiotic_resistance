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

    unique_genes_folder = os.path.join(drug_dirs.drug_dir, "unique_genes")

    # Main Organism
    organism_file = pd.read_csv(drug_dirs.target_phenotype_file, header=None)
    all_orgs = list(organism_file.iloc[:, 0])
    count = 0

    for org in all_orgs:
        count += 1
        print(count)
        res_matches = pd.read_csv(os.path.join(unique_genes_folder, f"{org}_unique.csv"), header=0)
        genes = res_matches["gene"]
        gene_descriptions = res_matches["gene_description"]
        gene_list = {}

        for i in range(len(genes)):
            description = gene_descriptions[i]
            gene = genes[i]
            formatted_description = description.replace(f"{gene.replace('.fasta', '')} ", "").replace(f" MULTISPECIES: ", "")
            print(formatted_description)
            if "hypothetical protein" in formatted_description:
                continue

            target_gene = Gene(org, gene)
            for comparison_org in all_orgs:
                if comparison_org == org:
                    continue
                compare_org_dirs = OrganismDirs(comparison_org, converted_ncbi_data=True)

                # First blast the first organism gene against the database of the gene in the list.
                blast_data = b.blast(target_gene, compare_org_dirs.database_dir, comparison_org)

                # if there is a hit
                if not blast_data:
                    if formatted_description not in gene_list:
                        gene_list[formatted_description] = [gene, formatted_description, 0, 0, 0, []]
                    continue

                if formatted_description in gene_list:
                    old_data = gene_list[formatted_description]
                    gene_name = old_data[0]
                    gene_description = old_data[1]
                    new_qcov = (old_data[2] + blast_data.qcov) / 2
                    new_pident = (old_data[3] + blast_data.pident) / 2
                    new_organism_count = old_data[4] + 1
                    new_organisms = old_data[5]
                    new_organisms.append(comparison_org)

                    gene_list[formatted_description] = [gene_name, gene_description, new_qcov, new_pident,
                                                        new_organism_count, new_organisms]
                else:
                    gene_list[formatted_description] = [gene, formatted_description, blast_data.qcov,
                                                        blast_data.pident, 1, [comparison_org]]

        with open(os.path.join(unique_genes_folder, f"{org}_unique_matches_to_res.csv"), "w") as total_matched_unique_genes:
            total_matched_unique_genes.write("gene,gene_description,qcov_average,pident_average,organism_count,organisms\n")
            for gene_desc, gene_data in gene_list.items():
                    total_matched_unique_genes.write(f"{gene_data[0]},{gene_data[1]},{str(gene_data[2])},"
                                                     f"{str(gene_data[3])},{str(gene_data[4])},"
                                                     f"{str(gene_data[5]).replace(',', '.')}\n")


PHENOTYPES = ["res"]
DRUGS = ["CIPRO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_uniques(DRUG, PHENOTYPE)
