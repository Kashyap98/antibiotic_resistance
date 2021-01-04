import os

import pandas as pd

from models import blast as b
from models.gene import Gene
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
    genes_to_skip = []
    final_output = []

    gene_list_by_organism = {}
    for org in all_orgs:
        res_matches = pd.read_csv(os.path.join(unique_genes_folder, f"{org}_unique.csv"), header=0)
        res_genes = []
        for res_match in res_matches["gene"]:
            res_genes.append(res_match)
        gene_list_by_organism[org] = res_genes

    count = 0
    for org in all_orgs:
        res_matches = pd.read_csv(os.path.join(unique_genes_folder, f"{org}_unique.csv"), header=0)
        res_genes = res_matches["gene"]
        org_output = []

        for res_gene_name in res_genes:
            if res_gene_name in genes_to_skip:
                continue
            else:
                genes_to_skip.append(res_gene_name)

            res_gene = Gene(org, res_gene_name, get_info=False)

            if "hypothetical protein" in res_gene.function_data:
                continue
            # print(res_gene_name)
            gene_output = []

            for compare_org in all_orgs:
                if compare_org == org:
                    continue
                compare_org_dirs = OrganismDirs(compare_org, converted_ncbi_data=True)

                # First blast the first organism gene against the database of the gene in the list.
                blast_data = b.blast(res_gene, compare_org_dirs.database_dir, compare_org)

                # if there is not a hit
                if not blast_data:
                    continue

                blast_gene_with_fasta = blast_data.gene_name + "fasta"

                if blast_gene_with_fasta in genes_to_skip:
                    continue

                if blast_gene_with_fasta in gene_list_by_organism[compare_org]:
                    genes_to_skip.append(blast_gene_with_fasta)

                    if len(gene_output) != 0:
                        old_data = gene_output.copy()
                        new_count = old_data[0] + 1
                        new_qcov = (old_data[1] + blast_data.qcov) / 2
                        new_pident = (old_data[2] + blast_data.pident) / 2

                        if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                            new_perfect_matches = old_data[3] + 1
                            new_perfect_matches_organisms = old_data[7]
                            new_perfect_matches_organisms.append(compare_org)
                            new_organisms = old_data[6]
                        else:
                            new_perfect_matches = old_data[3]
                            new_perfect_matches_organisms = old_data[7]
                            new_organisms = old_data[6]
                            new_organisms.append(compare_org)

                        gene_output = [new_count, new_qcov, new_pident, new_perfect_matches, res_gene.length,
                                       res_gene.function_data, new_organisms, new_perfect_matches_organisms, res_gene_name]

                    else:
                        hit_count = 1
                        if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                            perfect_matches = 1
                            perfect_match_organisms = [compare_org]
                            organisms = []
                        else:
                            perfect_matches = 0
                            perfect_match_organisms = []
                            organisms = [compare_org]
                        gene_output = [hit_count, blast_data.qcov, blast_data.pident, perfect_matches, res_gene.length,
                                       res_gene.function_data, organisms, perfect_match_organisms, res_gene_name]

            org_output.append(gene_output)

        with open(os.path.join(unique_genes_folder, f"{org}_unique_matches_to_res.csv"), "w") as organism_unique_genes:
            organism_unique_genes.write("gene,hit_count,qcov_average,pident_average,perfect_matches,gene_length,"
                                        "function,matched_organisms,perfect_matched_organisms\n")
            for match in org_output:
                if len(match) > 0:
                    organism_unique_genes.write(f"{match[8]},{str(match[0])},{str(match[1])},{str(match[2])},"
                                                f"{str(match[3])},{str(match[4])},{str(match[5])},"
                                                f"{str(match[6]).replace(',', '.')},{str(match[7]).replace(',', '.')}\n")

        final_output.extend(org_output)
        count += 1
        print(f"{count} / {len(all_orgs)}")

    with open(os.path.join(unique_genes_folder, "total_unique_matches_to_res.csv"), "w") as total_matched_unique_genes:
        total_matched_unique_genes.write("gene,hit_count,qcov_average,pident_average,perfect_matches,gene_length,"
                                         "function,matched_organisms,perfect_matched_organisms\n")
        for match in final_output:
            if len(match) > 0:
                total_matched_unique_genes.write(f"{match[8]},{str(match[0])},{str(match[1])},{str(match[2])},"
                                                 f"{str(match[3])},{str(match[4])},{str(match[5])},"
                                                 f"{str(match[6]).replace(',', '.')},{str(match[7]).replace(',', '.')}\n")


PHENOTYPES = ["res"]
DRUGS = ["CIPRO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_uniques(DRUG, PHENOTYPE)
