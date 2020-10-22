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
    drug_dirs.set_opposite_phenotype_file(op_phenotype)
    dir_utils.generate_dir(gene_finder_dir)

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
            if "beta-lactamase" in temp_gene.description:
                genes_found.append(temp_gene)
                count += 1
                print(f"{count} / {len(all_orgs)}")
                break

    final_output = []
    count = 0
    for res_gene in genes_found:
        gene_output = []

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
                if len(gene_output) == 0:
                    if res_gene.organism in op_orgs:
                        resistance = "sus"
                    else:
                        resistance = "res"
                    gene_output = [0, 0, 0, 0, res_gene.length,
                                   res_gene.nuc_sequence, res_gene.description, resistance, [],
                                   [], [],
                                   [], res_gene.organism]
                continue

            if len(gene_output) != 0:
                old_data = gene_output.copy()
                new_count = old_data[0] + 1
                new_qcov = (old_data[1] + blast_data.qcov) / 2
                new_pident = (old_data[2] + blast_data.pident) / 2
                resistance = old_data[7]
                new_perfect_matches = old_data[3]
                new_res_match_organisms = old_data[8]
                new_sus_match_organisms = old_data[9]
                new_res_perfect_match_organisms = old_data[10]
                new_sus_perfect_match_organisms = old_data[11]

                if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                    new_perfect_matches = old_data[3] + 1

                    if compare_org_resistance == "res":
                        new_res_perfect_match_organisms.append(compare_org)
                    else:
                        new_sus_perfect_match_organisms.append(compare_org)

                else:

                    if compare_org_resistance == "res":
                        new_res_match_organisms.append(compare_org)
                    else:
                        new_sus_match_organisms.append(compare_org)

                gene_output = [new_count, new_qcov, new_pident, new_perfect_matches, res_gene.length,
                               res_gene.nuc_sequence, res_gene.description, resistance, new_res_match_organisms,
                               new_sus_match_organisms, new_res_perfect_match_organisms,
                               new_sus_perfect_match_organisms, res_gene.organism]

            else:
                hit_count = 1

                if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                    perfect_matches = 1
                    res_match_organisms = []
                    sus_match_organisms = []
                    if compare_org_resistance == "res":
                        res_perfect_match_organisms = [compare_org]
                        sus_perfect_match_organisms = []
                    else:
                        res_perfect_match_organisms = []
                        sus_perfect_match_organisms = [compare_org]
                else:
                    perfect_matches = 0
                    res_perfect_match_organisms = []
                    sus_perfect_match_organisms = []
                    if compare_org_resistance == "res":
                        res_match_organisms = [compare_org]
                        sus_match_organisms = []
                    else:
                        res_match_organisms = []
                        sus_match_organisms = [compare_org]
                if res_gene.organism in op_orgs:
                    resistance = "sus"
                else:
                    resistance = "res"

                gene_output = [hit_count, blast_data.qcov, blast_data.pident, perfect_matches, res_gene.length,
                               res_gene.nuc_sequence, res_gene.description, resistance, res_match_organisms,
                               sus_match_organisms, res_perfect_match_organisms,
                               sus_perfect_match_organisms, res_gene.organism]

        final_output.append(gene_output)

        count += 1
        print(f"{count} / {len(genes_found)}")

    with open(os.path.join(gene_finder_dir, f"{drug} - beta-lactamase.csv"), "w") as total_matched_unique_genes:
        total_matched_unique_genes.write("organism,resistance,gene,hit_count,qcov_average,pident_average,perfect_matches,gene_length,"
                                         "aa_seq,res_matched,sus_matched,res_perfect_matched,sus_perfect_matched\n")
        for match in final_output:
            if len(match) > 0:
                total_matched_unique_genes.write(f"{match[12]},{match[7]},{match[6]},{str(match[0])},{str(match[1])},{str(match[2])},"
                                                 f"{str(match[3])},{str(match[4])},{str(match[5])},"
                                                 f"{str(match[8]).replace(',', '.')},{str(match[9]).replace(',', '.')},"
                                                 f"{str(match[10]).replace(',', '.')},{str(match[11]).replace(',', '.')}\n")


PHENOTYPES = ["res"]
DRUGS = ["AMOXO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_uniques(DRUG, PHENOTYPE)
