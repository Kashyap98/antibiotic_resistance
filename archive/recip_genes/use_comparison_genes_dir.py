import os

import pandas as pd

from utils.dir_utils import DrugDirs
from utils import output_util
import utils.gen_utils as gen_utils


def get_uniques(drug, phenotype):
    # create output file and write the headers
    drug_dirs = DrugDirs(drug, phenotype)
    res_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.res_file)
    genes_by_org = {}
    no_perfect_matches = {}
    potential_no_perfect_matches = {}
    thrown_out = {}
    organism = ""
    potentials_list = []
    for organism in res_organisms:
        file_path = os.path.join(drug_dirs.comparison_dir, f"{organism}_genes.csv")
        genes_by_org[organism] = pd.read_csv(file_path, header=0)

    count = len(genes_by_org[organism]["gene"])

    for i in range(count):
        unique_orgs = set()
        potential_orgs = []
        gene_info_row = {}
        gene_rows = {}
        perfect_match_tracker = []
        no_perfect_match_found = False
        for organism, organism_info in genes_by_org.items():
            gene_info_row = organism_info.loc[i]

            if "gyrase" in gene_info_row["gene_description"]:
                print("found")
            gene_rows[organism] = gene_info_row
            perfect_match_count = int(gene_info_row["perfect_match_count"])
            perfect_match_tracker.append(perfect_match_count)
            if perfect_match_count == 0:
                unique_orgs.add(organism)
                no_perfect_match_found = True
            else:
                potential_orgs.append(organism)

        potentials_list.append((gene_info_row["gene_description"], unique_orgs))

        if len(unique_orgs) == len(res_organisms):
            no_perfect_matches[gene_info_row["gene"]] = [gene_info_row["gene_description"]]
        else:
            output_row = [gene_info_row["gene_description"]]
            output_row.extend(perfect_match_tracker)
            if "subunit" in gene_info_row["gene_description"]:
                # for org, info_row in gene_rows.items():
                if no_perfect_match_found:
                    potential_no_perfect_matches[f"{gene_info_row['gene']}"] = output_row
            else:
                thrown_out[gene_info_row["gene"]] = output_row

    no_perfect_match_file = output_util.OutputFile(drug_dirs.res_to_sus_no_perfect_matches,
                                                   ["gene", "gene_description"])
    no_perfect_match_file.write_data_dict_to_output_file(no_perfect_matches)

    output_header = ["gene", "gene_description"]
    output_header.extend(res_organisms)

    potential_file = output_util.OutputFile(drug_dirs.res_to_sus_potential, output_header)
    potential_file.write_data_dict_to_output_file(potential_no_perfect_matches)

    thrown_out_file = output_util.OutputFile(drug_dirs.res_to_sus_thrown_out, output_header)
    thrown_out_file.write_data_dict_to_output_file(thrown_out)

    # potential_suspects = output_util.OutputFile(drug_dirs.potential_suspects, ["gene_1", "gene_2", "gene_3"])
    #
    # for i in range(len(potentials_list)):
    #     row = potentials_list[i]
    #     gene_1, gene_1_unique = row[0], row[1]
    #     print(f"Progress: {i}")
    #     for j in range(1, len(potentials_list)):
    #         compare_row = potentials_list[j]
    #         gene_2, gene_2_unique = compare_row[0], compare_row[1]
    #         temp_new_set = set()
    #         temp_new_set.update(gene_1_unique)
    #         temp_new_set.update(gene_2_unique)
    #         if len(temp_new_set) == len(res_organisms):
    #             potential_suspects.write_data_list_to_output_file([gene_1, gene_2])
    #             continue
    #
    #         for z in range(1, len(potentials_list)):
    #             new_compare_row = potentials_list[z]
    #             gene_3, gene_3_unique = new_compare_row[0], new_compare_row[1]
    #             newest_temp_new_set = set()
    #             newest_temp_new_set.update(temp_new_set)
    #             newest_temp_new_set.update(gene_3_unique)
    #             if len(newest_temp_new_set) == len(res_organisms):
    #                 potential_suspects.write_data_list_to_output_file([gene_1, gene_2, gene_3])


PHENOTYPES = ["res"]
DRUGS = ["SULF"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_uniques(DRUG, PHENOTYPE)
