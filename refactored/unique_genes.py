from collections import defaultdict
import os

from utils import dir_utils, gen_utils, output_util
from models import blast
from models import gene as gene_utils

DRUG = "SULF"
PHENOTYPE = "res"
MAX_PROCESSES = 16


def gather_filtered_potential_unique_genes(drug_dirs: dir_utils.DrugDirs, genes_to_collect: list,
                                           genes_to_filter: list, write_output: bool = False) -> dict:

    final_gene_output = {}
    potential_combinations = gather_potential_unique_combinations(drug_dirs, write_output)

    # check all possible unique genes
    for gene, gene_list in potential_combinations.items():
        for gene_to_collect in genes_to_collect:

            # make sure it is a gene we are interested in
            if gene_to_collect in gene:

                # check if gene should be filtered
                should_filter_gene = False
                for gene_to_filter in genes_to_filter:
                    if gene_to_filter in gene:
                        should_filter_gene = True
                        break

                if should_filter_gene:
                    continue

                final_gene_output[gene] = gene_list

    return final_gene_output


def investigate_potential_unique_combinations(drug_dirs: dir_utils.DrugDirs, genes_to_collect: list,
                                              genes_to_filter: list, write_output: bool = False):

    filtered_unique_combinations: dict = gather_filtered_potential_unique_genes(drug_dirs, genes_to_collect,
                                                                                genes_to_filter, write_output)
    res_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.res_file)

    final_gene_output = []
    for gene, gene_list in filtered_unique_combinations.items():
        unique_organisms = set(gene_utils.get_organisms_from_list_of_genes(gene_list))
        not_unique_organisms = res_organisms - unique_organisms

        print(f"Checking gene: {gene}")
        for potential_gene in gene_list:
            combined_result = blast.UniqueGeneCompareResult(potential_gene, unique_group=unique_organisms,
                                                            not_unique_group=not_unique_organisms)

            # check all other unique organisms for a perfect match
            for unique_organism in unique_organisms:
                # we don't want to compare the gene to itself.
                if unique_organism == potential_gene.organism:
                    continue

                unique_blast = blast.blast(potential_gene, unique_organism)
                if not unique_blast:
                    continue

                combined_result.add_new_result(unique_blast)

            # check all non unique organisms for a perfect match
            for not_unique_organism in not_unique_organisms:
                not_unique_blast = blast.blast(potential_gene, not_unique_organism)
                if not not_unique_blast:
                    continue

                combined_result.add_new_sus_result(not_unique_blast)

            final_gene_output.append(combined_result)
            print(f"Finished checking gene: {gene}")

    if write_output:
        # output final unique genes for organism
        output_file = output_util.OutputFile(file_path=os.path.join(drug_dirs.drug_dir,
                                                                    f"investigated_unique_genes.csv"),
                                             header_list=final_gene_output[0].header())
        for result in final_gene_output:
            output_file.write_data_list_to_output_file(result.data())

    return final_gene_output


def gather_potential_unique_combinations(drug_dirs: dir_utils.DrugDirs, write_output: bool = False):

    res_unique_genes_by_org = gen_utils.get_organism_and_all_genes_from_folder_csv(drug_dirs.unique_res_genes,
                                                                                   remove_hypothetical=True)
    final_gene_output = defaultdict(list)
    for organism, gene_list in res_unique_genes_by_org.items():
        print(f"Gathering unique genes from organism: {organism}")
        for gene in gene_list:
            if gene.description in final_gene_output:
                if organism not in final_gene_output[gene.description]:
                    final_gene_output[gene.description].append(gene)
            else:
                final_gene_output[gene.description].append(gene)

    if write_output:
        output_file = output_util.OutputFile(file_path=os.path.join(drug_dirs.drug_dir, f"potential_uniques.csv"),
                                             header_list=["gene", "res_organisms"])
        for gene, gene_list in final_gene_output.items():
            org_list = gene_utils.get_organisms_from_list_of_genes(gene_list)
            output_file.write_data_list_to_output_file([gene, org_list])

    return final_gene_output


def get_unique_genes_for_organism(res_organism, res_genes, sus_organisms, drug_dirs):
    print(f"Genes to check: {len(res_genes)}")
    for sus_organism in sus_organisms:
        unique_genes = []

        for res_gene in res_genes:
            blast_data = blast.blast(res_gene, sus_organism)

            # only add genes that did not perfectly match
            if not blast_data:
                unique_genes.append(res_gene)
                continue
            if not blast_data.is_perfect_match:
                unique_genes.append(res_gene)

        # replace all res_genes with only unique_genes
        print(f"Genes to check: {len(unique_genes)}")
        res_genes = unique_genes

    # output final unique genes for organism
    output_file = output_util.OutputFile(file_path=os.path.join(drug_dirs.unique_res_genes, f"{res_organism}.csv"),
                                         header_list=["gene_name", "gene_info"])
    for res_gene in res_genes:
        output_file.write_data_list_to_output_file([res_gene.name, res_gene.description])


if __name__ == '__main__':
    drug_dirs = dir_utils.DrugDirs(DRUG, PHENOTYPE)

    # gather_potential_unique_combinations(drug_dirs)

    # sulf genes
    genes_to_investigate = ["dihydrofolate", "dihydromonapterin"]
    genes_to_filter = ["resistant"]

    # cipro genes
    # genes_to_investigate = ["gyrase", "topoisomerase"]
    # genes_to_filter = ["type", "III", "inhibitor"]
    investigate_potential_unique_combinations(drug_dirs, genes_to_investigate, genes_to_filter, write_output=True)

    # dir_utils.generate_dir(drug_dirs.unique_res_genes)
    #
    # res_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.res_file)
    # sus_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.sus_file)
    # process_data = []
    #
    # # compare 1 resistant organism's genes to ALL susceptible organism's genes
    # for res_organism in res_organisms:
    #     res_genes = gen_utils.get_all_genes_for_organism(res_organism)
    #     process_data.append((res_organism, res_genes, sus_organisms, drug_dirs,))
    #
    # process_handler = gen_utils.MultiProcessHandler(max_processes=MAX_PROCESSES,
    #                                                 target=get_unique_genes_for_organism,
    #                                                 input_list=process_data)
    # process_handler.start()
