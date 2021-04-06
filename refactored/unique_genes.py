from collections import defaultdict
import os
from typing import List, Optional, Dict, Any

from models.gene import Gene
from utils import dir_utils, gen_utils, output_util
from models import blast
from models import gene as gene_utils

DRUG = "SULF"
PHENOTYPE = "res"
MAX_PROCESSES = 16


def gather_filtered_potential_unique_genes(drug_dirs: dir_utils.DrugDirs,
                                           genes_to_collect: Optional[List[str]] = None,
                                           genes_to_filter: Optional[List[str]] = None,
                                           write_output: bool = False) -> Dict[Any, List[Gene]]:
    """
    Takes in a dictionary of key gene with a list of genes that are related and unique to the key gene.
    This function will determine if any are needed to be filtered for collection or filtered for removal. Otherwise,
    all genes will eb taken.
    """
    final_gene_output = {}
    potential_combinations = gather_potential_unique_combinations(drug_dirs, write_output)

    # check all possible unique genes
    for gene, gene_list in potential_combinations.items():

        if len(genes_to_collect) > 0:
            should_collect_gene = False

            # make sure it is a gene we are interested in
            for gene_to_collect in genes_to_collect:
                if gene_to_collect in gene:
                    should_collect_gene = True
                    break
        else:
            should_collect_gene = True

        if should_collect_gene:
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


def investigate_potential_unique_combinations(drug_dirs: dir_utils.DrugDirs,
                                              genes_to_collect: Optional[List[str]] = None,
                                              genes_to_filter: Optional[List[str]] = None,
                                              write_output: bool = False):
    """
    This function will output filtered genes that are unique to the resistant set of organisms.
    """
    # determine if any of the genes need to be removed/filtered
    if genes_to_collect is None:
        genes_to_collect = []
    if genes_to_filter is None:
        genes_to_filter = []

    filtered_unique_combinations: dict = gather_filtered_potential_unique_genes(drug_dirs, genes_to_collect,
                                                                                genes_to_filter, write_output)
    res_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.res_file)

    final_gene_output = []
    # go through each gene in the filtered set of genes
    for gene, gene_list in filtered_unique_combinations.items():
        # get the unique organisms for the set of genes
        unique_organisms = set(gene_utils.get_organisms_from_list_of_genes(gene_list))
        # compare against other resistant organisms that do not contain unique copies of the gene of interest
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
        output_file = output_util.OutputFile(file_path=drug_dirs.investigated_unique_genes,
                                             header_list=final_gene_output[0].header())
        for result in final_gene_output:
            output_file.write_data_list_to_output_file(result.data())

    return final_gene_output


def gather_potential_unique_combinations(drug_dirs: dir_utils.DrugDirs,
                                         write_output: bool = False) -> Dict[Gene, List[Gene]]:
    """
    Takes all the unique genes for the resistant set and outputs the list with descriptions for each gene as well
    as the resistant organisms that are uniquely matched.
    """
    # gather all genes organized by organism
    res_unique_genes_by_org = gen_utils.get_organism_and_all_genes_from_folder_csv(drug_dirs.unique_res_genes,
                                                                                   remove_hypothetical=True)
    final_gene_output = defaultdict(list)
    for organism, gene_list in res_unique_genes_by_org.items():
        print(f"Gathering unique genes from organism: {organism}")
        for gene in gene_list:
            # check if gene name is in output
            if gene.description in final_gene_output:
                # add it if there is not already a copy in the output
                if organism not in final_gene_output[gene.description]:
                    final_gene_output[gene.description].append(gene)
            else:
                final_gene_output[gene.description].append(gene)

    # create output file for potential unique genes
    if write_output:
        output_file = output_util.OutputFile(file_path=drug_dirs.potential_uniques,
                                             header_list=["gene", "res_organisms"])
        for gene, gene_list in final_gene_output.items():
            org_list = gene_utils.get_organisms_from_list_of_genes(gene_list)
            output_file.write_data_list_to_output_file([gene, org_list])

    return final_gene_output


def get_unique_genes_for_organism(res_organism: str, res_genes: List[Gene],
                                  sus_organisms: List[str],  drug_dirs: dir_utils.DrugDirs):
    """
    This function compares the genes of one resistant organism to all susceptible organisms and their genes.
    The list of unique genes is then outputted.
    """
    print(f"Genes to check: {len(res_genes)}")
    for sus_organism in sus_organisms:
        unique_genes = []

        # go through every resistant gene for the organism
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
    drug_dirs_parent = dir_utils.DrugDirs(DRUG, PHENOTYPE)

    # gather_potential_unique_combinations(drug_dirs)

    # sulf genes
    drug_genes_to_investigate = ["dihydrofolate", "dihydromonapterin"]
    drug_genes_to_filter = ["resistant"]

    # cipro genes
    # drug_genes_to_investigate = ["gyrase", "topoisomerase"]
    # drug_genes_to_filter = ["type", "III", "inhibitor"]
    investigate_potential_unique_combinations(drug_dirs_parent, drug_genes_to_investigate, drug_genes_to_filter,
                                              write_output=True)

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
