
from utils import dir_utils, gen_utils
from refactored import unique_genes

DRUG = "CEPRO"
PHENOTYPE = "res"
MAX_PROCESSES = 4


if __name__ == '__main__':
    drug_dirs = dir_utils.DrugDirs(DRUG, PHENOTYPE)
    dir_utils.generate_dir(drug_dirs.unique_res_genes)

    res_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.res_file)
    sus_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.sus_file)

    # compare 1 resistant organism's genes to ALL susceptible organism's genes
    unique_gene_input_list = []
    for res_organism in res_organisms:
        res_genes = gen_utils.get_all_genes_for_organism(res_organism)
        unique_gene_input_list.append((res_organism, res_genes, sus_organisms, drug_dirs,))

    process_handler = gen_utils.MultiProcessHandler(max_processes=MAX_PROCESSES,
                                                    target=unique_genes.get_unique_genes_for_organism,
                                                    input_list=unique_gene_input_list)
    process_handler.start()
    try:
        process_handler.pool.join()
    except PermissionError as e:
        pass

    # res_unique_genes_by_org = gen_utils.get_organism_and_all_genes_from_folder_csv(drug_dirs.unique_res_genes,
    #                                                                                remove_hypothetical=True)
    #
    # minimum_count_organism, minimum_count, organism_gene_list = None, 0, list()
    # for organism, gene_list in res_unique_genes_by_org.items():
    #     gene_count = len(gene_list)
    #     if minimum_count_organism is None:
    #         minimum_count_organism, minimum_count, organism_gene_list = organism, gene_count, gene_list
    #         continue
    #
    #     if gene_count <= minimum_count:
    #         minimum_count_organism, minimum_count, organism_gene_list = organism, gene_count, gene_list
    #
    # # split_genes = [organism_gene_list[i:i + MAX_PROCESSES] for i in range(0, len(organism_gene_list), MAX_PROCESSES)]
    # process_data = []
    # for organism, gene_list in res_unique_genes_by_org.items():
    #     process_data.append((organism, gene_list, res_organisms, sus_organisms, drug_dirs,))
    #
    # process_handler = gen_utils.MultiProcessHandler(max_processes=MAX_PROCESSES,
    #                                                 target=reciprocal_genes.check_genes_from_minimum_unique_set,
    #                                                 input_list=process_data)
    # process_handler.start()
    # process_handler.pool.join()
    # print(f"Minimum count organism = {minimum_count_organism}")
