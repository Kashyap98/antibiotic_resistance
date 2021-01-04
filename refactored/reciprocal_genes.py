import multiprocessing
import os

from utils import dir_utils, gen_utils, output_util
from models import blast, gene

DRUG = "CIPRO"
PHENOTYPE = "res"


def check_if_unique(res_gene: gene.Gene, sus_organisms: list) -> bool:
    for sus_organism in sus_organisms:
        sus_blast_data = blast.blast(res_gene, sus_organism)
        if sus_blast_data:
            if sus_blast_data.is_perfect_match:
                return False

    return True


def get_reciprocal_genes_for_organism(organism, res_genes, res_organisms, sus_organisms, drug_dirs):
    all_results = []
    count = 0
    gene_list_length = len(res_genes)
    for res_gene in res_genes:
        count += 1
        combined_result = blast.CombinedResult(res_gene.name, res_gene.description)
        print(f"Organism: {organism} | Count: {count} / {gene_list_length}")

        if check_if_unique(res_gene, sus_organisms):
            for res_organism in res_organisms:
                if res_organism == organism:
                    continue

                blast_data = blast.blast(res_gene, res_organism)
                if not blast_data:
                    continue

                if not blast_data.is_homolog:
                    continue

                reciprocal_blast = blast.blast(blast_data.blast_gene, organism)
                if not reciprocal_blast:
                    continue

                if reciprocal_blast.is_homolog:
                    combined_result.add_new_result(blast_data)
                else:
                    continue

        if combined_result.results > 0:
            all_results.append(combined_result)

    # output final unique genes for organism
    output_file = output_util.OutputFile(file_path=os.path.join(drug_dirs.reciprocal_res_genes, f"{organism}.csv"),
                                         header_list=all_results[0].header())
    for result in all_results:
        output_file.write_data_list_to_output_file(result.data())


if __name__ == '__main__':
    drug_dirs = dir_utils.DrugDirs(DRUG, PHENOTYPE)
    dir_utils.generate_dir(drug_dirs.reciprocal_res_genes)

    res_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.res_file)
    sus_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.sus_file)
    for i in range(6):
        organism = res_organisms[i]
        res_genes = gen_utils.get_all_genes_for_organism(organism, remove_hypothetical=True)
        process = multiprocessing.Process(target=get_reciprocal_genes_for_organism,
                                          args=(organism, res_genes, res_organisms, sus_organisms, drug_dirs,))
        process.start()

    # res_unique_genes_by_org = gen_utils.get_organism_and_all_genes_from_folder_csv(drug_dirs.unique_res_genes,
    #                                                                                remove_hypothetical=True)

    # minimum_count_organism, minimum_count, organism_gene_list = None, 0, list()
    # for organism, gene_list in res_unique_genes_by_org.items():
    #     gene_count = len(gene_list)
    #     if minimum_count_organism is None:
    #         minimum_count_organism, minimum_count, organism_gene_list = organism, gene_count, gene_list
    #         continue
    #
    #     if gene_count <= minimum_count:
    #         minimum_count_organism, minimum_count, organism_gene_list = organism, gene_count, gene_list
