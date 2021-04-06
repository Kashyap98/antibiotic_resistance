import os

from utils import dir_utils, gen_utils, output_util
from models import blast, gene

DRUG = "SULF"
PHENOTYPE = "res"


def check_genes_from_minimum_unique_set(organism, res_gene_list, res_organisms, sus_organisms, drug_dirs):
    all_results = []
    count = 0
    gene_list_length = len(res_gene_list)
    for res_gene in res_gene_list:
        count += 1
        combined_result = blast.CombinedResult(res_gene.name)
        print(f"Gene: {res_gene.description} | Count: {count} / {gene_list_length}")

        if not gene.check_if_unique(res_gene, sus_organisms):
            all_results.append(combined_result)
            continue

        for res_organism in res_organisms:
            if res_organism == organism:
                continue

            blast_data = blast.blast(res_gene, res_organism)
            if not blast_data:
                continue

            if not blast_data.is_homolog:
                continue

            reciprocal_blast = blast.blast(blast_data.blast_gene, blast_data.blast_gene.organism)
            if not reciprocal_blast:
                continue

            if reciprocal_blast.is_homolog:
                if reciprocal_blast.blast_gene.name == blast_data.target_gene.name:
                    combined_result.add_new_result(blast_data)
            else:
                continue

        all_results.append(combined_result)

    # output final unique genes for organism
    output_file = output_util.OutputFile(file_path=os.path.join(drug_dirs.reciprocal_res_genes, f"{organism}.csv"),
                                         header_list=all_results[0].header())
    for result in all_results:
        output_file.write_data_list_to_output_file(result.data())


def get_reciprocal_genes_for_organism(organism, res_genes, res_organisms, sus_organisms, drug_dirs):
    all_results = []
    count = 0
    gene_list_length = len(res_genes)
    for res_gene in res_genes:
        count += 1
        combined_result = blast.CombinedResult(res_gene.name)
        print(f"Organism: {organism} | Count: {count} / {gene_list_length}")

        if gene.check_if_unique(res_gene, sus_organisms):
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
                    if reciprocal_blast.blast_gene.name == blast_data.target_gene.name:
                        combined_result.add_new_result(blast_data)
                else:
                    continue

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
    process_data = []
    for organism in res_organisms:
        res_genes = gen_utils.get_all_genes_for_organism(organism, remove_hypothetical=True)
        process_data.append((organism, res_genes, res_organisms, sus_organisms, drug_dirs,))

    process_handler = gen_utils.MultiProcessHandler(max_processes=4, target=get_reciprocal_genes_for_organism,
                                                    input_list=process_data)
    process_handler.start()
