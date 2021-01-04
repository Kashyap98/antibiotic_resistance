import multiprocessing
import os

from utils import dir_utils, gen_utils, output_util
from models import blast

DRUG = "CIPRO"
PHENOTYPE = "res"


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

        # replace all res_genes with unique ones
        print(f"Genes to check: {len(unique_genes)}")
        res_genes = unique_genes

    # output final unique genes for organism
    output_file = output_util.OutputFile(file_path=os.path.join(drug_dirs.unique_res_genes, f"{res_organism}.csv"),
                                         header_list=["gene_name", "gene_info"])
    for res_gene in res_genes:
        output_file.write_data_list_to_output_file([res_gene.name, res_gene.description])


if __name__ == '__main__':
    drug_dirs = dir_utils.DrugDirs(DRUG, PHENOTYPE)
    dir_utils.generate_dir(drug_dirs.unique_res_genes)

    res_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.res_file)
    sus_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.sus_file)

    # compare 1 resistant organism's genes to ALL susceptible organism's genes
    for i in range(6):
        res_organism = res_organisms[i]
        res_genes = gen_utils.get_all_genes_for_organism(res_organism)
        process = multiprocessing.Process(target=get_unique_genes_for_organism,
                                          args=(res_organism, res_genes, sus_organisms, drug_dirs,))
        process.start()
