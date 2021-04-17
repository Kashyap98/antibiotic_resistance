import os
import re
from typing import Optional, List

from Bio import SeqIO

from models import gene
from utils import dir_utils, gen_utils, output_util
import pandas as pd

DRUG = "CIPRO"
PHENOTYPE = "res"
MAX_PROCESSES = 4


def check_unique_clusters_for_genes_of_interest(drug_dirs: dir_utils.DrugDirs,
                                                organism_file_path: str,
                                                genes_to_collect: Optional[List[str]] = None,
                                                genes_to_filter: Optional[List[str]] = None):

    sus_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.sus_file)

    # get all unique clusters
    cluster_data = pd.read_csv(drug_dirs.unique_clusters, header=0)

    # filter clusters for genes of interest, removing any that should be filtered
    filtered_clusters = []
    for index, row in cluster_data.iterrows():
        gene_name = row["gene"]
        should_filter = gen_utils.check_if_gene_in_keyword_list(gene_name, genes_to_filter)
        # not filtering cluster
        if not should_filter:
            should_collect = gen_utils.check_if_gene_in_keyword_list(gene_name, genes_to_collect)

            # gene is in keyword list
            if should_collect:
                filtered_clusters.append(row)

    # get all genes of new organism
    organism_gene_list = gen_utils.get_list_of_genes_from_fasta_file(organism_file_path)

    # create fasta file and dir for each file (needed for blast)
    gene_object_list = []
    dir_utils.generate_dir(drug_dirs.new_organism_dir, overwrite_dir=True)
    for organism_gene in organism_gene_list:
        organism_gene.id = re.sub(r'[\\/*?:"<>|]', "", organism_gene.id).replace(".", "_")
        with open(os.path.join(drug_dirs.new_organism_dir, f"{organism_gene.id}.fasta"), "w") as output_handle:
            SeqIO.write(organism_gene, output_handle, "fasta")

        gene_object = gene.Gene("new_organism", organism_gene.id, new_organism=True)
        gene_object.description = organism_gene.description
        gene_object_list.append(gene_object)

    # check if gene is unique to resistant group
    unique_to_resistant: List[gene.Gene] = []
    for organism_gene in gene_object_list:
        print(organism_gene.description)
        is_unique = gene.check_if_unique(organism_gene, sus_organisms)
        if is_unique:
            unique_to_resistant.append(organism_gene)

    # check if gene belongs in cluster?

    # output genes that are unique/belong to cluster
    output_file = output_util.OutputFile(file_path=drug_dirs.investigated_unique_genes,
                                         header_list=["gene_name"])
    for result in unique_to_resistant:
        output_file.write_data_list_to_output_file(result.description)


if __name__ == '__main__':
    drug_genes_to_investigate = ["dihydrofolate", "dihydromonapterin"]
    drug_genes_to_filter = ["resistant"]
    drug_dirs_parent = dir_utils.DrugDirs(DRUG, PHENOTYPE)
    organism_path = os.path.join(dir_utils.MAIN_DIR, "ncbi_annotations.faa")

    check_unique_clusters_for_genes_of_interest(drug_dirs_parent,
                                                organism_file_path=organism_path,
                                                genes_to_collect=drug_genes_to_investigate,
                                                genes_to_filter=drug_genes_to_filter)


