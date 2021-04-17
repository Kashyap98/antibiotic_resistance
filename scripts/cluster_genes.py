import os
from typing import List, Dict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from models import gene
from models.gene import Gene
from utils import gen_utils, dir_utils, output_util

DRUG = "CEPRO"
PHENOTYPE = "res"
MAX_PROCESSES = 10


def _get_gene_clusters(folder_path: str) -> Dict[str, List[SeqRecord]]:
    """
    Get all of the gene clusters for unique genes.
    """
    clusters = os.listdir(folder_path)
    all_clusters = {}

    for cluster in clusters:
        print(f"Gathering genes for cluster: {cluster}")
        cluster_info = gen_utils.get_list_of_genes_from_fasta_file(os.path.join(folder_path, cluster))
        all_clusters[cluster] = cluster_info

    return all_clusters


def cluster_potential_unique_genes(drug_dirs: dir_utils.DrugDirs, write_sequence_file: bool = True):
    """
    Gather all the unique genes by organism. Write sequence record files if specified and use usearch to cluster the
    genes into groups based on their identity towards one another.
    """
    res_unique_genes_by_org = gen_utils.get_organism_and_all_genes_from_folder_csv(drug_dirs.unique_res_genes,
                                                                                   remove_hypothetical=True)
    # create the cluster directory
    dir_utils.generate_dir(drug_dirs.cluster_dir)
    if write_sequence_file:
        sequences = []
        # write sequence record files for each gene in the unique gene set
        for organism, gene_list in res_unique_genes_by_org.items():
            print(f"Processing {organism=}")
            for gene in gene_list:
                new_record = SeqRecord(seq=gene.nuc_sequence,
                                       id=gene.name,
                                       name=gene.name,
                                       description=f"{gene.description}~{gene.organism}")
                sequences.append(new_record)

        with open(drug_dirs.sequences_to_cluster, "w") as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")

    # run the usearch cluster organism
    os.system(f"{dir_utils.USEARCH} -cluster_fast {drug_dirs.sequences_to_cluster} -id 0.9"
              f" -centroids {drug_dirs.centroids} -uc {drug_dirs.cluster_uc} -clusters {drug_dirs.cluster_dir}/c_")


def _get_unique_clusters(cluster: str, unique_clusters_file: output_util.OutputFile,
                         genes: List[Gene], sus_organisms: List[str]):
    """
    This function takes in a cluster and its genes. It compares this clusters genes to the susceptible organisms to
    ensure that the genes in the cluster are unique. This then outputs to unique_clusters_file.
    """
    organisms = []
    gene_description = None
    is_unique = True
    # go though every gene in the cluster
    for cluster_gene in genes:
        if gene_description is None:
            gene_description = cluster_gene.description

        # check if the gene is unique
        gene_is_unique = gene.check_if_unique(cluster_gene, sus_organisms)
        organisms.append(cluster_gene.organism)
        if not gene_is_unique:
            is_unique = False
            break

    # write the unique gene to the file
    if is_unique:
        written = False
        while not written:
            try:
                unique_clusters_file.write_data_list_to_output_file(
                    [cluster, gene_description, len(genes), organisms])
                written = True
            except:
                pass


def analyze_unique_gene_clusters(drug_dirs: dir_utils.DrugDirs):
    """
    Parent function for gathering and analyzing gene clusters for resistant organisms.
    """
    sus_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.sus_file)
    cluster_info: dict = _get_gene_clusters(drug_dirs.cluster_dir)
    cluster_info_file = output_util.OutputFile(file_path=drug_dirs.cluster_info,
                                               header_list=["cluster", "gene", "count", "organisms"])
    unique_clusters_file = output_util.OutputFile(file_path=drug_dirs.unique_clusters,
                                                  header_list=["cluster", "gene", "count", "organisms"])
    process_data = []
    # go though each cluster as well the info related ot the cluster
    for cluster, info in cluster_info.items():
        print(f"Processing {cluster=}")

        first_gene: SeqRecord = info[0]
        organisms = []
        genes = []
        for gene_info in info:
            gene_description = gene_info.description
            gene_organism = gene_description.split("~")[1]
            organisms.append(gene_organism)
            genes.append(gene.Gene(organism=gene_organism, gene_name=gene_info.name))

        # write data for each cluster and prepare processes for analyzing if the clusters are unique
        cluster_info_file.write_data_list_to_output_file([cluster, first_gene.description, len(info), organisms])
        process_data.append((cluster, unique_clusters_file, genes, sus_organisms,))

    process_handler = gen_utils.MultiProcessHandler(max_processes=MAX_PROCESSES, target=_get_unique_clusters,
                                                    input_list=process_data)
    process_handler.start()


if __name__ == '__main__':
    drug_dirs_parent = dir_utils.DrugDirs(DRUG, PHENOTYPE)

    cluster_potential_unique_genes(drug_dirs_parent, write_sequence_file=True)
    analyze_unique_gene_clusters(drug_dirs_parent)
