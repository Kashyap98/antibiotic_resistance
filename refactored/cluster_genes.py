import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from models import gene

from utils import gen_utils, dir_utils, output_util

DRUG = "CEPRO"
PHENOTYPE = "res"


def _get_gene_clusters(folder_path: str) -> dict:
    clusters = os.listdir(folder_path)
    all_clusters = {}

    for cluster in clusters:
        print(f"Gathering genes for cluster: {cluster}")
        cluster_info = list(SeqIO.parse(open(os.path.join(folder_path, cluster), "r"), "fasta"))
        all_clusters[cluster] = cluster_info

    return all_clusters


def cluster_potential_unique_genes(drug_dirs: dir_utils.DrugDirs, write_sequence_file: bool = True):
    res_unique_genes_by_org = gen_utils.get_organism_and_all_genes_from_folder_csv(drug_dirs.unique_res_genes,
                                                                                   remove_hypothetical=True)
    dir_utils.generate_dir(drug_dirs.cluster_dir)
    if write_sequence_file:
        sequences = []
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

    os.system(f"{dir_utils.USEARCH} -cluster_fast {drug_dirs.sequences_to_cluster} -id 0.9"
              f" -centroids {drug_dirs.centroids} -uc {drug_dirs.cluster_uc} -clusters {drug_dirs.cluster_dir}/c_")


def _get_unique_clusters(cluster: str, unique_clusters_file: output_util.OutputFile, genes: list, sus_organisms: list):
    organisms = []
    gene_description = None
    is_unique = True
    for cluster_gene in genes:
        if gene_description is None:
            gene_description = cluster_gene.description

        gene_is_unique = gene.check_if_unique(cluster_gene, sus_organisms)
        organisms.append(cluster_gene.organism)
        if not gene_is_unique:
            is_unique = False
            break

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
    sus_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.sus_file)
    cluster_info: dict = _get_gene_clusters(drug_dirs.cluster_dir)
    cluster_info_file = output_util.OutputFile(file_path=os.path.join(drug_dirs.drug_dir, "cluster_info.csv"),
                                               header_list=["cluster", "gene", "count", "organisms"])
    unique_clusters_file = output_util.OutputFile(file_path=os.path.join(drug_dirs.drug_dir, "unique_clusters.csv"),
                                                  header_list=["cluster", "gene", "count", "organisms"])
    process_data = []
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

        cluster_info_file.write_data_list_to_output_file([cluster, first_gene.description, len(info), organisms])
        process_data.append((cluster, unique_clusters_file, genes, sus_organisms,))

    process_handler = gen_utils.MultiProcessHandler(max_processes=10, target=_get_unique_clusters,
                                                    input_list=process_data)
    process_handler.start()


if __name__ == '__main__':
    drug_dirs = dir_utils.DrugDirs(DRUG, PHENOTYPE)

    cluster_potential_unique_genes(drug_dirs, write_sequence_file=True)
    analyze_unique_gene_clusters(drug_dirs)
