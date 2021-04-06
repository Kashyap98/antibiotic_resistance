import csv
import multiprocessing
import os
import pandas as pd

from models import gene
from utils import dir_utils


class MultiProcessHandler:

    def __init__(self, max_processes: int, target, input_list: list):
        self.process_count = 0
        self.max_processes = max_processes
        self.target = target
        self.input_list = input_list
        self.pool: multiprocessing.Pool = None

    def start(self):
        self.pool = multiprocessing.Pool(self.max_processes)
        self.pool.starmap(self.target, self.input_list)


def get_op_phenotype(phenotype):
    op_phenotype = "res"
    if phenotype == "res":
        op_phenotype = "sus"

    return op_phenotype


def get_organisms_by_phenotype(organism_path):
    all_organisms_csv = list(csv.reader(open(organism_path), delimiter=','))
    all_organisms = set()
    for org in all_organisms_csv:
        all_organisms.add(org[0])
    return all_organisms


def get_organism_and_all_genes_from_folder_csv(folder_path: str, remove_hypothetical=False) -> dict:
    organisms = os.listdir(folder_path)
    all_genes = {}
    for organism in organisms:
        print(f"Gathering genes for organism: {organism}")
        organism_csv = pd.read_csv(os.path.join(folder_path, organism), header=0)
        organism = dir_utils.strip_extension(organism)
        organism_gene_names = list(organism_csv["gene_name"])
        organism_genes = gene.get_genes_from_list(organism=organism,
                                                  list_of_genes_names=organism_gene_names,
                                                  remove_hypothetical=remove_hypothetical)
        all_genes[organism] = organism_genes
    return all_genes


def get_all_genes_for_list_of_organisms(organisms: list, remove_hypothetical=False) -> dict:
    all_genes = {}
    for organism in organisms:
        all_genes[organism] = get_all_genes_for_organism(organism, remove_hypothetical)

    return all_genes


def get_all_genes_for_organism(organism: str, remove_hypothetical=False) -> list:
    organism_dir = dir_utils.OrganismDirs(organism)
    all_genes = gene.get_genes_from_list(organism=organism, list_of_genes_names=os.listdir(organism_dir.gene_folder),
                                         remove_hypothetical=remove_hypothetical)

    return all_genes


if __name__ == '__main__':
    genes = get_all_genes_for_organism("CP-49")
    for gene in genes:
        print(gene.description)
