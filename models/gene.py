import os

import pandas as pd
from Bio import SeqIO
from utils import dir_utils


def format_gene_name(gene_name):
    if ".fasta" in gene_name:
        name = gene_name.replace(".fasta", "")
    else:
        name = gene_name

    return name


def get_genes_from_list(organism: str, list_of_genes_names: list, remove_hypothetical=False):
    genes = []
    for organism_gene in list_of_genes_names:
        gene = Gene(organism, organism_gene, get_info=False)
        if remove_hypothetical:
            if "hypothetical" in gene.description:
                continue

        genes.append(gene)

    return genes


# contains all the properties we care about from the fasta files/csv full of gene information
class Gene:

    def __init__(self, organism, gene_name, get_info=False):

        organism_dirs = dir_utils.OrganismDirs(organism)
        self.organism = organism
        self.name = format_gene_name(gene_name)
        self.fasta_file = os.path.join(organism_dirs.gene_folder, f"{self.name}.fasta")
        self.txt_file = os.path.join(organism_dirs.gene_info_folder, f"{self.name}.txt")

        record = list(SeqIO.parse(self.fasta_file, "fasta"))
        self.nuc_sequence = record[0].seq
        self.description = record[0].description
        self.length = len(self.nuc_sequence)
        self.upper_length = int(self.length * 1.10)
        self.lower_length = int(self.length * 0.90)
        self.function_data = ""
        self.match_list = []

        if get_info:
            gene_info = pd.read_csv(self.txt_file, error_bad_lines=False)
            gene_data = list(gene_info.iloc[:, 1])
            self.location = gene_data[2]
            self.type_data = gene_data[1]
            self.start = gene_data[3]
            self.stop = gene_data[4]
            self.strand = gene_data[5]
            self.function_data = gene_data[6]
            self.aliases = gene_data[7]
            self.figfam = gene_data[8]
            self.evidence_codes = gene_data[9]
            self.aa_sequence = gene_data[11]
