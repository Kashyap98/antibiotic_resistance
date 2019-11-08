import os

import pandas as pd
from Bio import SeqIO
from utils.dir_utils import OrganismDirs


def format_gene_name(gene_name):
    if ".fasta" in gene_name:
        name = gene_name.replace(".fasta", "")
    else:
        name = gene_name

    return name


# contains all the properties we care about from the fasta files/csv full of gene information
class Gene:

    def __init__(self, organism, gene_name, get_info=False):

        organism_dirs = OrganismDirs(organism)
        self.organism = organism
        self.name = format_gene_name(gene_name)
        self.fasta_file = os.path.join(organism_dirs.gene_folder, self.name + ".fasta")
        self.txt_file = os.path.join(organism_dirs.gene_info_folder, self.name + ".txt")

        record = list(SeqIO.parse(self.fasta_file, "fasta"))
        self.nuc_sequence = record[0].seq
        self.length = len(self.nuc_sequence)
        self.upper_length = int(self.length * 1.10)
        self.lower_length = int(self.length * 0.90)

        if get_info:
            geneInfo = pd.read_csv(self.txt_file, error_bad_lines=False)
            geneData = list(geneInfo.iloc[:, 1])
            self.location = geneData[2]
            self.type_data = geneData[1]
            self.start = geneData[3]
            self.stop = geneData[4]
            self.strand = geneData[5]
            self.function_data = geneData[6]
            self.aliases = geneData[7]
            self.figfam = geneData[8]
            self.evidence_codes = geneData[9]
            self.aa_sequence = geneData[11]
