from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
from models.Gene import Gene
import pandas as pd


def get_op_phenotype(phenotype):
    op_phenotype = "res"
    if phenotype == "res":
        op_phenotype = "sus"

    return op_phenotype


def get_organisms_genes_data_by_phenotype(drug, phenotype):
    all_organisms = pd.read_csv(os.path.join(os.getcwd(), "..", "sorted_data", drug, phenotype + ".csv"), header=None)
    all_organisms_list = list(all_organisms[0])
    recip_genes = pd.read_csv(os.path.join(os.getcwd(), "..", "sorted_data", drug, phenotype + "_RecipGenes.csv"))
    recip_genes_data = pd.read_csv(os.path.join(os.getcwd(), "..", "sorted_data", drug, phenotype + "_RecipBlastInfo.csv"))

    return all_organisms_list, recip_genes, recip_genes_data


def generate_recip_genes_mappings(recip_genes, original_organism, recip_genes_data):
    all_recip_genes = {}
    for recip_gene in recip_genes:
        all_recip_genes[recip_gene] = {}
        current_original_genes = all_recip_genes[recip_gene]
        current_original_genes[recip_gene] = (original_organism, 0)
        all_recip_genes[recip_gene] = current_original_genes

        gene_matches = recip_genes_data.loc[recip_genes_data[' feature_id'] == recip_gene]
        for index, row in gene_matches.iterrows():
            current_matches = all_recip_genes[row[" feature_id"]]
            current_matches[row[" recip_feature_id"]] = (row[" recip_organism"], 0)
            all_recip_genes[row[" feature_id"]] = current_matches
    return all_recip_genes


def get_records(all_genes, records, phenotype):
    for gene, info in all_genes.items():
        organism = info[0]
        new_gene = Gene(organism, gene.replace(".fasta", ""), get_info=True)
        records.append(SeqRecord(Seq(new_gene.aa_sequence, IUPAC.protein), id=f"{phenotype}_{organism}_{new_gene.name}"))
    return records


def generate_temp_total_fasta_file(drug, genes, opp_genes, gene):
    in_file = os.path.join(os.getcwd(), "..", "sorted_data", drug, "variations", "alignments", "temp_genes.fasta")
    out_file = os.path.join(os.getcwd(), "..", "sorted_data", drug, "variations", "alignments", f"{gene}")
    with open(in_file, "w") as temp_file:
        records = []
        records = get_records(genes, records, phenotype="res")
        records = get_records(opp_genes, records, phenotype="sus")
        SeqIO.write(records, temp_file, "fasta")

    dir_copy = os.getcwd()
    os.chdir(os.path.join(os.path.join(os.getcwd(), "..", "sorted_data", drug, "variations", "alignments")))
    os.system(f"wsl.exe mafft --quiet --auto --clustalout temp_genes.fasta > {gene}")
    os.chdir(dir_copy)
    msa = MultipleSeqAlignment(AlignIO.read(out_file, "clustal"))
    return msa
