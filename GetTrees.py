import os

import pandas as pd
from Bio import SeqIO, AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from models.Gene import Gene


def get_organisms_genes_data_by_phenotype(drug, phenotype):
    all_organisms = pd.read_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + ".csv"), header=None)
    all_organisms_list = list(all_organisms[0])
    recip_genes = pd.read_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_RecipGenes.csv"))
    recip_genes_data = pd.read_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_RecipBlastInfo.csv"))

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


def generate_temp_total_fasta_file(genes, opp_genes):
    in_file = os.path.join(os.getcwd(), "sorted_data", drug, "variations", "temp_genes.fasta")
    out_file = os.path.join(os.getcwd(), "sorted_data", drug, "variations", "aligned_temp_genes.fasta")
    with open(in_file, "w") as temp_file:
        records = []
        records = get_records(genes, records, phenotype="res")
        records = get_records(opp_genes, records, phenotype="sus")
        SeqIO.write(records, temp_file, "fasta")

    dir_copy = os.getcwd()
    os.chdir(os.path.join(os.path.join(os.getcwd(), "sorted_data", drug, "variations")))
    os.system("wsl.exe mafft --quiet --auto --clustalout temp_genes.fasta > aligned_temp_genes.fasta")
    os.chdir(dir_copy)
    msa = MultipleSeqAlignment(AlignIO.read(out_file, "clustal"))
    return msa


constructor = DistanceTreeConstructor()
calculator = DistanceCalculator('identity')


def get_trees(drug, phenotype):
    op_phenotype = "res"
    if phenotype == "res":
        op_phenotype = "sus"

    # get all the organism data.
    df = pd.read_csv(os.path.join(os.getcwd(), "sorted_data", drug, "variations", f"{phenotype}_both_var.csv"), index_col=0)
    all_organisms_list, recip_genes_df, recip_genes_data = get_organisms_genes_data_by_phenotype(drug, phenotype)
    recip_genes = list(recip_genes_df["0"])
    original_organism = all_organisms_list.pop(0)

    opposite_organisms_list, opposite_genes_df, opposite_genes_data = get_organisms_genes_data_by_phenotype(drug,
                                                                                                            op_phenotype)
    opposite_genes = list(opposite_genes_df["0"])
    opposite_organism = opposite_organisms_list.pop(0)
    all_recip_genes = generate_recip_genes_mappings(recip_genes, original_organism, recip_genes_data)
    all_opposite_genes = generate_recip_genes_mappings(opposite_genes, opposite_organism, opposite_genes_data)

    for gene in df["res_gene"]:
        if gene in all_recip_genes:
            opp_gene_name = df["sus_gene"].loc[df["res_gene"] == gene]
            # get the msa
            msa = generate_temp_total_fasta_file(all_recip_genes[gene], all_opposite_genes[opp_gene_name[0]])
            dm = calculator.get_distance(msa)
            # Create the tree
            upgmatree = constructor.upgma(dm)
            with open(os.path.join(os.getcwd(), "sorted_data", drug, "variations", "trees",
                                   gene.replace(".fasta", ".txt")), "w") as temp_tree:
                Phylo.write(upgmatree, temp_tree, 'newick')


phenotypes = ["res"]
drugs = ["CIPRO"]

for drug in drugs:
    for phenotype in phenotypes:
        get_trees(drug, phenotype)
