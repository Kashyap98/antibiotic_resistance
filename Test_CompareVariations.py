import os

import pandas as pd
from Bio import AlignIO
from Bio.Align import AlignInfo

from models import Blast as b
from models.BlastResult import BlastResult
from models.Gene import Gene


def lcs(X, Y, m, n):
    L = [[0 for x in range(n + 1)] for x in range(m + 1)]

    # Following steps build L[m+1][n+1] in bottom up fashion. Note
    # that L[i][j] contains length of LCS of X[0..i-1] and Y[0..j-1]
    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0 or j == 0:
                L[i][j] = 0
            elif X[i - 1] == Y[j - 1]:
                L[i][j] = L[i - 1][j - 1] + 1
            else:
                L[i][j] = max(L[i - 1][j], L[i][j - 1])

                # Following code is used to print LCS
    index = L[m][n]

    # Create a character array to store the lcs string
    lcs = [""] * (index + 1)
    lcs[index] = ""

    # Start from the right-most-bottom-most corner and
    # one by one store characters in lcs[]
    i = m
    j = n
    while i > 0 and j > 0:

        # If current character in X[] and Y are same, then
        # current character is part of LCS
        if X[i - 1] == Y[j - 1]:
            lcs[index - 1] = X[i - 1]
            i -= 1
            j -= 1
            index -= 1

        # If not same, then find the larger of two and
        # go in the direction of larger value
        elif L[i - 1][j] > L[i][j - 1]:
            i -= 1
        else:
            j -= 1

    return len(lcs)


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


def find_genes_same_phenotype(value, gene_file):
    for gene, info in value.items():
        organism = info[0]
        variation_score = info[1]
        new_gene = Gene(organism, gene.replace(".fasta", ""), get_info=True)
        gene_file.write("> " + new_gene.name + "\n")
        gene_file.write(new_gene.aa_sequence + "\n")


def find_lcs(original_gene, new_gene, final_df, final_percent_df, original_gene_name, organism, add_to_df=False):
    variation_score = lcs(original_gene.aa_sequence, new_gene.aa_sequence, len(original_gene.aa_sequence),
                          len(new_gene.aa_sequence))

    if add_to_df:
        final_df.set_value(original_gene_name, organism, round(variation_score))
        final_percent_df.set_value(original_gene_name, organism,
                                   round(variation_score / len(original_gene.aa_sequence), 2))


def find_lcs_same_phenotype_and_lsc_consensus(value, original_gene, original_gene_name, final_df, final_percent_df,
                                              opposite_organism, opposite_genes, opposite_blast_gene, final_consensus_df,
                                              op_phenotype, original_organism):
    orignal_genes_file = open(os.path.join(os.getcwd(), "sorted_data", drug, "original_genes.fasta"), "w")
    opposite_genes_file = open(os.path.join(os.getcwd(), "sorted_data", drug, "opposite_genes.fasta"), "w")
    opposite_gene_name = opposite_blast_gene.name + ".fasta"
    opposite_value = opposite_genes[opposite_gene_name]
    find_genes_same_phenotype(value, orignal_genes_file)
    find_genes_same_phenotype(opposite_value, opposite_genes_file)
    orignal_genes_file.close()
    opposite_genes_file.close()

    dir_copy = os.getcwd()
    os.chdir(os.path.join(os.path.join(os.getcwd(), "sorted_data", drug)))
    os.system("wsl.exe mafft --quiet original_genes.fasta > original_output.fasta")
    os.system("wsl.exe mafft --quiet opposite_genes.fasta > opposite_output.fasta")

    consensus_gene_records = AlignIO.read(open(os.path.join(os.getcwd(), "original_output.fasta")), format="fasta")
    opposite_consensus_gene_records = AlignIO.read(open(os.path.join(os.getcwd(), "opposite_output.fasta")), format="fasta")

    consensus_gene = AlignInfo.SummaryInfo(consensus_gene_records).gap_consensus()
    opposite_consensus_gene = AlignInfo.SummaryInfo(opposite_consensus_gene_records).gap_consensus()
    os.chdir(dir_copy)

    consensus_variation = lcs(consensus_gene, opposite_consensus_gene, len(consensus_gene), len(opposite_consensus_gene))
    final_consensus_df = final_consensus_df.append({f"{phenotype}_organism": original_organism,
                                                    f"{phenotype}_gene": original_gene_name,
                                                    f"{op_phenotype}_organism": opposite_organism,
                                                    f"{op_phenotype}_gene": opposite_gene_name,
                                                    f"{phenotype}_length": len(original_gene.aa_sequence),
                                                    f"{op_phenotype}_length": len(opposite_blast_gene.aa_sequence),
                                                    f"lcs": consensus_variation}, ignore_index=True)
    final_consensus_df.to_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_ConsensusVariation.csv"))

    return final_df, final_percent_df, final_consensus_df


def find_variations(drug, phenotype):
    op_phenotype = "res"
    if phenotype == "res":
        op_phenotype = "sus"

    all_organisms_list, recip_genes_df, recip_genes_data = get_organisms_genes_data_by_phenotype(drug, phenotype)
    recip_genes = list(recip_genes_df["0"])
    all_genes_length = len(recip_genes)

    opposite_organisms_list, opposite_genes_df, opposite_genes_data = get_organisms_genes_data_by_phenotype(drug,
                                                                                                         op_phenotype)
    opposite_genes = list(opposite_genes_df["0"])
    opposite_genes_length = len(opposite_genes)

    final_df = pd.DataFrame()
    final_percent_df = pd.DataFrame()
    original_organism = all_organisms_list.pop(0)
    opposite_organism = opposite_organisms_list.pop(0)
    opposite_organism_db_path = os.path.join(os.getcwd(), "converted_data", opposite_organism, opposite_organism)
    all_recip_genes = generate_recip_genes_mappings(recip_genes, original_organism, recip_genes_data)
    all_opposite_genes = generate_recip_genes_mappings(opposite_genes, opposite_organism, opposite_genes_data)

    headers = f"{phenotype}_organism,{phenotype}_gene,{op_phenotype}_organism,{op_phenotype}_gene,{phenotype}_length," \
              f"{op_phenotype}_length,lcs"
    final_consensus_df = pd.DataFrame(columns=headers.split(","))

    for key, value in all_recip_genes.items():
        original_gene_name = key
        original_gene = Gene(original_organism, original_gene_name.replace(".fasta", ""), get_info=True)
        print(f"{len(final_df)} / {all_genes_length}")

        blast_data = b.blast(original_gene.fasta_file, opposite_organism_db_path)
        if len(blast_data) > 1:
            blast_result = BlastResult(blast_data)
            opposite_blast_gene = Gene(opposite_organism, blast_result.gene_name, get_info=True)
            opposite_gene_name = opposite_blast_gene.name + ".fasta"
            if opposite_gene_name in opposite_genes:

                final_df, final_percent_df, final_consensus_df = \
                    find_lcs_same_phenotype_and_lsc_consensus(value, original_gene, original_gene_name, final_df,
                                                              final_percent_df, opposite_organism, all_opposite_genes,
                                                              opposite_blast_gene, final_consensus_df, op_phenotype,
                                                              original_organism)
                opposite_genes.remove(opposite_gene_name)

    final_df.to_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_Variations.csv"))
    final_percent_df.to_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_VariationsPercents.csv"))
    final_consensus_df.to_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_ConsensusVariation.csv"))


phenotypes = ["res"]
drugs = ["CIPRO"]

for drug in drugs:
    for phenotype in phenotypes:
        find_variations(drug, phenotype)