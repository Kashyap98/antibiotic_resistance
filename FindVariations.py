import os

import pandas as pd
from Bio import AlignIO
from Bio.Align import AlignInfo

from models import Blast as b
from models.BlastResult import BlastResult
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


def find_lcs_same_phenotype(value, original_gene, original_gene_name, final_df, final_percent_df, gene_file, add_to_df=False):
    for gene, info in value.items():
        organism = info[0]
        variation_score = info[1]
        new_gene = Gene(organism, gene.replace(".fasta", ""), get_info=True)
        gene_file.write("> " + new_gene.name + "\n")
        gene_file.write(new_gene.aa_sequence + "\n")
        variation_score = lcs(original_gene.aa_sequence, new_gene.aa_sequence, len(original_gene.aa_sequence),
                              len(new_gene.aa_sequence))

        if add_to_df:
            final_df.set_value(original_gene_name, organism, round(variation_score))
            final_percent_df.set_value(original_gene_name, organism,
                                       round(variation_score / len(original_gene.aa_sequence), 2))


def generate_consensus(records):
    consensus = ""
    for i in range(0, len(str(records[0].seq))):
        aa_dict = {}
        for record in records:
            aa = str(record.seq[i])
            if aa in aa_dict:
                aa_dict[aa] += 1
            else:
                aa_dict[aa] = 1
        consensus += max(aa_dict, key=aa_dict.get)
    return consensus


def find_lcs_same_phenotype_and_lsc_consensus(value, original_gene, original_gene_name, final_df, final_percent_df,
                                              opposite_organism, opposite_genes, opposite_blast_gene, final_consensus_df,
                                              op_phenotype, original_organism):
    orignal_genes_file = open(os.path.join(os.getcwd(), "sorted_data", drug, "original_genes.fasta"), "w")
    opposite_genes_file = open(os.path.join(os.getcwd(), "sorted_data", drug, "opposite_genes.fasta"), "w")

    opposite_gene_name = opposite_blast_gene.name + ".fasta"
    opposite_value = opposite_genes[opposite_gene_name]
    find_lcs_same_phenotype(value, original_gene, original_gene_name, final_df,
                                                         final_percent_df, orignal_genes_file, add_to_df=True)
    find_lcs_same_phenotype(opposite_value, opposite_blast_gene, opposite_gene_name,
                                                         final_df, final_percent_df, opposite_genes_file)

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
