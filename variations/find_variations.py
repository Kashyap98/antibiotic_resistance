import os

import pandas as pd

from models import Blast as b
from models.BlastResult import BlastResult
from models.Gene import Gene
from utils import dir_utils, gen_utils, output_util


def findstem(arr):
    n = len(arr)
    s = arr[0]
    l = len(s)
    res = ""
    for i in range(l):
        for j in range(i + 1, l + 1):
            stem = s[i:j]
            k = 1
            for k in range(1, n):
                if stem not in arr[k]:
                    break
            if k + 1 == n and len(res) < len(stem):
                res = stem
    return res


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


def get_lcs(original_gene_data, opposite_gene_data):
    res_genes = []
    sus_genes = []

    for gene, info in original_gene_data.items():
        organism = info[0]
        variation_score = info[1]
        new_gene = Gene(organism, gene.replace(".fasta", ""), get_info=True)
        res_genes.append(new_gene.aa_sequence)

    for gene, info in opposite_gene_data.items():
        organism = info[0]
        variation_score = info[1]
        opp_gene = Gene(organism, gene.replace(".fasta", ""), get_info=True)
        sus_genes.append(opp_gene.aa_sequence)

    res_lcs = findstem(res_genes)
    sus_lcs = findstem(sus_genes)
    combined_lcs = findstem([res_lcs, sus_lcs])

    return len(res_lcs), len(sus_lcs), len(combined_lcs)


def find_interesting(drug, phenotype):
    op_phenotype = gen_utils.get_op_phenotype(phenotype)
    drug_dirs = dir_utils.DrugDirs(drug, phenotype)
    all_organisms_list, recip_genes_df, recip_genes_data = gen_utils.get_organisms_genes_data_by_phenotype(drug,
                                                                                                           phenotype)
    recip_genes = list(recip_genes_df["0"])
    all_genes_length = len(recip_genes)

    opposite_organisms_list, opposite_genes_df, opposite_genes_data = gen_utils.get_organisms_genes_data_by_phenotype(
        drug, op_phenotype)
    opposite_genes = list(opposite_genes_df["0"])

    original_organism = all_organisms_list.pop(0)
    opposite_organism = opposite_organisms_list.pop(0)
    opposite_organism_dirs = dir_utils.OrganismDirs(opposite_organism)
    all_recip_genes = generate_recip_genes_mappings(recip_genes, original_organism, recip_genes_data)
    all_opposite_genes = generate_recip_genes_mappings(opposite_genes, opposite_organism, opposite_genes_data)

    headers = f"{phenotype}_organism,{phenotype}_gene,{op_phenotype}_organism,{op_phenotype}_gene,{phenotype}_length," \
              f"{op_phenotype}_length,{phenotype}_lcs,{op_phenotype}_lcs,combined_lcs,{phenotype}_function," \
              f"{op_phenotype}_function"
    final_df = pd.DataFrame(columns=headers.split(","))

    for key, value in all_recip_genes.items():
        original_gene_name = key
        original_gene = Gene(original_organism, original_gene_name.replace(".fasta", ""), get_info=True)
        print(f"{len(final_df)} / {all_genes_length}")

        blast_data = b.blast(original_gene.fasta_file, opposite_organism_dirs.database_dir, opposite_organism)
        if not blast_data:
            continue

        opposite_blast_gene = Gene(opposite_organism, blast_data.blast_gene.name, get_info=True)
        opposite_gene_name = opposite_blast_gene.name + ".fasta"
        if opposite_gene_name in opposite_genes:
            res_lcs, sus_lcs, combined_lcs = get_lcs(all_recip_genes[original_gene_name], all_opposite_genes[opposite_gene_name])
            final_df = final_df.append({f"{phenotype}_organism": original_organism,
                                        f"{phenotype}_gene": original_gene_name,
                                        f"{op_phenotype}_organism": opposite_organism,
                                        f"{op_phenotype}_gene": opposite_gene_name,
                                        f"{phenotype}_length": len(original_gene.aa_sequence),
                                        f"{op_phenotype}_length": len(opposite_blast_gene.aa_sequence),
                                        f"{phenotype}_lcs": round(res_lcs / len(original_gene.aa_sequence), 2),
                                        f"{op_phenotype}_lcs": round(sus_lcs / len(opposite_blast_gene.aa_sequence), 2),
                                        "combined_lcs": round(combined_lcs / len(original_gene.aa_sequence), 2),
                                        f"{phenotype}_function": original_gene.function_data,
                                        f"{op_phenotype}_function": opposite_blast_gene.function_data},
                                       ignore_index=True)

            opposite_genes.remove(opposite_gene_name)

    final_df.to_csv(os.path.join(drug_dirs.drug_dir, f"{phenotype}_Variations.csv"))


PHENOTYPES = ["res"]
DRUGS = ["CIPRO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        find_interesting(DRUG, PHENOTYPE)
