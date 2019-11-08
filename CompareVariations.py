import os
import statistics

import pandas as pd

from models import Blast as b
from models.BlastResult import BlastResult
from models.Gene import Gene


def get_organism_and_genes_by_phenotype(phenotype, drug):
    # Main Organism
    variations = pd.read_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_VariationsPercents.csv"))
    variations = variations.rename(columns={'Unnamed: 0': 'Genes'})
    genes = variations["Genes"]
    original_organism = variations.columns[1]
    variations = variations.drop(original_organism, axis=1)
    organism_db_path = os.path.join(os.getcwd(), "converted_data", original_organism, original_organism)

    print("Drug: ", drug)
    print("Phenotype: ", phenotype)
    print("Target Organism: ", original_organism)
    print("Number of res Genes: ", len(genes))

    return original_organism, variations, genes, organism_db_path


## This is used to find the genes that are unique to each phenotype after being recip blasted to the same phenotype
def compare_variations(drug, phenotype):
    op_phenotype = "res"
    if phenotype == "res":
        op_phenotype = "sus"

    # get data for all the genes of each phenotype
    target_organism, target_gene_data, target_genes, target_db_path = get_organism_and_genes_by_phenotype(phenotype, drug)
    opposite_organism, opposite_gene_data, opposite_genes, opposite_db_path = get_organism_and_genes_by_phenotype(op_phenotype, drug)

    headers = f"{phenotype}_organism,{phenotype}_gene,{op_phenotype}_organism,{op_phenotype}_gene,{phenotype}_length," \
              f"{op_phenotype}_length,{phenotype}_average,{op_phenotype}_average,{phenotype}_sd,{op_phenotype}_sd"
    final_df = pd.DataFrame(columns=headers.split(","))
    opposite_genes_list = list(opposite_genes)
    count = 0
    for suspect_gene in target_genes:
        target_gene = Gene(target_organism, str(suspect_gene.replace(".fasta", "")), get_info=False)
        blast_data = b.blast(target_gene.fasta_file, opposite_db_path)
        # if there is a match
        if len(blast_data) > 1:
            ## Now blast the returned gene against the original organism
            blast_result = BlastResult(blast_data)
            opposite_blast_gene = Gene(opposite_organism, blast_result.gene_name, get_info=False)
            opposite_gene_name = opposite_blast_gene.name + ".fasta"
            if opposite_gene_name in opposite_genes_list:
                print(f"Matches found: {count}")
                target_variation_values = target_gene_data[target_gene_data["Genes"] == suspect_gene]
                list_target_variation_values = list(list(target_variation_values.iloc[0]))
                list_target_variation_values.pop(0)
                converted_target_variation_values = [float(i) for i in list_target_variation_values]

                opposite_variation_values = opposite_gene_data[opposite_gene_data["Genes"] == opposite_gene_name]
                list_opposite_variation_values = list(list(opposite_variation_values.iloc[0]))
                list_opposite_variation_values.pop(0)
                converted_opposite_variation_values = [float(i) for i in list_opposite_variation_values]

                # get information for the differences in the variation of each set of genes
                target_average = statistics.mean(converted_target_variation_values)
                target_sd = statistics.stdev(converted_target_variation_values, target_average)
                opposite_average = statistics.mean(converted_opposite_variation_values)
                opposite_sd = statistics.stdev(converted_opposite_variation_values, opposite_average)

                # output it into the dataframe.
                final_df = final_df.append({f"{phenotype}_organism": target_organism, f"{phenotype}_gene": suspect_gene,
                                            f"{op_phenotype}_organism": opposite_organism, f"{op_phenotype}_gene": opposite_gene_name,
                                            f"{phenotype}_length": round(target_gene.length/3), f"{op_phenotype}_length": round(opposite_blast_gene.length/3),
                                            f"{phenotype}_average": round(target_average, 2), f"{op_phenotype}_average": round(opposite_average, 2),
                                            f"{phenotype}_sd": target_sd, f"{op_phenotype}_sd": opposite_sd}, ignore_index=True)

                count += 1

    final_df.to_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_VariationsCompared.csv"), index=False)


phenotypes = ["sus", "res"]
drugs = ["CIPRO"]

for drug in drugs:
    for phenotype in phenotypes:
        compare_variations(drug, phenotype)
