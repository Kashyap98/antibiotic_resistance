import os

import pandas as pd
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import utils.gen_utils as gen_utils
import utils.dir_utils as dir_utils

constructor = DistanceTreeConstructor()
calculator = DistanceCalculator('identity')


def get_trees(drug, phenotype):
    op_phenotype = gen_utils.get_op_phenotype(phenotype)
    dir_utils.generate_dir(os.path.join(os.getcwd(), "..", "sorted_data", drug, "variations", "alignments"))
    drug_dirs = dir_utils.DrugDirs(drug, phenotype)

    # get all the organism data.
    df = pd.read_csv(os.path.join(drug_dirs.variations_dir, f"{phenotype}_both_var.csv"), index_col=0)
    all_organisms_list, recip_genes_df, recip_genes_data = gen_utils.get_organisms_genes_data_by_phenotype(drug,
                                                                                                           phenotype)
    recip_genes = list(recip_genes_df["0"])
    original_organism = all_organisms_list.pop(0)

    opposite_organisms_list, opposite_genes_df, opposite_genes_data = gen_utils.get_organisms_genes_data_by_phenotype(
        drug, op_phenotype)
    opposite_genes = list(opposite_genes_df["0"])
    opposite_organism = opposite_organisms_list.pop(0)
    all_recip_genes = gen_utils.generate_recip_genes_mappings(recip_genes, original_organism, recip_genes_data)
    all_opposite_genes = gen_utils.generate_recip_genes_mappings(opposite_genes, opposite_organism, opposite_genes_data)

    for gene in df["res_gene"]:
        if gene in all_recip_genes:
            opp_gene_name = df["sus_gene"].loc[df["res_gene"] == gene]
            # get the msa
            msa = gen_utils.generate_temp_total_fasta_file(drug, all_recip_genes[gene],
                                                           all_opposite_genes[opp_gene_name[0]], gene)
            dm = calculator.get_distance(msa)
            # Create the tree
            upgmatree = constructor.upgma(dm)
            with open(os.path.join(drug_dirs.variations_dir, "trees", gene.replace(".fasta", ".txt")), "w") \
                    as temp_tree:
                Phylo.write(upgmatree, temp_tree, 'newick')


phenotypes = ["res"]
drugs = ["CIPRO"]

for drug in drugs:
    for phenotype in phenotypes:
        get_trees(drug, phenotype)
