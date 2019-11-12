import csv
import os

import pandas as pd

from models import Blast as b
from models.Gene import Gene
from models.Mutation import Mutation
from utils import gen_utils, dir_utils, output_util


def get_organism_and_genes_by_phenotype(drug_dirs, phenotype):
    # Main Organism
    organism_file = pd.read_csv(os.path.join(drug_dirs.drug_dir, f"{phenotype}.csv"), header=None)
    all_orgs = list(organism_file.iloc[:, 0])
    target_organism = str(all_orgs[0])
    organism_dir = dir_utils.OrganismDirs(target_organism)
    # Get file with all the Genes
    gene_file = pd.read_csv(os.path.join(drug_dirs.drug_dir, phenotype + f"{phenotype}_RecipGenes.csv"))
    all_genes = list(gene_file.iloc[:, 0])

    return organism_dir.organism, all_genes, organism_dir.database_dir


# This is used to find the genes that are unique to each phenotype after being recip blasted to the same phenotype
def get_genes_in_both_organisms(drug, phenotype):
    op_phenotype = gen_utils.get_op_phenotype(phenotype)
    drug_dirs = dir_utils.DrugDirs(drug, phenotype)
    output_file = output_util.DataOutput(f"{phenotype}_{op_phenotype}_CommonGeneData.csv", drug, phenotype)
    output_file.add_headers(f"{phenotype}_organism, {op_phenotype}_organism, {phenotype}_gene, {op_phenotype}_gene, "
                            f"bitscore, recip_bitscore, {phenotype}_function, {op_phenotype}_function, nonpolar_polar, "
                            f"nonpolar_acidic, nonpolar_basic, polar_nonpolar, polar_acidic, polar_basic, "
                            f"acidic_nonpolar, acidic_polar, acidic_basic, basic_nonpolar, basic_polar, basic_acidic, "
                            f"gap_reference, gap_matched \n")

    target_organism, target_genes, target_db_path = get_organism_and_genes_by_phenotype(drug_dirs, phenotype)
    opposite_organism, opposite_genes, opposite_db_path = get_organism_and_genes_by_phenotype(drug_dirs, op_phenotype)
    final_genes = []

    for suspect_gene in target_genes:
        target_gene = Gene(target_organism, suspect_gene, get_info=True)
        blast_data = b.blast(target_gene.fasta_file, opposite_db_path, opposite_organism)

        if not blast_data:
            target_genes.remove(suspect_gene)
            continue

        # Now blast the returned gene against the original organism
        if not blast_data.bitscore >= 1000:
            target_genes.remove(suspect_gene)
            continue

        if not blast_data.target_gene.lower_length <= blast_data.match_length <= blast_data.target_gene.upper_length:
            target_genes.remove(suspect_gene)
            continue

        if not blast_data.gene_name + ".fasta" in opposite_genes:
            target_genes.remove(suspect_gene)
            continue

        # Blast new gene back to original to make sure it is reciprocal.
        recip_blast_data = b.blast(blast_data.blast_gene.fasta_file, target_db_path, blast_data.blast_gene.organism)

        if not recip_blast_data:
            target_genes.remove(suspect_gene)
            continue

        if not recip_blast_data.bitscore >= 1000:
            target_genes.remove(suspect_gene)
            continue

        if not recip_blast_data.target_gene.lower_length <= recip_blast_data.match_length <= \
               recip_blast_data.target_gene.upper_length:
            target_genes.remove(suspect_gene)
            continue

        if not recip_blast_data.gene_name + ".fasta" in opposite_genes:
            target_genes.remove(suspect_gene)
            continue

        if str(recip_blast_data.gene_name) != str(suspect_gene.replace(".fasta", "")):
            target_genes.remove(suspect_gene)
            continue

        final_genes.append(target_gene.name)
        mutation_info = Mutation(target_gene.aa_sequence, blast_data.blast_gene.aa_sequence)

        output_file.write_mutation_data(target_organism, opposite_organism, target_gene.name, blast_data.gene_name,
                                        blast_data.bitscore, recip_blast_data.bitscore, target_gene.function_data,
                                        blast_data.blast_gene.function_data, mutation_info.nonpolar_polar,
                                        mutation_info.nonpolar_acidic, mutation_info.nonpolar_basic,
                                        mutation_info.polar_nonpolar, mutation_info.polar_acidic,
                                        mutation_info.polar_basic, mutation_info.acidic_nonpolar,
                                        mutation_info.acidic_polar, mutation_info.acidic_basic,
                                        mutation_info.basic_nonpolar, mutation_info.basic_polar,
                                        mutation_info.basic_acidic, mutation_info.gap_reference,
                                        mutation_info.gap_matched)

        print("LENGTH OF GENES: ", len(target_genes))

    with open(os.path.join(os.getcwd(), "sorted_data", drug, f"{phenotype}_{op_phenotype}_CommonGenes.csv"), 'w') \
            as my_file:
        wr = csv.writer(my_file, quoting=csv.QUOTE_ALL)
        wr.writerow(final_genes)


PHENOTYPES = ["sus", "res"]
DRUGS = ["CIPRO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_genes_in_both_organisms(DRUG, PHENOTYPE)
