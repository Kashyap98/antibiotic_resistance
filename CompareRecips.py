import csv
import os

import pandas as pd

from models import Blast as b
from models.BlastResult import BlastResult
from models.Gene import Gene
from models.Mutation import Mutation


def write_cell(outputFile, data, isStart=False):
    if isStart:
        outputFile.write(data + ",")
    else:
        outputFile.write(str(data) + ",")


def gene_info_writer(target_organism, opposite_organism, target_gene, opposite_gene, bitscore, recip_bitscore, outputFile,
                     target_function, opposite_function, nonpolar_polar, nonpolar_acidic, nonpolar_basic,
                     polar_nonpolar, polar_acidic, polar_basic, acidic_nonpolar, acidic_polar, acidic_basic,
                     basic_nonpolar, basic_polar, basic_acidic, gap_reference, gap_matched):
    write_cell(outputFile, target_organism, isStart=True)
    write_cell(outputFile, opposite_organism)
    write_cell(outputFile, target_gene)
    write_cell(outputFile, opposite_gene)
    write_cell(outputFile, bitscore)
    write_cell(outputFile, recip_bitscore)
    write_cell(outputFile, target_function)
    write_cell(outputFile, opposite_function)
    write_cell(outputFile, nonpolar_polar)
    write_cell(outputFile, nonpolar_acidic)
    write_cell(outputFile, nonpolar_basic)
    write_cell(outputFile, polar_nonpolar)
    write_cell(outputFile, polar_acidic)
    write_cell(outputFile, polar_basic)
    write_cell(outputFile, acidic_nonpolar)
    write_cell(outputFile, acidic_polar)
    write_cell(outputFile, acidic_basic)
    write_cell(outputFile, basic_nonpolar)
    write_cell(outputFile, basic_polar)
    write_cell(outputFile, basic_acidic)
    write_cell(outputFile, gap_reference)
    write_cell(outputFile, gap_matched)
    outputFile.write("\n")
    
    
def get_organism_and_genes_by_phenotype(phenotype, drug):
    # Main Organism
    organism_file = pd.read_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + ".csv"), header=None)
    all_orgs = list(organism_file.iloc[:, 0])
    targetOrganism = str(all_orgs[0])
    organism_db_path = os.path.join(os.getcwd(), "converted_data", targetOrganism, targetOrganism)
    # Get file with all the Genes
    gene_file = pd.read_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_RecipGenes.csv"))
    all_genes = list(gene_file.iloc[:, 0])
    print("Drug: ", drug)
    print("Phenotype: ", phenotype)
    print("Target Organism: ", targetOrganism)
    print("Number of res Genes: ", len(all_genes))
    
    return targetOrganism, all_genes, organism_db_path


## This is used to find the genes that are unique to each phenotype after being recip blasted to the same phenotype
def get_genes_in_both_organisms(drug, phenotype):
    op_phenotype = "res"
    if phenotype == "res":
        op_phenotype = "sus"
    outputFile = open(os.path.join(os.getcwd(), "sorted_data", drug, "{}_{}_CommonGeneData.csv".format(phenotype, op_phenotype)), "w")
    outputFile.write("{}_organism, {}_organism, {}_gene, {}_gene, bitscore, recip_bitscore, {}_function, {}_function, "
                     "nonpolar_polar, nonpolar_acidic, nonpolar_basic, polar_nonpolar, polar_acidic, polar_basic, "
                     "acidic_nonpolar, acidic_polar, acidic_basic, basic_nonpolar, basic_polar, basic_acidic,"
                     "gap_reference, gap_matched \n".format(phenotype, op_phenotype, phenotype, op_phenotype, phenotype, op_phenotype))

    target_organism, target_genes, target_db_path = get_organism_and_genes_by_phenotype(phenotype, drug)
    opposite_organism, opposite_genes, opposite_db_path = get_organism_and_genes_by_phenotype(op_phenotype, drug)
    final_genes = []
    
    for suspect_gene in target_genes:
        target_gene = Gene(target_organism, str(suspect_gene.replace(".fasta", "")), get_info=True)
        blast_data = b.blast(target_gene.fasta_file, opposite_db_path)
        if len(blast_data) > 1:
            ## Now blast the returned gene against the original organism
            blast_result = BlastResult(blast_data)
            opposite_blast_gene = Gene(opposite_organism, blast_result.gene_name, get_info=True)
            if blast_result.bitscore >= 1000 \
                    and target_gene.lower_length <= blast_result.match_length <= target_gene.upper_length \
                    and blast_result.gene_name + ".fasta" in opposite_genes:
                # Blast new gene back to original to make sure it is reciprocal.
                recip_blast_data = b.blast(opposite_blast_gene.fasta_file, target_db_path)
                if len(recip_blast_data) > 1:
                    recip_blast_result = BlastResult(recip_blast_data)
                    if recip_blast_result.bitscore >= 1000 \
                            and opposite_blast_gene.lower_length <= recip_blast_result.match_length <= opposite_blast_gene.upper_length\
                            and recip_blast_result.gene_name + ".fasta" in target_genes:
                        if str(recip_blast_result.gene_name) != str(suspect_gene.replace(".fasta", "")):
                            ## Remove the gene from the matched list
                            target_genes.remove(suspect_gene)
                        else:
                            matchRatio = float(round(int(blast_result.match_length) / int(target_gene.length), 2))
                            final_genes.append(target_gene.name)
                            mutation_info = Mutation(target_gene.aa_sequence, opposite_blast_gene.aa_sequence)

                            gene_info_writer(target_organism, opposite_organism, target_gene.name, blast_result.gene_name,
                                             blast_result.bitscore, recip_blast_result.bitscore, outputFile,
                                             target_gene.function_data, opposite_blast_gene.function_data,
                                             mutation_info.nonpolar_polar, mutation_info.nonpolar_acidic, 
                                             mutation_info.nonpolar_basic, mutation_info.polar_nonpolar,
                                             mutation_info.polar_acidic, mutation_info.polar_basic,
                                             mutation_info.acidic_nonpolar, mutation_info.acidic_polar, 
                                             mutation_info.acidic_basic, mutation_info.basic_nonpolar, 
                                             mutation_info.basic_polar, mutation_info.basic_acidic, 
                                             mutation_info.gap_reference, mutation_info.gap_matched)

                    # Remove if anything does not match because we are looking for recips here.
                    else:
                        target_genes.remove(suspect_gene)
                else:
                    target_genes.remove(suspect_gene)
            else:
                target_genes.remove(suspect_gene)
        else:
            target_genes.remove(suspect_gene)
        print("LENGTH OF GENES: ", len(target_genes))

    outputFile.close()
    with open(os.path.join(os.getcwd(), "sorted_data", drug, "{}_{}_CommonGenes.csv".format(phenotype, op_phenotype)), 'w') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        wr.writerow(final_genes)


phenotypes = ["sus", "res"]
drugs = ["CIPRO"]

for drug in drugs:
    for phenotype in phenotypes:
        get_genes_in_both_organisms(drug, phenotype)
