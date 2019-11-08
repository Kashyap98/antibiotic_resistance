import csv
import os

import pandas as pd

from models import Blast as b
from models.BlastResult import BlastResult
from models.Gene import Gene
from utils.output_util import DataOutput
from utils.dir_utils import OrganismDirs, DrugDirs


def gene_info_writer(target_organism, gene_name, bitscore, match_ratio, recip_gene_name, recip_organism_name,
                     outputFile):
    outputFile.write(target_organism + "," + gene_name + "," + recip_organism_name + "," + recip_gene_name + ","
                     + str(bitscore) + "," + str(match_ratio))
    outputFile.write("\n")


def get_recips(drug, phenotype):
    # Create output file
    output_file = DataOutput(f"{phenotype}_RecipBlastInfo.csv", drug, phenotype)
    output_file.add_headers("organism, feature_id, recip_organism, recip_feature_id, bitscore, match_ratio \n")
    drug_dirs = DrugDirs(drug, phenotype)

    # Get all the organisms for the phenotype
    gene_count = []
    organism = csv.reader(drug_dirs.target_phenotype_file, delimiter=',')
    all_organisms = []
    for row in organism:
        all_organisms.append(row[0])
    target_org = str(all_organisms.pop(0))
    match_count = len(all_organisms)
    final_gene_counter = {}
    target_org_dirs = OrganismDirs(target_org)

    all_genes_path = target_org_dirs.gene_folder
    organismDatabasePath = target_org_dirs.database_dir

    # get all the genes for the target organism
    all_genes = os.listdir(all_genes_path)
    print("All Organisms: ", all_organisms)
    print("Number of all Genes: ", len(all_genes))
    print("Originial Organism: ", target_org)
    genes_removed = {}
    count = 0
    for organism in all_organisms:
        ## For each other genome we want to compare to, Recip BLAST those genes and only keep the ones are the same
        organism_name = str(organism)
        print(organism_name)
        print(str(count) + "/" + str(len(all_organisms)))
        recip_organism_dirs = OrganismDirs(organism_name)
        database_path = recip_organism_dirs.database_dir

        for gene in all_genes:
            ## First blast the first organism gene against the database of the gene in the list.

            current_gene = Gene(target_org, gene)
            blastData = b.blast(current_gene, database_path, organism_name)

            if len(blastData) > 1:
                ## Now blast the returned gene against the original organism
                blastResult = BlastResult(blastData)
                recipGene = Gene(organism_name, blastResult.gene_name, get_info=False)
                if blastResult.bitscore >= 1000 \
                        and current_gene.lower_length <= blastResult.match_length <= current_gene.upper_length:
                    # Blast new gene back to original to make sure it is reciprocal.
                    recipBlastResultData = b.blast(recipGene.fasta_file, organismDatabasePath)
                    if len(recipBlastResultData) > 1:
                        recipBlastResult = BlastResult(recipBlastResultData)
                        if recipBlastResult.bitscore >= 1000 \
                                and recipGene.lower_length <= recipBlastResult.match_length <= recipGene.upper_length:
                            if str(recipBlastResult.gene_name) != str(gene.replace(".fasta", "")):
                                ## Remove the gene from the matched list
                                # all_genes.remove(gene)
                                genes_removed[gene] = "Recip Blast was not the same name!"
                            else:
                                matchRatio = float(round(int(blastResult.match_length) / int(current_gene.length), 2))
                                gene_info_writer(target_org, gene, blastResult.bitscore,
                                                 matchRatio, blastResult.gene_name, organism, output_file)
                                if gene in final_gene_counter:
                                    final_gene_counter[gene] += 1
                                else:
                                    final_gene_counter[gene] = 1
                        # Remove if anything does not match because we are looking for recips here.
                        else:
                            # all_genes.remove(gene)
                            print(gene)
                            genes_removed[gene] = "Recip Blast bitscore was < 1000!"
                    else:
                        # all_genes.remove(gene)
                        print(gene)
                        genes_removed[gene] = "Recip Blast did not match"
                else:
                    # all_genes.remove(gene)
                    print(gene)
                    genes_removed[gene] = "Match bitscore was < 1000!"
            else:
                # all_genes.remove(gene)
                print(gene)
                genes_removed[gene] = "There was no match!"

        print("Organism: " + organism_name)
        print("Matched Genes Length: ", len(all_genes))
        gene_count.append(len(all_genes))
        count += 1

    final_genes = []
    for key, value in final_gene_counter.items():
        if value == match_count:
            final_genes.append(key)
    total = pd.DataFrame(final_genes)
    total.to_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_RecipGenes.csv"), index=False)
    geneCountdf = pd.DataFrame(gene_count)
    geneCountdf.to_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_GeneCount.csv"), index=False)
    removed_genes_df = pd.DataFrame(genes_removed, index=[0])
    removed_genes_df.to_csv(os.path.join(os.getcwd(), "sorted_data", drug, f"{phenotype}_GenesNotMatched.csv"))


phenotypes = ["res"]
drugs = ["CIPRO"]

for drug in drugs:
    for phenotype in phenotypes:
        get_recips(drug, phenotype)
