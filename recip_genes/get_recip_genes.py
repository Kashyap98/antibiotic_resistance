import csv
import os

import pandas as pd

from models import Blast as b
from models.Gene import Gene
from utils.output_util import DataOutput
from utils.dir_utils import OrganismDirs, DrugDirs


def get_recips(drug, phenotype):
    # Create output file
    output_file = DataOutput(f"{phenotype}_RecipBlastInfo.csv", drug, phenotype)
    output_file.add_headers("organism, feature_id, recip_organism, recip_feature_id, bitscore, match_ratio \n")
    drug_dirs = DrugDirs(drug, phenotype)

    # Get all the organisms for the phenotype
    gene_count = []
    all_organisms_raw = list(csv.reader(open(drug_dirs.target_phenotype_file), delimiter=','))
    all_organisms = []
    for org in all_organisms_raw:
        all_organisms.append(org[0])
    target_org = str(all_organisms.pop(0))
    match_count = len(all_organisms)
    final_gene_counter = {}
    target_org_dirs = OrganismDirs(target_org)

    all_genes_path = target_org_dirs.gene_folder
    organism_database_path = target_org_dirs.database_dir

    # get all the genes for the target organism
    all_genes = os.listdir(all_genes_path)
    print("All Organisms: ", all_organisms)
    print("Number of all Genes: ", len(all_genes))
    print("Original Organism: ", target_org)
    genes_removed = {}
    count = 1
    for organism in all_organisms:
        # For each other genome we want to compare to, Recip BLAST those genes and only keep the ones are the same
        organism_name = str(organism)
        print(organism_name)
        print(str(count) + "/" + str(len(all_organisms) + 1))
        recip_organism_dirs = OrganismDirs(organism_name)
        database_path = recip_organism_dirs.database_dir

        for gene in all_genes:
            # First blast the first organism gene against the database of the gene in the list.
            current_gene = Gene(target_org, gene)
            blast_data = b.blast(current_gene, database_path, organism_name)

            if not blast_data:
                # all_genes.remove(gene)
                genes_removed[gene] = "There was no match!"
                continue

            if not blast_data.blast_in_threshold:
                # all_genes.remove(gene)
                genes_removed[gene] = "Match bitscore was < 1000 or Lengths did not match!"
                continue

            recip_blast = b.blast(blast_data.blast_gene, organism_database_path, target_org)

            if not recip_blast:
                # all_genes.remove(gene)
                genes_removed[gene] = "Recip Blast did not match"
                continue

            if not recip_blast.blast_in_threshold:
                # all_genes.remove(gene)
                genes_removed[gene] = "Recip Blast bitscore was < 1000!"
                continue

            if current_gene.name == recip_blast.gene_name:
                output_file.write_recip_info(target_org, gene, blast_data.bitscore, blast_data.match_ratio,
                                             blast_data.gene_name, organism)
                if gene in final_gene_counter:
                    final_gene_counter[gene] += 1
                else:
                    final_gene_counter[gene] = 1
            else:
                # all_genes.remove(gene)
                genes_removed[gene] = "Recip Blast was not the same name!"
                continue

        print("Organism: " + organism_name)
        print("Matched Genes Length: ", len(all_genes))
        gene_count.append(len(all_genes))
        count += 1

    final_genes = []
    for key, value in final_gene_counter.items():
        if value == match_count:
            final_genes.append(key)
    # total = pd.DataFrame(final_genes)
    # total.to_csv(os.path.join(drug_dirs.drug_dir, f"{phenotype}_RecipGenes.csv"), index=False)
    with open(os.path.join(drug_dirs.drug_dir, f"{phenotype}_RecipGenes.csv"), "w") as total_recip_genes:
        for f_gene in final_genes:
            total_recip_genes.write(f"{f_gene}\n")

    with open(os.path.join(drug_dirs.drug_dir, f"{phenotype}_GeneCount.csv"), "w") as gene_count_file:
        for g_count in gene_count:
            gene_count_file.write(f"{g_count}\n")

    with open(os.path.join(drug_dirs.drug_dir, f"{phenotype}_GenesNotMatched.csv"), "w") as not_matched:
        for key, value in genes_removed.items():
            not_matched.write(f"{key},{value}\n")


PHENOTYPES = ["res"]
DRUGS = ["AMOXO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_recips(DRUG, PHENOTYPE)
