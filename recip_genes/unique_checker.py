import os

import pandas as pd

from models import Blast as b
from models.Gene import Gene
from utils.output_util import DataOutput
from utils.dir_utils import OrganismDirs, DrugDirs
import utils.gen_utils as gen_utils


# This is used to find the genes that are unique to each phenotype after being recip blasted to the same phenotype
def get_uniques(drug, phenotype):
    # create output file and write the headers
    output_file = DataOutput(f"{phenotype}_UniqueBlastInfo.csv", drug, phenotype)
    output_file.add_headers("organism, feature_id, recip_organism, recip_feature_id, bitscore, match_ratio,"
                            "contig_id, type_data, location, start, stop, strand, function_data,"
                            "aliases, figfam, evidence_codes \n")
    drug_dirs = DrugDirs(drug, phenotype)
    op_phenotype = gen_utils.get_op_phenotype(phenotype)
    drug_dirs.set_opposite_phenotype_file(op_phenotype)

    # Main Organism
    organism_file = pd.read_csv(drug_dirs.target_phenotype_file, header=None)
    all_orgs = list(organism_file.iloc[:, 0])
    target_organism_dirs = OrganismDirs(all_orgs[0])

    # Get file with all the Genes
    gene_file = pd.read_csv(drug_dirs.res_recip_genes_file, header=None)
    all_genes = list(gene_file.iloc[:, 0])

    print("Number of all Genes: ", len(all_genes))

    unique_matched = {}
    op_organism_file = pd.read_csv(drug_dirs.op_phenotype_file, header=None)
    op_orgs = list(op_organism_file.iloc[:, 0])
    count = 0
    for op_organism in op_orgs:
        # For each other genome we want to compare to, Recip BLAST those genes and only keep the ones are the same
        op_organism = str(op_organism)
        print(str(count) + "/" + str(len(op_orgs)))
        op_organism_dirs = OrganismDirs(op_organism)

        for gene in all_genes:
            # First blast the first organism gene against the database of the gene in the list.
            target_gene = Gene(target_organism_dirs.organism, gene, get_info=True)
            blast_data = b.blast(target_gene, op_organism_dirs.database_dir, op_organism)
            # print(blast_data)
            matched_organisms = []
            matched_genes = []

            # If we have already found it in unique_matched
            if gene in unique_matched:
                matched_organisms = unique_matched[gene][2]
                matched_genes = unique_matched[gene][3]

            # if there is a hit
            if not blast_data:
                all_genes.remove(gene)
                continue

            output_file.write_unique_info(target_gene, blast_data.blast_gene, blast_data.bitscore,
                                          blast_data.match_ratio)

            # if the result is unique enough
            if not blast_data.bitscore < 1600:
                all_genes.remove(gene)
                continue

            matched_genes.append(blast_data.gene_name)
            matched_organisms.append(op_organism)
            unique_matched[gene] = (target_organism_dirs.organism, len(matched_organisms), matched_organisms,
                                    matched_genes)

        print("Organism: " + op_organism)
        print("Matched Genes Length: ", len(all_genes))
        count += 1

    # output the results as a csv, for some reason I started doing it this way for a while
    with open(os.path.join(drug_dirs.drug_dir, f"{phenotype}_UniqueMatches.csv"), "w") as unique_matches:
        for key in unique_matched.keys():
            row = unique_matched[key]
            unique_matches.write(f"{key},{row[0]},{row[1]},{row[2]},{row[3]}\n")

    with open(drug_dirs.unique_genes_file, "w") as unique_genes:
        for gene in all_genes:
            unique_genes.write(f"{gene}\n")


PHENOTYPES = ["res"]
DRUGS = ["AMOXO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_uniques(DRUG, PHENOTYPE)
