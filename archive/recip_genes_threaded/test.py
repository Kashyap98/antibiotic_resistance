import csv
import os
import threading
import time
from queue import Queue


from models import blast as b
from models.gene import Gene
from utils import gen_utils, output_util, dir_utils


def get_recips(drug, phenotype):
    # Create output file
    drug_dirs = dir_utils.DrugDirs(drug, phenotype)

    # Get all the organisms for the phenotype
    res_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.res_file)
    target_org = str(res_organisms.pop(0))
    target_org_dirs = dir_utils.OrganismDirs(target_org, converted_ncbi_data=True)

    # get all the genes for the target organism
    all_genes = os.listdir(target_org_dirs.gene_folder)
    final_gene_info = {}
    print("All Organisms: ", res_organisms)
    print("Number of all Genes: ", len(all_genes))
    print("Original Organism: ", target_org)

    for organism in res_organisms:
        # For each other genome we want to compare to, Recip BLAST those genes and only keep the ones are the same
        organism_name = str(organism)
        recip_organism_dirs = dir_utils.OrganismDirs(organism_name, converted_ncbi_data=True)
        database_path = recip_organism_dirs.database_dir

        for gene in all_genes:
            # First blast the first organism gene against the database of the gene in the list.
            current_gene = Gene(target_org, gene)
            blast_data = b.blast(current_gene, database_path, organism_name)

            if not blast_data:
                break

            res_blast_data = b.blast(blast_data.blast_gene,
                                     dir_utils.OrganismDirs(blast_data.blast_gene.organism,
                                                            converted_ncbi_data=True).database_dir,
                                     blast_data.blast_gene.organism)

            if not res_blast_data:
                break

            # if not blast_data.above_bitscore_threshold_check():
            #     # with print_lock:
            #     #     print(threading.currentThread().getName(), f"Thrown Out bitscore {current_gene.description}")
            #     break

            # if not blast_data.within_size_constraints_check():
            #     # with print_lock:
            #     #     print(threading.currentThread().getName(), f"Thrown Out size {current_gene.description}")
            #     break

            if blast_data.qcov > 100:
                print(f"QCOV > 100 {current_gene.description}")

            if gene in final_gene_info:
                old_data = final_gene_info[gene]
                gene_description = old_data[0]
                new_count = old_data[1] + 1
                new_bitscore = (old_data[2] + blast_data.bitscore) / 2
                new_qcov = (old_data[3] + blast_data.qcov) / 2
                new_pident = (old_data[4] + blast_data.pident) / 2
                new_genes = old_data[6]
                new_genes.append(f"{organism}~{blast_data.blast_gene.name}")
                if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                    new_perfect_matches = old_data[5] + 1
                else:
                    new_perfect_matches = old_data[5]

                final_gene_info[gene] = [gene_description, new_count, new_bitscore, new_qcov, new_pident,
                                              new_perfect_matches, new_genes]
            else:
                hit_count = 1
                if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                    perfect_matches = 1
                else:
                    perfect_matches = 0
                final_gene_info[gene] = [current_gene.description, hit_count, blast_data.qcov, blast_data.bitscore,
                                              blast_data.pident, perfect_matches,
                                              [f"{target_org}~{gene}",
                                               f"{organism}~{blast_data.blast_gene.name}"]]


PHENOTYPES = ["res"]
DRUGS = ["CIPRO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_recips(DRUG, PHENOTYPE)
