import csv
import os
import threading
import time
from queue import Queue

import pandas as pd

from models import blast as b
from models.gene import Gene
from utils.dir_utils import OrganismDirs, DrugDirs
from utils import dir_utils, gen_utils

print_lock = threading.Lock()


class MyThread(threading.Thread):
    def __init__(self, input_queue, output_queue, sus_organisms, drug_dirs):
        super().__init__()
        self.queue = input_queue
        self.output_queue = output_queue
        self.daemon = True
        self.target_organism = ""
        self.sus_organisms = sus_organisms
        self.final_gene_info = {}
        self.drug_dirs = drug_dirs

    def run(self):
        while True:
            organism = self.queue.get()
            if organism is None:  # If you send `None`, the thread will exit.
                self.output_queue.put(self.final_gene_info)
                return
            self.do_work(organism)

    def do_work(self, organism):
        self.target_organism = organism
        self.final_gene_info = {}
        target_org_dirs = OrganismDirs(organism, converted_ncbi_data=True)
        all_genes_path = target_org_dirs.gene_folder
        completed_gene = ""

        # get all the genes for the target organism
        all_genes = os.listdir(all_genes_path)
        count = 0
        all_genes_count = len(all_genes)

        with open(os.path.join(self.drug_dirs.drug_dir, "complete_comparison", f"{organism}_complete_comparison.csv"),
                  "w") as recip_detailed:
            recip_detailed.write("gene,gene_description,sus_qcov_average,sus_pident_average,sus_match_count,"
                                 "sus_perfect_match_count,sus_matches,sus_perfect_matches\n")

        for gene in all_genes:
            count += 1
            completed_gene = gene
            # For each other genome we want to compare to, Recip BLAST those genes and only keep the ones are the same
            current_gene = Gene(self.target_organism, gene)

            for sus_organism in self.sus_organisms:
                sus_organism_dirs = OrganismDirs(str(sus_organism), converted_ncbi_data=True)
                sus_database_path = sus_organism_dirs.database_dir
                blast_data = b.blast(current_gene, sus_database_path, str(sus_organism))

                if not blast_data:
                    continue

                if f"{organism} - {gene}" in self.final_gene_info:
                    old_data = self.final_gene_info[f"{organism} - {gene}"]
                    gene_description = old_data[0]
                    new_qcov = (old_data[1] + blast_data.qcov) / 2
                    new_pident = (old_data[2] + blast_data.pident) / 2
                    if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                        new_count = old_data[3]
                        new_perfect_matches = old_data[4] + 1
                        new_match_organisms = old_data[5]
                        new_perfect_match_organisms = old_data[6]
                        new_perfect_match_organisms.append(sus_organism)
                    else:
                        new_count = old_data[3] + 1
                        new_perfect_matches = old_data[4]
                        new_perfect_match_organisms = old_data[6]
                        new_match_organisms = old_data[5]
                        new_match_organisms.append(sus_organism)

                    self.final_gene_info[f"{organism} - {gene}"] = [gene_description, new_qcov, new_pident, new_count,
                                                                    new_perfect_matches, new_match_organisms,
                                                                    new_perfect_match_organisms]
                else:
                    if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                        perfect_matches = 1
                        hit_count = 0
                        res_perfect_match_organisms = [sus_organism]
                        res_match_organisms = []
                    else:
                        hit_count = 1
                        perfect_matches = 0
                        res_perfect_match_organisms = []
                        res_match_organisms = [sus_organism]

                    self.final_gene_info[f"{organism} - {gene}"] = [current_gene.description, blast_data.qcov,
                                                                    blast_data.pident, hit_count, perfect_matches,
                                                                    res_match_organisms, res_perfect_match_organisms]

            with open(os.path.join(self.drug_dirs.drug_dir, "complete_comparison", f"{organism}_complete_comparison.csv"),
                      "a") as recip_detailed:
                data = self.final_gene_info[f"{organism} - {completed_gene}"]
                recip_detailed.write(f"{organism} - {completed_gene},{str(data[0])},{str(data[1])},{str(data[2])},{str(data[3])},"
                                     f"{str(data[4])},{str(data[5]).replace(',', '.')},"
                                     f"{str(data[6]).replace(',', '.')}\n")

            with print_lock:
                print(threading.currentThread().getName(), f"Genes Compared {count} / {all_genes_count}")

            # for gene, data in self.final_gene_info.items():
            #     recip_detailed.write(f"{gene},{str(data[0])},{str(data[1])},{str(data[2])},{str(data[3])},"
            #                          f"{str(data[4])},{str(data[5]).replace(',', '.')},"
            #                          f"{str(data[6]).replace(',', '.')}\n")


def get_recips(drug, phenotype):
    # Create output file
    drug_dirs = DrugDirs(drug, phenotype)
    dir_utils.generate_dir(os.path.join(drug_dirs.drug_dir, "complete_comparison"))
    op_phenotype = gen_utils.get_op_phenotype(phenotype)
    drug_dirs.set_opposite_phenotype_file(op_phenotype)

    # Get all the organisms for the phenotype
    all_organisms_raw = list(csv.reader(open(drug_dirs.target_phenotype_file), delimiter=','))
    res_organisms = []
    for org in all_organisms_raw:
        res_organisms.append(org[0])

    op_organism_file = pd.read_csv(drug_dirs.op_phenotype_file, header=None)
    op_orgs = list(op_organism_file.iloc[:, 0])

    # get all the genes for the target organism

    print(f"{phenotype} Organisms: ", res_organisms)
    print(f"{op_phenotype} Organisms: ", op_orgs)

    threads = []
    for t in range(6):
        q = Queue()
        output_queue = Queue()
        threads.append(MyThread(q, output_queue, op_orgs, drug_dirs))
        threads[t].start()
        time.sleep(0.1)

    thread_number = 0
    for x in res_organisms:
        threads[thread_number].queue.put(x)

        thread_number += 1
        if thread_number == 6:
            thread_number = 0

    for t in threads:
        t.queue.put(None)

    for t in threads:
        while t.queue.qsize() > 0:
            pass

    for t in threads:
        t.queue.put(None)
        t.join()


PHENOTYPES = ["res"]
DRUGS = ["AMOXO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_recips(DRUG, PHENOTYPE)
