import csv
import glob
import os
import threading
import time
from queue import Queue

import pandas as pd

from models import blast as b
from models.gene import Gene
from utils.output_util import DataOutput
from utils.dir_utils import OrganismDirs, DrugDirs


print_lock = threading.Lock()


class MyThread(threading.Thread):
    def __init__(self, input_queue, output_queue, target_organism, organisms):
        super().__init__()
        self.queue = input_queue
        self.output_queue = output_queue
        self.daemon = True
        self.target_organism = target_organism
        self.organisms = organisms
        self.final_gene_info = {}

    def run(self):
        while True:
            gene = self.queue.get()
            if gene is None:  # If you send `None`, the thread will exit.
                self.output_queue.put(self.final_gene_info)
                return
            self.do_work(gene)

    def do_work(self, gene):
        for organism in self.organisms:
            # For each other genome we want to compare to, Recip BLAST those genes and only keep the ones are the same
            organism_name = str(organism)
            recip_organism_dirs = OrganismDirs(organism_name, converted_ncbi_data=True)
            database_path = recip_organism_dirs.database_dir

            # First blast the first organism gene against the database of the gene in the list.
            current_gene = Gene(self.target_organism, gene)
            blast_data = b.blast(current_gene, database_path, organism_name)

            if not blast_data:
                continue

            if gene in self.final_gene_info:
                old_data = self.final_gene_info[gene]
                gene_description = old_data[0]
                new_count = old_data[1] + 1
                new_qcov = (old_data[2] + blast_data.qcov) / 2
                new_pident = (old_data[3] + blast_data.pident) / 2
                if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                    new_perfect_matches = old_data[4] + 1
                else:
                    new_perfect_matches = old_data[4]

                self.final_gene_info[gene] = [gene_description, new_count, new_qcov, new_pident, new_perfect_matches]
            else:
                hit_count = 1
                if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                    perfect_matches = 1
                else:
                    perfect_matches = 0
                self.final_gene_info[gene] = [current_gene.description, hit_count, blast_data.qcov, blast_data.pident,
                                              perfect_matches]

        with print_lock:
            print(threading.currentThread().getName(), f"Genes left {self.queue.qsize()}")


def get_recips(drug, phenotype):
    # Create output file
    drug_dirs = DrugDirs(drug, phenotype)
    with open(os.path.join(drug_dirs.drug_dir, f"{phenotype}_aa_recip_detailed_info.csv"), "w") as recip_detailed:
        recip_detailed.write("gene,gene_description,hit_count,qcov_average,pident_average,perfect_matches\n")
    # Get all the organisms for the phenotype
    all_organisms_raw = list(csv.reader(open(drug_dirs.target_phenotype_file), delimiter=','))
    all_organisms = []
    for org in all_organisms_raw:
        all_organisms.append(org[0])
    target_org = str(all_organisms.pop(0))
    target_org_dirs = OrganismDirs(target_org, converted_ncbi_data=True)

    all_genes_path = target_org_dirs.gene_folder

    # get all the genes for the target organism
    all_genes = os.listdir(all_genes_path)

    print("All Organisms: ", all_organisms)
    print("Number of all Genes: ", len(all_genes))
    print("Original Organism: ", target_org)

    threads = []
    for t in range(10):
        q = Queue()
        output_queue = Queue()
        threads.append(MyThread(q, output_queue, target_org, all_organisms))
        threads[t].start()
        time.sleep(0.1)

    thread_number = 0
    for x in all_genes:
        threads[thread_number].queue.put(x)

        thread_number += 1
        if thread_number == 10:
            thread_number = 0

    for t in threads:
        t.queue.put(None)

    for t in threads:
        while t.queue.qsize() > 0:
            pass

    for t in threads:
        while not t.output_queue.empty():
            thread_data = t.output_queue.get()
            with open(os.path.join(drug_dirs.drug_dir, f"{phenotype}_aa_recip_detailed_info.csv"), "a") as recip_detailed:
                for gene, data in thread_data.items():
                    recip_detailed.write(f"{gene},{str(data[0])},{str(data[1])},{str(data[2])},{str(data[3])},"
                                         f"{str(data[4])}\n")
            print(data)

        t.queue.put(None)
        t.join()


PHENOTYPES = ["res"]
DRUGS = ["AMOXO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_recips(DRUG, PHENOTYPE)
