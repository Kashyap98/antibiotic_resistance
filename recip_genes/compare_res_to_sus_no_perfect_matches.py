import csv
import glob
import os
import threading
import time
from queue import Queue

import pandas as pd

from models import Blast as b
from models.Gene import Gene
from utils import gen_utils, dir_utils
from utils.output_util import DataOutput
from utils.dir_utils import OrganismDirs, DrugDirs


print_lock = threading.Lock()


class MyThread(threading.Thread):
    def __init__(self, input_queue, output_queue, unique_folder, op_orgs):
        super().__init__()
        self.queue = input_queue
        self.output_queue = output_queue
        self.unique_genes_folder = unique_folder
        self.daemon = True
        self.op_orgs = op_orgs
        self.final_gene_info = {}

    def run(self):
        while True:
            org = self.queue.get()
            if org is None:  # If you send `None`, the thread will exit.
                self.output_queue.put(self.final_gene_info)
                return
            self.do_work(org)

    def do_work(self, org):
        with print_lock:
            print(threading.currentThread().getName(), f"Org {org}")

        res_org_dirs = OrganismDirs(org, converted_ncbi_data=True)
        res_genes_copy = {}
        res_genes = os.listdir(res_org_dirs.gene_folder)
        count = 0
        with open(os.path.join(self.unique_genes_folder, f"{org}_unique.csv"), "w") as organism_unique_genes:
            organism_unique_genes.write(
                "gene,gene_description,hit_count,qcov_average,pident_average,perfect_matches,gene_length,function\n")

        for res_gene_name in res_genes:
            res_gene = Gene(org, res_gene_name)

            for op_organism in self.op_orgs:
                op_organism = str(op_organism)
                op_organism_dirs = OrganismDirs(op_organism, converted_ncbi_data=True)

                # First blast the first organism gene against the database of the gene in the list.
                blast_data = b.blast(res_gene, op_organism_dirs.database_dir, op_organism)

                # if there is a hit
                if not blast_data:
                    if res_gene_name not in res_genes_copy:
                        res_genes_copy[res_gene_name] = [res_gene.description, 0, 0, 0, 0, res_gene.length]
                    continue

                if res_gene_name in res_genes_copy:
                    old_data = res_genes_copy[res_gene_name]
                    new_count = old_data[1] + 1
                    new_qcov = (old_data[2] + blast_data.qcov) / 2
                    new_pident = (old_data[3] + blast_data.pident) / 2
                    if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                        new_perfect_matches = old_data[4] + 1
                    else:
                        new_perfect_matches = old_data[4]
                    res_genes_copy[res_gene_name] = [res_gene.description, new_count, new_qcov, new_pident,
                                                     new_perfect_matches, res_gene.length]
                else:
                    hit_count = 0
                    if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                        perfect_matches = 1
                    else:
                        perfect_matches = 0
                    res_genes_copy[res_gene_name] = [res_gene.description, hit_count, blast_data.qcov, blast_data.pident,
                                                     perfect_matches, res_gene.length]
            count += 1
            with print_lock:
                print(threading.currentThread().getName(), f"Org: {org} - Compare Gene: {count} / {len(res_genes)}")

        with open(os.path.join(self.unique_genes_folder, f"{org}_unique.csv"), "a") as organism_unique_genes:
            for gene, info in res_genes_copy.items():
                # if info[4] == 0:
                organism_unique_genes.write(f"{gene},{str(info[0])},{str(info[1])},{str(info[2])},{str(info[3])},"
                                            f"{str(info[4])},{str(info[5])}\n")

        with print_lock:
            print(threading.currentThread().getName(), f"Orgs left {self.queue.qsize()}")


def get_recips(drug, phenotype):
    # create output file and write the headers
    drug_dirs = DrugDirs(drug, phenotype)
    op_phenotype = gen_utils.get_op_phenotype(phenotype)
    drug_dirs.set_opposite_phenotype_file(op_phenotype)

    unique_genes_folder = os.path.join(drug_dirs.drug_dir, "unique_genes")
    dir_utils.generate_dir(unique_genes_folder)

    # Main Organism
    organism_file = pd.read_csv(drug_dirs.target_phenotype_file, header=None)
    all_orgs = list(organism_file.iloc[:, 0])
    op_organism_file = pd.read_csv(drug_dirs.op_phenotype_file, header=None)
    op_orgs = list(op_organism_file.iloc[:, 0])

    threads = []
    for t in range(10):
        q = Queue()
        output_queue = Queue()
        threads.append(MyThread(q, output_queue, unique_genes_folder, op_orgs))
        threads[t].start()
        time.sleep(0.1)

    thread_number = 0
    for x in all_orgs:
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
        t.queue.put(None)
        t.join()


PHENOTYPES = ["res"]
DRUGS = ["CIPRO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_recips(DRUG, PHENOTYPE)
