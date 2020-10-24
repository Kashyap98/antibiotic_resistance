import csv
import glob
import os
import threading
import time
from queue import Queue

import pandas as pd

from models import blast as b
from models.gene import Gene
from utils import gen_utils, output_util
from utils.dir_utils import OrganismDirs, DrugDirs


print_lock = threading.Lock()


class MyThread(threading.Thread):
    def __init__(self, input_queue, output_queue, target_organism, organisms, op_orgs):
        super().__init__()
        self.queue = input_queue
        self.output_queue = output_queue
        self.daemon = True
        self.target_organism = target_organism
        self.organisms = organisms
        self.op_orgs = op_orgs
        self.final_gene_info = {}

    def run(self):
        while True:
            gene_list = self.queue.get()
            if gene_list is None:  # If you send `None`, the thread will exit.
                self.output_queue.put(self.final_gene_info)
                return
            self.do_work(gene_list)

    def do_work(self, gene_list):

        cleaned_gene_list = gene_list.replace("[", "").replace("]", "").replace("'", "").replace(" ", "")
        final_gene_list = cleaned_gene_list.split("+")

        for gene_info in final_gene_list:

            res_organism, res_gene_name = gene_info.split("~")
            current_gene = Gene(res_organism, res_gene_name, get_info=False)

            for op_org in self.op_orgs:

                sus_org_dirs = OrganismDirs(op_org, converted_ncbi_data=True)
                sus_database_path = sus_org_dirs.database_dir

                blast_data = b.blast(current_gene, sus_database_path, op_org)

                if not blast_data:
                    continue

                if gene_info in self.final_gene_info:
                    old_data = self.final_gene_info[gene_info]
                    gene_description = old_data[0]
                    new_qcov = (old_data[1] + blast_data.qcov) / 2
                    new_pident = (old_data[2] + blast_data.pident) / 2

                    if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                        new_count = old_data[3]
                        new_perfect_matches = old_data[4] + 1
                        new_match_organisms = old_data[5]
                        new_perfect_match_organisms = old_data[6]
                        new_perfect_match_organisms.append(op_org)
                    else:
                        new_count = old_data[3] + 1
                        new_perfect_matches = old_data[4]
                        new_perfect_match_organisms = old_data[6]
                        new_match_organisms = old_data[5]
                        new_match_organisms.append(op_org)

                    self.final_gene_info[gene_info] = [gene_description, new_qcov, new_pident, new_count,
                                                  new_perfect_matches,
                                                  new_match_organisms, new_perfect_match_organisms]
                else:
                    if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                        perfect_matches = 1
                        hit_count = 0
                        res_match_organisms = []
                        res_perfect_match_organisms = [op_org]
                        
                    else:
                        perfect_matches = 0
                        hit_count = 1
                        res_perfect_match_organisms = []
                        res_match_organisms = [op_org]

                    self.final_gene_info[gene_info] = [current_gene.description, blast_data.qcov,
                                                  blast_data.pident, hit_count, perfect_matches, res_match_organisms,
                                                  res_perfect_match_organisms]

        with print_lock:
            print(threading.currentThread().getName(), f"Genes left {self.queue.qsize()}")


def get_recips(drug, phenotype):
    # Create output file
    drug_dirs = DrugDirs(drug, phenotype)
    output_file = output_util.OutputFile(drug_dirs.res_recip_to_all_sus, header_list=["gene", "gene_description",
                                                                                      "qcov_average", "pident_average",
                                                                                      "hit_count",
                                                                                      "perfect_match_count",
                                                                                      "hit_orgs", "perfect_match_orgs"])

    res_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.res_file)
    sus_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.sus_file)
    target_org = str(res_organisms.pop(0))

    all_genes_df = pd.read_csv(drug_dirs.res_recip_genes_file, header=0)
    filtered_df = all_genes_df[all_genes_df["hit_count"] == len(res_organisms)]

    # get all the genes for the target organism
    all_genes = list(filtered_df["recip_genes"])

    print("All Organisms: ", sus_organisms)
    print("Number of all Genes: ", len(all_genes))
    print("Original Organism: ", target_org)

    threads = []
    for t in range(10):
        q = Queue()
        output_queue = Queue()
        threads.append(MyThread(q, output_queue, target_org, res_organisms, sus_organisms))
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
            output_file.write_data_dict_to_output_file(thread_data)

        t.queue.put(None)
        t.join()


PHENOTYPES = ["res"]
DRUGS = ["CIPRO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_recips(DRUG, PHENOTYPE)
