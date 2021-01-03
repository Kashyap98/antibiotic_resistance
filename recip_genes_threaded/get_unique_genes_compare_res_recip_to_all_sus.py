import os
import threading
import time
from queue import Queue

import pandas as pd

from models import blast as b
from models.gene import Gene
from utils import gen_utils, output_util, dir_utils
from utils.dir_utils import OrganismDirs, DrugDirs


print_lock = threading.Lock()


class MyThread(threading.Thread):
    def __init__(self, input_queue, output_queue, op_orgs, drug_dirs):
        super().__init__()
        self.queue = input_queue
        self.output_queue = output_queue
        self.daemon = True
        self.op_orgs = op_orgs
        self.drug_dirs = drug_dirs
        self.final_gene_info = {}

    def run(self):
        while True:
            gene_list = self.queue.get()
            if gene_list is None:  # If you send `None`, the thread will exit.
                self.output_queue.put(self.final_gene_info)
                return
            self.do_work(gene_list)

    def do_work(self, gene_list):
        count = 0
        gene_list_length = len(gene_list)
        organism = gene_list[0].organism
        output_file = output_util.OutputFile(os.path.join(self.drug_dirs.comparison_dir, f"{organism}_genes.csv"),
                                             header_list=["gene", "gene_description", "qcov_average", "pident_average",
                                                          "hit_count", "perfect_match_count", "hit_orgs",
                                                          "perfect_match_orgs"])

        for current_gene in gene_list:
            count += 1
            # self.final_gene_info[f"{current_gene.name}-{count}"] = [current_gene.description, 0, 0, 0, 0, [], []]
            for op_org in self.op_orgs:

                sus_org_dirs = OrganismDirs(op_org, converted_ncbi_data=True)
                sus_database_path = sus_org_dirs.database_dir

                blast_data = b.blast(current_gene, sus_database_path, op_org)

                if not blast_data:
                    if f"{current_gene.name}-{count}" not in self.final_gene_info:
                        self.final_gene_info[f"{current_gene.name}-{count}"] = [current_gene.description, 0, 0, 0, 0, [], []]
                    continue

                if f"{current_gene.name}-{count}" in self.final_gene_info:
                    old_data = self.final_gene_info[f"{current_gene.name}-{count}"]
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

                    self.final_gene_info[f"{current_gene.name}-{count}"] = [gene_description, new_qcov, new_pident, new_count,
                                                               new_perfect_matches, new_match_organisms,
                                                               new_perfect_match_organisms]
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

                    self.final_gene_info[f"{current_gene.name}-{count}"] = [current_gene.description, blast_data.qcov,
                                                               blast_data.pident, hit_count, perfect_matches,
                                                               res_match_organisms, res_perfect_match_organisms]

            with print_lock:
                print(threading.currentThread().getName(), f"Genes Completed {count} / {gene_list_length}")

        output_file.write_data_dict_to_output_file(self.final_gene_info)


def get_recips(drug, phenotype):
    # Create output file
    drug_dirs = DrugDirs(drug, phenotype)
    dir_utils.generate_dir(drug_dirs.comparison_dir)

    res_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.res_file)
    sus_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.sus_file)

    all_genes_df = pd.read_csv(drug_dirs.res_recip_genes_file, header=0)
    filtered_df = all_genes_df[all_genes_df["hit_count"] == len(res_organisms) - 1]

    # get all the genes for the target organism
    all_genes = list(filtered_df["recip_genes"])

    genes_by_org = {}
    for org in res_organisms:
        genes_by_org[org] = list()

    for gene_list in all_genes:

        cleaned_gene_list = gene_list.replace("[", "").replace("]", "").replace("'", "").replace(" ", "")
        final_gene_list = cleaned_gene_list.split("+")

        for gene_info in final_gene_list:
            res_organism, res_gene_name = gene_info.split("~")
            current_gene = Gene(res_organism, res_gene_name, get_info=False)
            temp_gene_list = genes_by_org[res_organism]
            temp_gene_list.append(current_gene)
            genes_by_org[res_organism] = temp_gene_list

    for organism, gene_list in genes_by_org.items():

        temp_set = set()
        for gene in gene_list:
            temp_set.add(gene.name)

        print(len(temp_set))

    print("All Organisms: ", sus_organisms)
    print("Number of all Genes: ", len(all_genes))

    threads = []
    for t in range(5):
        q = Queue()
        output_queue = Queue()
        threads.append(MyThread(q, output_queue, sus_organisms, drug_dirs))
        threads[t].start()
        time.sleep(0.1)

    thread_number = 0
    for organism, gene_list in genes_by_org.items():
        threads[thread_number].queue.put(gene_list)

        thread_number += 1
        if thread_number == 5:
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
DRUGS = ["SULF"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_recips(DRUG, PHENOTYPE)
