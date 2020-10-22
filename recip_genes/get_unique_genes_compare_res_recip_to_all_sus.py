import csv
import glob
import os
import threading
import time
from queue import Queue

import pandas as pd

from models import blast as b
from models.gene import Gene
from utils import gen_utils
from utils.output_util import DataOutput
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
            gene = self.queue.get()
            if gene is None:  # If you send `None`, the thread will exit.
                self.output_queue.put(self.final_gene_info)
                return
            self.do_work(gene)

    def do_work(self, gene):

        organism_name = str(self.target_organism)
        recip_organism_dirs = OrganismDirs(organism_name, converted_ncbi_data=True)
        database_path = recip_organism_dirs.database_dir

        # First blast the first organism gene against the database of the gene in the list.
        target_gene = Gene(self.target_organism, gene)

        for res_organism in self.organisms:
            res_organism_dirs = OrganismDirs(res_organism, converted_ncbi_data=True)
            res_database_path = res_organism_dirs.database_dir

            res_blast_data = b.blast(target_gene, res_database_path, res_organism)
            current_gene = res_blast_data.blast_gene

            if not res_blast_data:
                continue

            for op_org in self.op_orgs:

                sus_org_dirs = OrganismDirs(op_org, converted_ncbi_data=True)
                sus_database_path = sus_org_dirs.database_dir

                blast_data = b.blast(current_gene, sus_database_path, op_org)

                if not blast_data:
                    continue

                if gene in self.final_gene_info:
                    old_data = self.final_gene_info[gene]
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

                    self.final_gene_info[gene] = [gene_description, new_qcov, new_pident, new_count,
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

                    self.final_gene_info[gene] = [current_gene.description, blast_data.qcov,
                                                  blast_data.pident, hit_count, perfect_matches, res_match_organisms,
                                                  res_perfect_match_organisms]

        with print_lock:
            print(threading.currentThread().getName(), f"Genes left {self.queue.qsize()}")


def get_recips(drug, phenotype):
    # Create output file
    drug_dirs = DrugDirs(drug, phenotype)
    with open(os.path.join(drug_dirs.drug_dir, f"{phenotype}_aa_sus_detailed_info_ALL_TEST.csv"), "w") as recip_detailed:
        recip_detailed.write("gene,gene_description,qcov_average,pident_average,hit_count,perfect_match_count,hit_orgs,perfect_match_orgs\n")
    # Get all the organisms for the phenotype
    all_organisms_raw = list(csv.reader(open(drug_dirs.target_phenotype_file), delimiter=','))
    all_organisms = []
    for org in all_organisms_raw:
        all_organisms.append(org[0])
    target_org = str(all_organisms.pop(0))

    op_phenotype = gen_utils.get_op_phenotype(phenotype)
    drug_dirs.set_opposite_phenotype_file(op_phenotype)
    op_organism_file = pd.read_csv(drug_dirs.op_phenotype_file, header=None)
    op_orgs = list(op_organism_file.iloc[:, 0])

    all_genes_df = pd.read_csv(os.path.join(drug_dirs.drug_dir, f"{phenotype}_aa_recip_detailed_info.csv"), header=0)
    filtered_df = all_genes_df[all_genes_df["hit_count"] == 11]

    # get all the genes for the target organism
    all_genes = list(filtered_df["gene"])

    print("All Organisms: ", op_orgs)
    print("Number of all Genes: ", len(all_genes))
    print("Original Organism: ", target_org)

    threads = []
    for t in range(10):
        q = Queue()
        output_queue = Queue()
        threads.append(MyThread(q, output_queue, target_org, all_organisms, op_orgs))
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
            with open(os.path.join(drug_dirs.drug_dir, f"{phenotype}_aa_sus_detailed_info_ALL_TEST.csv"), "a") as recip_detailed:
                for gene, data in thread_data.items():
                    recip_detailed.write(f"{gene},{str(data[0])},{str(data[1])},{str(data[2])},{str(data[3])},"
                                         f"{str(data[4])},{str(data[5]).replace(',', '.')},"
                                         f"{str(data[6]).replace(',', '.')}\n")

        t.queue.put(None)
        t.join()


PHENOTYPES = ["res"]
DRUGS = ["CIPRO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_recips(DRUG, PHENOTYPE)
