import csv
import glob
import os
import threading
import time
from queue import Queue

import pandas as pd

from models import Blast as b
from models.Gene import Gene
from utils.output_util import DataOutput
from utils.dir_utils import OrganismDirs, DrugDirs


print_lock = threading.Lock()


class MyThread(threading.Thread):
    def __init__(self, input_queue, output_queue, target_organism, res_organisms, sus_organisms):
        super().__init__()
        self.queue = input_queue
        self.output_queue = output_queue
        self.daemon = True
        self.target_organism = target_organism
        self.res_organisms = res_organisms
        self.sus_organisms = sus_organisms
        self.final_gene_info = {}

    def run(self):
        while True:
            gene = self.queue.get()
            if gene is None:  # If you send `None`, the thread will exit.
                self.output_queue.put(self.final_gene_info)
                return
            self.do_work(gene)

    def do_work(self, gene):
        current_gene = Gene(self.target_organism, gene)
        res_genes = [current_gene]

        for organism in self.res_organisms:
            # For each other genome we want to compare to, Recip BLAST those genes and only keep the ones are the same
            organism_name = str(organism)
            recip_organism_dirs = OrganismDirs(organism_name, converted_ncbi_data=True)
            database_path = recip_organism_dirs.database_dir

            # First blast the first organism gene against the database of the gene in the list.
            blast_data = b.blast(current_gene, database_path, organism_name)

            if not blast_data:
                self.final_gene_info.pop(gene)
                break

            res_genes.append(blast_data.blast_gene)

            if gene in self.final_gene_info:
                old_data = self.final_gene_info[gene]
                gene_description = old_data[0]
                new_qcov = (old_data[1] + blast_data.qcov) / 2
                new_pident = (old_data[2] + blast_data.pident) / 2
                new_count = old_data[3] + 1
                if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                    new_perfect_matches = old_data[4] + 1
                    new_match_organisms = old_data[5]
                    new_perfect_match_organisms = old_data[6]
                    new_perfect_match_organisms.append(organism)
                else:
                    new_perfect_matches = old_data[4]
                    new_perfect_match_organisms = old_data[6]
                    new_match_organisms = old_data[5]
                    new_match_organisms.append(organism)

                self.final_gene_info[gene] = [gene_description, new_qcov, new_pident, new_count, new_perfect_matches,
                                              new_match_organisms, new_perfect_match_organisms, 0, 0, 0, 0, [], []]
            else:
                hit_count = 1
                if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                    perfect_matches = 1
                    res_match_organisms = []
                    res_perfect_match_organisms = [organism_name]
                else:
                    perfect_matches = 0
                    res_perfect_match_organisms = []
                    res_match_organisms = [organism_name]

                self.final_gene_info[gene] = [current_gene.description, blast_data.qcov, blast_data.pident, hit_count,
                                              perfect_matches, res_match_organisms, res_perfect_match_organisms, 0, 0,
                                              0, 0, [], []]

        if gene in self.final_gene_info.keys():

            for current_gene in res_genes:
                for organism in self.sus_organisms:
                    organism_name = str(organism)
                    recip_organism_dirs = OrganismDirs(organism_name, converted_ncbi_data=True)
                    database_path = recip_organism_dirs.database_dir

                    # First blast the first organism gene against the database of the gene in the list.
                    blast_data = b.blast(current_gene, database_path, organism_name)

                    if not blast_data:
                        continue

                    if gene in self.final_gene_info:
                        old_data = self.final_gene_info[gene]
                        gene_description = old_data[0]
                        new_qcov = (old_data[1] + blast_data.qcov)
                        new_pident = (old_data[2] + blast_data.pident)
                        new_count = old_data[3]
                        new_perfect_matches = old_data[4]
                        new_match_organisms = old_data[5]
                        new_perfect_match_organisms = old_data[6]

                        if old_data[7] == 0:
                            new_sus_qcov = blast_data.qcov
                        else:
                            new_sus_qcov = (old_data[7] + blast_data.qcov) / 2

                        if old_data[8] == 0:
                            new_sus_pident = blast_data.pident
                        else:
                            new_sus_pident =  (old_data[8] + blast_data.pident) / 2

                        new_sus_matches = old_data[9] + 1

                        if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                            new_sus_perfect_matches = old_data[10] + 1
                            new_sus_match_organisms = old_data[11]
                            new_sus_perfect_match_organisms = old_data[12]
                            new_sus_perfect_match_organisms.append(organism)
                        else:
                            new_sus_perfect_matches = old_data[10]
                            new_sus_perfect_match_organisms = old_data[12]
                            new_sus_match_organisms = old_data[11]
                            new_sus_match_organisms.append(organism)

                        self.final_gene_info[gene] = [gene_description, new_qcov, new_pident, new_count,
                                                      new_perfect_matches,
                                                      new_match_organisms, new_perfect_match_organisms, new_sus_qcov,
                                                      new_sus_pident, new_sus_matches, new_sus_perfect_matches,
                                                      new_sus_match_organisms, new_sus_perfect_match_organisms]

        with print_lock:
            print(threading.currentThread().getName(), f"Genes left {self.queue.qsize()}")


def get_recips(drug, phenotype):
    # Create output file
    drug_dirs = DrugDirs(drug, phenotype)
    with open(os.path.join(drug_dirs.drug_dir, f"{phenotype}_complete_comparison.csv"), "w") as recip_detailed:
        recip_detailed.write("gene,gene_description,res_qcov_average,res_pident_average,res_match_count,res_perfect_match_count,res_matches,res_perfect_matches,"
                             "sus_qcov_average,sus_pident_average,sus_match_count,sus_perfect_match_count,sus_matches,sus_perfect_matches\n")
    # Get all the organisms for the phenotype
    all_organisms_raw = list(csv.reader(open(drug_dirs.target_phenotype_file), delimiter=','))
    res_organisms = []
    for org in all_organisms_raw:
        res_organisms.append(org[0])
    target_org = str(res_organisms.pop(0))
    target_org_dirs = OrganismDirs(target_org, converted_ncbi_data=True)
    op_organism_file = pd.read_csv(drug_dirs.op_phenotype_file, header=None)
    op_orgs = list(op_organism_file.iloc[:, 0])
    all_genes_path = target_org_dirs.gene_folder

    # get all the genes for the target organism
    all_genes = os.listdir(all_genes_path)

    print("All Organisms: ", res_organisms)
    print("Number of all Genes: ", len(all_genes))
    print("Original Organism: ", target_org)

    threads = []
    for t in range(6):
        q = Queue()
        output_queue = Queue()
        threads.append(MyThread(q, output_queue, target_org, res_organisms, op_orgs))
        threads[t].start()
        time.sleep(0.1)

    thread_number = 0
    for x in all_genes:
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
        while not t.output_queue.empty():
            thread_data = t.output_queue.get()
            with open(os.path.join(drug_dirs.drug_dir, f"{phenotype}_complete_comparison.csv"), "a") as recip_detailed:
                for gene, data in thread_data.items():
                    recip_detailed.write(f"{gene},{str(data[0])},{str(data[1])},{str(data[2])},{str(data[3])},"
                                         f"{str(data[4])},{str(data[5]).replace(',', '.')},"
                                         f"{str(data[6]).replace(',', '.')},{str(data[7])},{str(data[8])},"
                                         f"{str(data[9])},{str(data[10])},{str(data[11]).replace(',', '.')},"
                                         f"{str(data[12]).replace(',', '.')}\n")
            print(data)

        t.queue.put(None)
        t.join()


PHENOTYPES = ["res"]
DRUGS = ["AMOXO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_recips(DRUG, PHENOTYPE)
