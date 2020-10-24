import csv
import os
import threading
import time
from queue import Queue


from models import blast as b
from models.gene import Gene
from utils import gen_utils, output_util, dir_utils


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
            recip_organism_dirs = dir_utils.OrganismDirs(organism_name, converted_ncbi_data=True)
            database_path = recip_organism_dirs.database_dir

            # First blast the first organism gene against the database of the gene in the list.
            current_gene = Gene(self.target_organism, gene)
            blast_data = b.blast(current_gene, database_path, organism_name)

            if not blast_data:
                break

            if not blast_data.above_bitscore_threshold_check():
                break

            if not blast_data.within_size_constraints_check():
                break

            if gene in self.final_gene_info:
                old_data = self.final_gene_info[gene]
                gene_description = old_data[0]
                new_count = old_data[1] + 1
                new_qcov = (old_data[2] + blast_data.qcov) / 2
                new_pident = (old_data[3] + blast_data.pident) / 2
                new_genes = old_data[5]
                new_genes.append(f"{organism}~{blast_data.blast_gene.name}")
                if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                    new_perfect_matches = old_data[4] + 1
                else:
                    new_perfect_matches = old_data[4]

                self.final_gene_info[gene] = [gene_description, new_count, new_qcov, new_pident,
                                              new_perfect_matches, new_genes]
            else:
                hit_count = 1
                if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                    perfect_matches = 1
                else:
                    perfect_matches = 0
                self.final_gene_info[gene] = [current_gene.description, hit_count, blast_data.qcov,
                                              blast_data.pident, perfect_matches,
                                              [f"{self.target_organism}~{gene}",
                                               f"{organism}~{blast_data.blast_gene.name}"]]

        with print_lock:
            print(threading.currentThread().getName(), f"Genes left {self.queue.qsize()}")


def get_recips(drug, phenotype):
    # Create output file
    drug_dirs = dir_utils.DrugDirs(drug, phenotype)
    output_file = output_util.OutputFile(drug_dirs.res_recip_genes_file, header_list=["gene", "gene_description",
                                                                                      "hit_count", "qcov_average",
                                                                                      "pident_average",
                                                                                      "perfect_matches", "recip_genes"])
    # Get all the organisms for the phenotype
    res_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.res_file)
    target_org = str(res_organisms.pop(0))
    target_org_dirs = dir_utils.OrganismDirs(target_org, converted_ncbi_data=True)

    # get all the genes for the target organism
    all_genes = os.listdir(target_org_dirs.gene_folder)

    print("All Organisms: ", res_organisms)
    print("Number of all Genes: ", len(all_genes))
    print("Original Organism: ", target_org)

    threads = []
    for t in range(10):
        q = Queue()
        output_queue = Queue()
        threads.append(MyThread(q, output_queue, target_org, res_organisms))
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
