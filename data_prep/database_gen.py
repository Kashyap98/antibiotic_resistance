import os
import utils.dir_utils as dir_utils


def create_blast_databases_for_organisms(database_type="prot"):
    # Used to create the BLAST Database for each organism.
    for organism in os.listdir(dir_utils.CONVERTED_DATA_DIR):
        print(organism)

        current_organism = os.path.join(dir_utils.CONVERTED_DATA_DIR, organism)
        os.system(f"type {os.path.join(current_organism, 'genes', '*.fasta')} > "
                  f"{os.path.join(current_organism, 'all_fasta.fasta')}")
        os.system(f'makeblastdb -in all_fasta.fasta -title {organism} -out {organism} -dbtype {database_type}')
