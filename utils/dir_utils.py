import os
import shutil


MAIN_DIR = os.path.join(os.getcwd(), "..")
ANNOTATION_DIR = os.path.join(MAIN_DIR, "annotations")
CONVERTED_DATA_DIR = os.path.join(MAIN_DIR, "converted_data")
SORTED_DATA_DIR = os.path.join(MAIN_DIR, "sorted_data")
NCBI_DATA_DIR = os.path.join(MAIN_DIR, "ncbi_data")

# Drug constants
FOSF_DIR = os.path.join(SORTED_DATA_DIR, "FOSF")
CIPRO_DIR = os.path.join(SORTED_DATA_DIR, "CIPRO")
AMOXO_DIR = os.path.join(SORTED_DATA_DIR, "AMOXO")
SULF_DIR = os.path.join(SORTED_DATA_DIR, "SULF")
CEPRO_DIR = os.path.join(SORTED_DATA_DIR, "CEPRO")

ORGANISM_LABELS_FILE = os.path.join(MAIN_DIR, "organism_labels.csv")
DATA_SOURCE_FILE = os.path.join(MAIN_DIR, "antibiotic_resistance_csv.csv")
NEW_DATA_SOURCE_FILE = os.path.join(MAIN_DIR, "ResistanceInfo.csv")
USEARCH = os.path.join(MAIN_DIR, "usearch.exe")


def generate_dir(in_dir, overwrite_dir=False):
    if overwrite_dir:
        if os.path.exists(in_dir):
            shutil.rmtree(in_dir)

    try:
        os.makedirs(in_dir)
    except FileExistsError:
        print(f"Folder already made! - {in_dir}")


def cp_umb_dir_fixer(cp):
    if str(cp).isnumeric():
        return f"CP-{cp}"
    else:
        return cp


def strip_extension(file_name: str) -> str:
    return file_name.split(".")[0]


class OrganismDirs:

    def __init__(self, organism: str, new_organism: bool = False):
        self.organism = organism

        if new_organism:
            self.organism_folder = os.path.join(CIPRO_DIR, "new_org_dir")
            self.gene_folder = self.organism_folder
            self.gene_info_folder = self.organism_folder
            self.database_dir = self.organism_folder
        else:
            self.organism_folder = os.path.join(CONVERTED_DATA_DIR, cp_umb_dir_fixer(organism))
            self.gene_folder = os.path.join(self.organism_folder, "genes")
            self.gene_info_folder = os.path.join(self.organism_folder, "gene_info")
            self.database_dir = os.path.join(self.organism_folder, cp_umb_dir_fixer(organism))


class DrugDirs:

    def __init__(self, drug, phenotype):
        self.drug_dir = os.path.join(SORTED_DATA_DIR, drug)
        self.new_organism_dir = os.path.join(self.drug_dir, "new_org_dir")
        self.unique_res_genes = os.path.join(self.drug_dir, "unique_res_genes")
        self.reciprocal_res_genes = os.path.join(self.drug_dir, "reciprocal_res_genes")
        self.res_file = os.path.join(self.drug_dir, "res.csv")
        self.sus_file = os.path.join(self.drug_dir, "sus.csv")
        self.ind_file = os.path.join(self.drug_dir, "ind.csv")

        self.potential_uniques = os.path.join(self.drug_dir, f"potential_uniques.csv")
        self.investigated_unique_genes = os.path.join(self.drug_dir, "investigated_unique_genes.csv")
        self.cluster_info = os.path.join(self.drug_dir, "cluster_info.csv")
        self.unique_clusters = os.path.join(self.drug_dir, "unique_clusters.csv")

        self.sequences_to_cluster = os.path.join(self.drug_dir, f"sequences_to_cluster.fasta")
        self.centroids = os.path.join(self.drug_dir, f"centroids.fasta")
        self.cluster_uc = os.path.join(self.drug_dir, f"clusters.uc")
        self.cluster_dir = os.path.join(self.drug_dir, f"clusters")
