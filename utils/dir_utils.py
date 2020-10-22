import os

from utils import gen_utils

MAIN_DIR = os.path.join(os.getcwd(), "..")
ANNOTATION_DIR = os.path.join(MAIN_DIR, "annotations")
CONVERTED_DATA_DIR = os.path.join(MAIN_DIR, "converted_data")
SORTED_DATA_DIR = os.path.join(MAIN_DIR, "sorted_data")
NCBI_DATA_DIR = os.path.join(MAIN_DIR, "ncbi_data")
CONVERTED_NCBI_DATA_DIR = os.path.join(NCBI_DATA_DIR, "converted_ncbi_data")

# Drug constants
FOSF_DIR = os.path.join(SORTED_DATA_DIR, "FOSF")
CIPRO_DIR = os.path.join(SORTED_DATA_DIR, "CIPRO")
AMOXO_DIR = os.path.join(SORTED_DATA_DIR, "AMOXO")
SULF_DIR = os.path.join(SORTED_DATA_DIR, "SULF")
CEPRO_DIR = os.path.join(SORTED_DATA_DIR, "CEPRO")

ORGANISM_LABELS_FILE = os.path.join(MAIN_DIR, "organism_labels.csv")
DATA_SOURCE_FILE = os.path.join(MAIN_DIR, "antibiotic_resistance_csv.csv")


def generate_dir(in_dir):
    try:
        os.makedirs(in_dir)
    except FileExistsError:
        print(f"Folder already made! - {in_dir}")


def cp_umb_dir_fixer(cp):
    if str(cp).isnumeric():
        return f"CP-{cp}"
    else:
        return cp


class OrganismDirs:

    def __init__(self, organism, converted_ncbi_data=False):
        self.organism = organism
        if converted_ncbi_data:
            self.organism_folder = os.path.join(CONVERTED_NCBI_DATA_DIR, cp_umb_dir_fixer(organism))
        else:
            self.organism_folder = os.path.join(CONVERTED_DATA_DIR, cp_umb_dir_fixer(organism))

        self.gene_folder = os.path.join(self.organism_folder, "genes")
        self.gene_info_folder = os.path.join(self.organism_folder, "gene_info")
        self.database_dir = os.path.join(self.organism_folder, cp_umb_dir_fixer(organism))


class DrugDirs:

    def __init__(self, drug, phenotype):
        self.drug_dir = os.path.join(SORTED_DATA_DIR, drug)
        self.res_file = os.path.join(self.drug_dir, "res.csv")
        self.sus_file = os.path.join(self.drug_dir, "sus.csv")
        self.ind_file = os.path.join(self.drug_dir, "ind.csv")

        self.target_phenotype = phenotype
        self.opposite_phenotype = gen_utils.get_op_phenotype(phenotype)

        self.target_phenotype_file = os.path.join(self.drug_dir, f"{self.target_phenotype}.csv")
        self.op_phenotype_file = os.path.join(self.drug_dir, f"{self.opposite_phenotype}.csv")

        self.res_recip_genes_file = os.path.join(self.drug_dir, f"res_recip_genes.csv")
        self.res_unique_genes_file = os.path.join(self.drug_dir, f"res_unique_genes.csv")
