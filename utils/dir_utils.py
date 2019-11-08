import os

MAIN_DIR = os.path.join(os.getcwd(), "..")
ANNOTATION_DIR = os.path.join(MAIN_DIR, "annotations")
CONVERTED_DATA_DIR = os.path.join(MAIN_DIR, "converted_data")
SORTED_DATA_DIR = os.path.join(MAIN_DIR, "sorted_data")

# Drug constants
FOSF_DIR = os.path.join(SORTED_DATA_DIR, "FOSF")
CIPRO_DIR = os.path.join(SORTED_DATA_DIR, "CIPRO")
AMOXO_DIR = os.path.join(SORTED_DATA_DIR, "AMOXO")
SULF_DIR = os.path.join(SORTED_DATA_DIR, "SULF")
CEPRO_DIR = os.path.join(SORTED_DATA_DIR, "CEPRO")


class OrganismDirs:

    def __init__(self, organism):
        self.organism_folder = os.path.join(CONVERTED_DATA_DIR, organism)
        self.gene_folder = os.path.join(self.organism_folder, "genes")
        self.gene_info_folder = os.path.join(self.organism_folder, "gene_info")
        self.database_dir = os.path.join(self.organism_folder, organism)


class DrugDirs:

    def __init__(self, drug, phenotype):
        self.drug_dir = os.path.join(SORTED_DATA_DIR, drug, phenotype)
        self.res_file = os.path.join(self.drug_dir, "res.csv")
        self.sus_file = os.path.join(self.drug_dir, "sus.csv")
        self.ind_file = os.path.join(self.drug_dir, "ind.csv")
        self.target_phenotype_file = os.path.join(self.drug_dir, f"{phenotype}.csv")
