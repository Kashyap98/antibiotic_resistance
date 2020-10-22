import csv
import os
import utils.dir_utils as dir_utils
import pandas as pd


def cp_fixer(cp):
    if cp.startswith("UMB"):
        return cp
    else:
        return f"CP-{cp}"


def organize_organisms(cp, antibiotic, antibiotic_type):
    print(f"{cp}: {antibiotic_type} - {antibiotic}")

    if antibiotic_type == "fosf":
        DIR = dir_utils.FOSF_DIR
    elif antibiotic_type == "cipro":
        DIR = dir_utils.CIPRO_DIR
    elif antibiotic_type == "amoxo":
        DIR = dir_utils.AMOXO_DIR
    elif antibiotic_type == "cepro":
        DIR = dir_utils.CEPRO_DIR
    else:
        DIR = dir_utils.SULF_DIR

    # Sort organism by phenotype
    if antibiotic == "R":
        res_file = open(os.path.join(DIR, "res.csv"), "a")
        res_file.write(f"{cp}\n")
        res_file.close()
    elif antibiotic == "S":
        sus_file = open(os.path.join(DIR, "sus.csv"), "a")
        sus_file.write(f"{cp}\n")
        sus_file.close()
    elif antibiotic == "I":
        ind_file = open(os.path.join(DIR, "ind.csv"), "a")
        ind_file.write(f"{cp}\n")
        ind_file.close()


def sort_organisms_by_phenotype():
    dir_utils.generate_dir(dir_utils.FOSF_DIR)
    dir_utils.generate_dir(dir_utils.CIPRO_DIR)
    dir_utils.generate_dir(dir_utils.AMOXO_DIR)
    dir_utils.generate_dir(dir_utils.CEPRO_DIR)
    dir_utils.generate_dir(dir_utils.SULF_DIR)

    data = list(csv.reader(open(dir_utils.DATA_SOURCE_FILE, 'r'), delimiter=','))
    df = pd.DataFrame.from_records(data)

    folders = len(df.index)
    count = 0
    while count < folders:
        # organism each organism
        organism = df.iloc[count]
        cp = cp_fixer(organism[1])
        fosf = organism[2]
        cipro = organism[3]
        amoxo = organism[4]
        sulf = organism[5]
        cepro = organism[6]
        organize_organisms(cp, fosf, "fosf")
        organize_organisms(cp, cipro, "cipro")
        organize_organisms(cp, amoxo, "amoxo")
        organize_organisms(cp, sulf, "sulf")
        organize_organisms(cp, cepro, "cepro")
        count += 1
        print(f"Progress: {count}")
