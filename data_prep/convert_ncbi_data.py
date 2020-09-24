import os
import shutil

import pandas as pd
from Bio import SeqIO

from utils import dir_utils


organism_labels = pd.read_csv(dir_utils.ORGANISM_LABELS_FILE, header=0)
CONVERSIONS = dict(zip(list(organism_labels.iloc[:, 2]), list(organism_labels.iloc[:, 0])))
dir_utils.generate_dir(dir_utils.CONVERTED_NCBI_DATA_DIR)

for org_file in os.listdir(dir_utils.NCBI_DATA_DIR):
    org_file_path = os.path.join(dir_utils.NCBI_DATA_DIR, org_file)
    ncbi_label = org_file.split(".")[0]
    if ncbi_label == "converted_ncbi_data":
        continue

    print(CONVERSIONS[ncbi_label])
    if "UMB" in CONVERSIONS[ncbi_label]:
        converted_name = CONVERSIONS[ncbi_label]
    else:
        converted_name = f"CP-{CONVERSIONS[ncbi_label]}"

    converted_org_path = os.path.join(dir_utils.CONVERTED_NCBI_DATA_DIR, converted_name)
    gene_folder_path = os.path.join(converted_org_path, "genes")
    dir_utils.generate_dir(converted_org_path)
    dir_utils.generate_dir(gene_folder_path)

    os.chdir(dir_utils.NCBI_DATA_DIR)
    shutil.copy(org_file_path, os.path.join(converted_org_path, f"{converted_name}.faa"))
    os.chdir(converted_org_path)
    os.system(f'makeblastdb -in {os.path.join(converted_org_path, f"{converted_name}.faa")} '
              f'-title {converted_name} -out {converted_name} -dbtype prot')

    for record in SeqIO.parse(org_file_path, "fasta"):
        with open(os.path.join(gene_folder_path, f"{record.id}.fasta"), "w") as new_record_file:
            new_record_file.write(f"> {record.description.replace(',', '_')}\n")
            new_record_file.write(f"{record.seq}")
