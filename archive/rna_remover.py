import os
import shutil
import utils.dir_utils as dir_utils


# Used to create the BLASTn Database for each organism.
for organism in os.listdir(dir_utils.CONVERTED_DATA_DIR):
    print(organism)
    currentOrganism = os.path.join(dir_utils.CONVERTED_DATA_DIR, organism, "genes")
    for gene in os.listdir(currentOrganism):
        if "rna" in gene:
            os.remove(os.path.join(currentOrganism, gene))
