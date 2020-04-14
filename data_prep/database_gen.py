import os
import shutil
import utils.dir_utils as dir_utils


# Used to create the BLASTn Database for each organism.
for organism in os.listdir(dir_utils.CONVERTED_DATA_DIR):
    print(organism)
    currentOrganism = os.path.join(dir_utils.CONVERTED_DATA_DIR, organism)
    if os.path.exists(os.path.join(dir_utils.CONVERTED_DATA_DIR, organism, "AllFasta.fasta")):
        os.remove(os.path.join(dir_utils.CONVERTED_DATA_DIR, organism, "AllFasta.fasta"))
    if os.path.exists(os.path.join(dir_utils.CONVERTED_DATA_DIR, organism, f"{organism}.nhr")):
        os.remove(os.path.join(dir_utils.CONVERTED_DATA_DIR, organism, f"{organism}.nhr"))
    if os.path.exists(os.path.join(dir_utils.CONVERTED_DATA_DIR, organism, f"{organism}.nin")):
        os.remove(os.path.join(dir_utils.CONVERTED_DATA_DIR, organism, f"{organism}.nin"))
    if os.path.exists(os.path.join(dir_utils.CONVERTED_DATA_DIR, organism, f"{organism}.nsq")):
        os.remove(os.path.join(dir_utils.CONVERTED_DATA_DIR, organism, f"{organism}.nsq"))
    os.chdir(os.path.join(currentOrganism, "genes"))
    os.system("type *.fasta > AllFasta.fasta")
    shutil.move(os.path.join(os.getcwd(), "AllFasta.fasta"), currentOrganism)
    os.chdir(currentOrganism)
    os.system(f'makeblastdb -in AllFasta.fasta -title {organism} -out {organism} -dbtype nucl')
