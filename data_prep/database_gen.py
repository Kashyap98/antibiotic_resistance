import os
import shutil
import utils.dir_utils as dir_utils


# Used to create the BLASTn Database for each organism.
for organism in os.listdir(dir_utils.CONVERTED_DATA_DIR):
    print(organism)
    currentOrganism = os.path.join(dir_utils.CONVERTED_DATA_DIR, organism)
    os.chdir(os.path.join(currentOrganism, "genes"))
    os.system("type *.fasta> AllFasta.fasta")
    shutil.move(os.path.join(os.getcwd(), "AllFasta.fasta"), currentOrganism)
    os.chdir(currentOrganism)
    os.system('makeblastdb -in AllFasta.fasta -title ' + organism + ' -out ' + organism + ' -dbtype nucl')
