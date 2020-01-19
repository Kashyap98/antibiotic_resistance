import os
import csv
import utils.dir_utils as dir_utils


def sort_phage_data():
    phage_data = list(csv.reader(open(os.path.join(dir_utils.MAIN_DIR, "converted_phage_data.csv")), delimiter=','))
    org_dirs = []
    first_row = True
    for row in phage_data:
        if first_row:
            first_row = False
            orgs = row[2:]
            for org in orgs:
                org_dir = os.path.join(dir_utils.SORTED_DATA_DIR, org)
                dir_utils.generate_dir(org_dir)
                org_dirs.append(org_dir)
            continue

        cp = row[1]
        pheno_data = row[2:]
        count = 0
        for data in pheno_data:
            org_dir = org_dirs[count]
            if data == "R":
                with open(os.path.join(org_dir, "res.csv"), "a") as res_file:
                    res_file.write(f"{cp}\n")
            else:
                with open(os.path.join(org_dir, "sus.csv"), "a") as sus_file:
                    sus_file.write(f"{cp}\n")
            count += 1


sort_phage_data()