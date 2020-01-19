import os
import csv
import pandas as pd
import utils.dir_utils as dir_utils


class PhageData:

    def __init__(self, umb, is_infected):
        self.umb = umb
        if is_infected == "no":
            self.is_infected = 'R'
        else:
            self.is_infected = 'S'


class AntibioticData:

    def __init__(self, umb, cp):
        self.umb = umb
        self.cp = cp


def convert_antibiotic_data():
    phage_data = list(csv.reader(open(os.path.join(dir_utils.MAIN_DIR, "phage_resistance.csv")), delimiter=','))
    antibiotic_data = list(csv.reader(open(os.path.join(dir_utils.MAIN_DIR, "antibiotic_resistance_csv.csv")), delimiter=','))
    first_phage = True
    first_antibiotic = True
    phage_data_list = {}
    antibiotic_data_list = {}

    for row in phage_data:
        if first_phage:
            first_phage = False
            continue

        phage_data = PhageData(row[0], row[1])
        phage_data_list[phage_data.umb] = phage_data

    for row in antibiotic_data:
        if first_antibiotic:
            first_antibiotic = False
            continue

        antibiotic_data = AntibioticData(row[0], row[1])
        antibiotic_data_list[antibiotic_data.umb] = antibiotic_data

    print(len(phage_data_list))
    print(len(antibiotic_data_list))

    with open(os.path.join(dir_utils.MAIN_DIR, "converted_phage_data.csv"), "w") as converted_data:
        converted_data.write("UMB,CP #,phage\n")
        for key, data in phage_data_list.items():
            if key in antibiotic_data_list.keys():
                antibiotic_data = antibiotic_data_list[key]
                converted_data.write(f"{key},{antibiotic_data.cp},{data.is_infected}\n")


convert_antibiotic_data()
