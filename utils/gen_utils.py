import csv


def get_op_phenotype(phenotype):
    op_phenotype = "res"
    if phenotype == "res":
        op_phenotype = "sus"

    return op_phenotype


def get_organisms_by_phenotype(organism_path):
    all_organisms_csv = list(csv.reader(open(organism_path), delimiter=','))
    all_organisms = []
    for org in all_organisms_csv:
        all_organisms.append(org[0])
    return all_organisms
