import os

import pandas as pd


def gene_info_writer(path, gene_name, organism_name, outputFile):
    geneInfo = pd.read_csv(path + "..\\gene_info\\" + gene_name + ".txt", error_bad_lines=False)
    geneData = list(geneInfo.iloc[:, 1])
    outputFile.write(organism_name + "," + gene_name + "," + geneInfo.columns.values.tolist()[1]
                     + ",")

    for x in range(1, 9):
        outputFile.write(str(geneData[x]) + ",")

    outputFile.write("\n")


# This is used to get all the genes in common for a phenotype
def get_gene_info(drug, phenotype):
    outputFile = open(os.getcwd() + "\\sorted_data\\" + drug + "\\" + phenotype + "_UniqueGeneInfo.csv", "w")
    outputFile.write("organism, feature_id, contig_id, type_data, location, start, stop, strand, function_data,"
                     "aliases, figfam, evidence_codes \n")

    # Main Organism
    organismFile = pd.read_csv(os.getcwd() + "\\sorted_data\\" + drug + "\\" + phenotype + ".csv")
    allOrgs = list(organismFile.iloc[:, 0])
    targetOrganism = str(allOrgs.pop(0))

    ## Get file with all the Genes
    geneFile = pd.read_csv(os.getcwd() + "\\sorted_data\\" + drug + "\\" + phenotype + "_UniqueGenes.csv")
    allGenes = list(geneFile.iloc[:, 0])
    print("Target Organism: ", targetOrganism)

    if targetOrganism.startswith("UMB"):
        allGenesPath = os.getcwd() + "\\converted_data\\" + targetOrganism + "\\genes\\"
    else:
        allGenesPath = os.getcwd() + "\\converted_data\\CP-" + targetOrganism + "\\genes\\"

    print("Number of all Genes: ", len(allGenes))

    for gene in allGenes:
        gene_info_writer(allGenesPath, gene.replace(".fasta", ""), targetOrganism, outputFile)


phenotypes = ["sus"]
drugs = ["AMOXO"]

for drug in drugs:
    for phenotype in phenotypes:
        get_gene_info(drug, phenotype)
