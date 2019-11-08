import os

import pandas as pd

from models import Blast as b
from models.BlastResult import BlastResult
from models.Gene import Gene


def write_cell(outputFile, data, isStart=False):
    # write to the file based on the type of data being inputted
    if isStart:
        outputFile.write(data + ",")
    else:
        outputFile.write(str(data) + ",")


def gene_info_writer(targetGene, recipGene, bitscore, matchRatio, outputFile):
    # write each cell of the csv
    write_cell(outputFile, targetGene.organism, isStart=True)
    write_cell(outputFile, targetGene.name)
    write_cell(outputFile, recipGene.organism)
    write_cell(outputFile, recipGene.name)
    write_cell(outputFile, bitscore)
    write_cell(outputFile, matchRatio)
    write_cell(outputFile, "Remember to fix")
    write_cell(outputFile, recipGene.type_data)
    write_cell(outputFile, recipGene.location)
    write_cell(outputFile, recipGene.start)
    write_cell(outputFile, recipGene.stop)
    write_cell(outputFile, recipGene.strand)
    write_cell(outputFile, recipGene.function_data)
    write_cell(outputFile, recipGene.aliases)
    write_cell(outputFile, recipGene.figfam)
    write_cell(outputFile, recipGene.evidence_codes)
    outputFile.write("\n")


## This is used to find the genes that are unique to each phenotype after being recip blasted to the same phenotype
def get_uniques(drug, phenotype):
    # create output file and write the headers
    outputFile = open(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_UniqueBlastInfo.csv"), "w")
    outputFile.write("organism, feature_id, recip_organism, recip_feature_id, bitscore, match_ratio,"
                     "contig_id, type_data, location, start, stop, strand, function_data,"
                     "aliases, figfam, evidence_codes \n")
    op_phenotype = "res"
    if phenotype == "res":
        op_phenotype = "sus"

    # Main Organism
    organismFile = pd.read_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + ".csv"))
    allOrgs = list(organismFile.iloc[:, 0])
    targetOrganism = str(allOrgs[0])

    # Get file with all the Genes
    geneFile = pd.read_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_RecipGenes.csv"))
    allGenes = list(geneFile.iloc[:, 0])

    print("Number of all Genes: ", len(allGenes))

    uniqueMatched = {}
    op_organismFile = pd.read_csv(os.path.join(os.getcwd(), "sorted_data", drug, op_phenotype + ".csv"))
    opOrgs = list(op_organismFile.iloc[:, 0])
    count = 0
    for organism in opOrgs:
        # For each other genome we want to compare to, Recip BLAST those genes and only keep the ones are the same
        organism = str(organism)
        print(str(count) + "/" + str(len(opOrgs)))

        databasePath = os.path.join(os.getcwd(), "converted_data", organism, organism)

        for gene in allGenes:
            # First blast the first organism gene against the database of the gene in the list.
            targetGene = Gene(targetOrganism, str(gene.replace(".fasta", "")), get_info=True)
            blastData = b.blast(targetGene.fasta_file, databasePath)
            # print(blastData)
            matchedOrganisms = []
            matchedGenes = []

            # If we have already found it in unique_matched
            if gene in uniqueMatched:
                matchedOrganisms = uniqueMatched[gene][2]
                matchedGenes = uniqueMatched[gene][3]

            # if there is a hit
            if len(blastData) > 1:
                blastResult = BlastResult(blastData)
                recipGene = Gene(organism, blastResult.gene_name, get_info=True)
                matchRatio = float(round(int(blastResult.match_length) / int(targetGene.length), 2))
                gene_info_writer(targetGene, recipGene, blastResult.bitscore, matchRatio, outputFile)

                # if the result is unique enought
                if blastResult.bitscore < 1600:
                    matchedGenes.append(blastResult.gene_name)
                    matchedOrganisms.append(organism)
                    uniqueMatched[gene] = (targetOrganism, len(matchedOrganisms), matchedOrganisms, matchedGenes)
                else:
                    allGenes.remove(gene)
            else:
                allGenes.remove(gene)

        print("Organism: " + organism)
        print("Matched Genes Length: ", len(allGenes))
        count += 1

    # output the results as a csv.
    unique_df = pd.DataFrame()
    temp_count = 0
    for key in uniqueMatched.keys():
        row = uniqueMatched[key]
        series = pd.Series([key, row[0], row[1], row[2], row[3]])
        unique_df = unique_df.append(series, ignore_index=True)
        temp_count += 1

    unique_df.to_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_UniqueMatches.csv"), index=False)


phenotypes = ["sus", "res"]
drugs = ["CIPRO"]

for drug in drugs:
    for phenotype in phenotypes:
        get_uniques(drug, phenotype)
