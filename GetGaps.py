import os

import pandas as pd

from models.Gene import Gene
from models.Mutation import Mutation


# Same as get mutations but is only gaps
def get_mutations(drug, phenotype):

    # Main Organism
    dataFile = pd.read_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_UniqueMatches.csv"))
    dataFile = dataFile.loc[dataFile['2'] == dataFile['2'].max()]
    originalGenes = list(dataFile.iloc[:, 0])
    originalOrganism = dataFile.iloc[1][1]
    matchedOrganisms = dataFile.iloc[3][3]
    print(originalOrganism)
    print(matchedOrganisms)
    matchedOrganisms = matchedOrganisms.strip("[").strip("]").replace("'", "").replace(" ", "")
    matchedOrganisms = matchedOrganisms.split(",")

    output = {}
    count = 1
    for gene in originalGenes:

        originalGene = Gene(originalOrganism, gene.replace(".fasta", ""), get_info=True)
        print(str(count) + "/" + str(len(originalGenes)))
        for row in dataFile["4"].loc[dataFile["0"] == gene]:
            matchedGenes = row.strip("[").strip("]").replace("'", "").replace(" ", "")
            matchedGenes = matchedGenes.split(",")

            for x in range(0, len(matchedOrganisms)):

                matchedGene = Gene(matchedOrganisms[x], matchedGenes[x], get_info=True)
                mutationInfo = Mutation(originalGene.aa_sequence, matchedGene.aa_sequence)
                output[matchedGenes[x]] = [matchedOrganisms[x], gene, originalOrganism, mutationInfo.score,
                                           mutationInfo.percentage]

                for mutation in mutationInfo.get_gaps():
                    output[matchedGenes[x]].append(mutation)
        count += 1
    unique_df = pd.DataFrame()

    for key in output.keys():
        row = output[key]
        series = pd.Series([key, row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10]
                               , row[11], row[12]])
        unique_df = unique_df.append(series, ignore_index=True)

    unique_df.rename(columns={0: "Matched Gene", 1: "Matched Organism", 2: "Original Gene", 3: "Original Organism",
                              4: "Mutation Score", 5: "Mutation Percentage", 6: "nonpolar_reference",
                              7: "polar_reference", 8: "acidic_reference", 9: "basic_refernce", 10: "nonpolar_matched",
                              11: "polar_matched", 12: "acidic_matched", 13: "baic_matched"}, inplace=True)
    unique_df.to_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_Gaps.csv"), index=False)


phenotypes = ["res", "sus"]
drugs = ["FOSF"]

for drug in drugs:
    for phenotype in phenotypes:
        get_mutations(drug, phenotype)
