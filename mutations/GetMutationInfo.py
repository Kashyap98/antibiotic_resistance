import os

import pandas as pd

from models.Gene import Gene


def get_mutations(drug, phenotype):

    # Main Organism
    dataFile = pd.read_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_Mutations.csv"))
    originalOrganism = dataFile["Original Organism"][0]
    print(dataFile["Original Organism"][0])
    dataFile.insert(loc=len(dataFile.columns), column="Original Gene Function", value=None)
    # dataFile.insert(loc=len(dataFile.columns), column="Matched Gene Function", value=None)

    for index, row in dataFile.iterrows():

        originalGene = Gene(row["Original Organism"], row["Original Gene"].replace(".fasta", ""), get_info=True)
        # matchedGene = Gene(row["Matched Organism"], row["Matched Gene"], get_info=True)

        dataFile.set_value(index, "Original Gene Function", originalGene.function_data)
        # dataFile.set_value(index, "Matched Gene Function", matchedGene.function_data)

    dataFile.to_csv(os.path.join(os.getcwd(), "sorted_data", drug, phenotype + "_InfoMutations.csv"), index=False)


phenotypes = ["res"]
drugs = ["CIPRO"]

for drug in drugs:
    for phenotype in phenotypes:
        get_mutations(drug, phenotype)
