import os
import utils.dir_utils as dir_utils


class DataOutput:

    def __init__(self, file_name, drug, phenotype):
        self.drug_dirs = dir_utils.DrugDirs(drug, phenotype)
        self.output_file = os.path.join(self.drug_dirs.drug_dir, file_name)

    def add_headers(self, header_string):
        with open(self.output_file, "w+") as output_file:
            output_file.write(header_string)

    def write_recip_info(self, target_organism, gene_name, bitscore, match_ratio, recip_gene_name, recip_organism_name,
                         outputFile):
        with open(self.output_file, "a+") as output_file:
            output_file.write(target_organism + "," + gene_name + "," + recip_organism_name + "," + recip_gene_name + ","
                             + str(bitscore) + "," + str(match_ratio))
            output_file.write("\n")
