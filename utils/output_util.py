import os
import utils.dir_utils as dir_utils


class DataOutput:

    def __init__(self, file_name, drug, phenotype):
        self.drug_dirs = dir_utils.DrugDirs(drug, phenotype)
        self.output_file = os.path.join(self.drug_dirs.drug_dir, file_name)

    def add_headers(self, header_string):
        with open(self.output_file, "w+") as output_file:
            output_file.write(header_string)

    def write_recip_info(self, target_organism, gene_name, bitscore, match_ratio, recip_gene_name, recip_organism_name):
        with open(self.output_file, "a+") as output_file:
            output_file.write(f"{target_organism},{gene_name},{recip_organism_name},{recip_gene_name},{str(bitscore)},"
                              f"{match_ratio}\n")

    def write_cell(self, data, is_start=False, is_end=False):
        # write to the file based on the type of data being inputted
        with open(self.output_file, "a+") as output_file:
            if is_start:
                output_file.write(str(data) + ",")
            else:
                output_file.write(str(data) + ",")
            if is_end:
                output_file.write("\n")

    def write_unique_info(self, target_gene, recip_gene, bitscore, match_ratio, qcov):
        # write each cell of the csv
        self.write_cell(target_gene.organism, is_start=True)
        self.write_cell(target_gene.name)
        self.write_cell(recip_gene.organism)
        self.write_cell(recip_gene.name)
        self.write_cell(bitscore)
        self.write_cell(match_ratio)
        self.write_cell(qcov)
        self.write_cell("Remember to fix")
        self.write_cell(recip_gene.type_data)
        self.write_cell(recip_gene.location)
        self.write_cell(recip_gene.start)
        self.write_cell(recip_gene.stop)
        self.write_cell(recip_gene.strand)
        self.write_cell(recip_gene.function_data)
        self.write_cell(recip_gene.aliases)
        self.write_cell(recip_gene.figfam)
        self.write_cell(recip_gene.evidence_codes, is_end=True)
