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

    def write_unique_info(self, target_gene, recip_gene, bitscore, match_ratio):
        # write each cell of the csv
        self.write_cell(target_gene.organism, is_start=True)
        self.write_cell(target_gene.name)
        self.write_cell(recip_gene.organism)
        self.write_cell(recip_gene.name)
        self.write_cell(bitscore)
        self.write_cell(match_ratio)
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

    def write_mutation_data(self, target_organism, opposite_organism, target_gene, opposite_gene, bitscore,
                            recip_bitscore, target_function, opposite_function, nonpolar_polar, nonpolar_acidic,
                            nonpolar_basic,  polar_nonpolar, polar_acidic, polar_basic, acidic_nonpolar, acidic_polar,
                            acidic_basic,  basic_nonpolar, basic_polar, basic_acidic, gap_reference, gap_matched):
        self.write_cell(target_organism, is_start=True)
        self.write_cell(opposite_organism)
        self.write_cell(target_gene)
        self.write_cell(opposite_gene)
        self.write_cell(bitscore)
        self.write_cell(recip_bitscore)
        self.write_cell(target_function)
        self.write_cell(opposite_function)
        self.write_cell(nonpolar_polar)
        self.write_cell(nonpolar_acidic)
        self.write_cell(nonpolar_basic)
        self.write_cell(polar_nonpolar)
        self.write_cell(polar_acidic)
        self.write_cell(polar_basic)
        self.write_cell(acidic_nonpolar)
        self.write_cell(acidic_polar)
        self.write_cell(acidic_basic)
        self.write_cell(basic_nonpolar)
        self.write_cell(basic_polar)
        self.write_cell(basic_acidic)
        self.write_cell(gap_reference)
        self.write_cell(gap_matched, is_end=True)
