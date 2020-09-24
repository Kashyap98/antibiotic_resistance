import os

import pandas as pd

from models import Blast as b
from models.Gene import Gene
from utils.dir_utils import OrganismDirs, DrugDirs
from utils import dir_utils
import utils.gen_utils as gen_utils


def get_uniques(drug, phenotype):
    # create output file and write the headers
    drug_dirs = DrugDirs(drug, phenotype)
    all_strain_analysis_dir = os.path.join(drug_dirs.drug_dir, "all_strains")
    op_phenotype = gen_utils.get_op_phenotype(phenotype)
    drug_dirs.set_opposite_phenotype_file(op_phenotype)
    dir_utils.generate_dir(all_strain_analysis_dir)

    # Main Organism
    organism_file = pd.read_csv(drug_dirs.target_phenotype_file, header=None)
    all_orgs = list(organism_file.iloc[:, 0])
    target_organism_dirs = OrganismDirs(all_orgs[0])
    op_organism_file = pd.read_csv(drug_dirs.op_phenotype_file, header=None)
    op_orgs = list(op_organism_file.iloc[:, 0])
    # all_orgs.extend(op_orgs)

    # Get file with all the Genes
    gene_file = pd.read_csv(os.path.join(drug_dirs.drug_dir, f"{phenotype}_sus_aa_no_perfect_matches.csv"))
    res_genes = list(gene_file["gene"])
    res_genes_length = len(res_genes)
    print("Number of all Genes: ", res_genes_length)

    # key = gene, value = list of genes found in order of all_orgs
    complete_gene_dict = {}
    current_gene_progress_count = 0
    for res_gene in res_genes:
        res_gene_object = Gene(target_organism_dirs.organism, res_gene, get_info=True)
        gene_list = []
        current_gene_progress_count += 1

        for organism in all_orgs:
            organism_dirs = OrganismDirs(organism)
            # First blast the first organism gene against the database of the gene in the list.
            blast_data = b.blast(res_gene_object, organism_dirs.database_dir, organism)
            gene_list.append(blast_data.blast_gene)

        print(f"{str(current_gene_progress_count)} / {str(res_genes_length)}")
        complete_gene_dict[res_gene] = gene_list

    with open(os.path.join(all_strain_analysis_dir, f"{phenotype}_aa_all_strains.csv"), "w") as full_strain:
        full_strain.write("organism,gene,blast_organism,blast_gene,qcov,pident,perfect_match,function\n")

    with open(os.path.join(all_strain_analysis_dir, f"{phenotype}_aa_no_perfect.csv"), "w") as no_perfect:
        no_perfect.write("gene,function\n")

    with open(os.path.join(all_strain_analysis_dir, f"{phenotype}_aa_culprit.csv"), "w") as culprit:
        culprit.write("gene,function,perfect_matches")
        for organism in all_orgs:
            culprit.write(f",{organism}")
        culprit.write("\n")

    complete_list_progress = 0
    for res_gene, list_of_genes in complete_gene_dict.items():
        complete_list_progress += 1
        perfect_matches = 0
        res_gene_object = Gene(target_organism_dirs.organism, res_gene, get_info=True)
        culprit_organisms_by_gene = []
        for i in range(len(list_of_genes)):
            culprit_organisms = []
            strain_gene = list_of_genes[i]
            result_strings = []
            for organism in op_orgs:
                organism_dirs = OrganismDirs(organism)
                if strain_gene.organism == organism:
                    continue
                blast_data = b.blast(strain_gene, organism_dirs.database_dir, organism)
                perfect_match = 0
                if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                    perfect_match = 1
                    perfect_matches += 1
                    culprit_organisms.append(organism)
                    result_strings.append(f"{strain_gene.organism},{strain_gene.name},{blast_data.blast_gene.organism},"
                                          f"{blast_data.blast_gene.name},{blast_data.qcov},{blast_data.pident},"
                                          f"{perfect_match},"
                                          f"{strain_gene.function_data}\n")
            with open(os.path.join(all_strain_analysis_dir, f"{phenotype}_aa_all_strains.csv"), "a") as full_strain:
                for string in result_strings:
                    full_strain.write(string)
            culprit_organisms_by_gene.append(culprit_organisms)

        if perfect_matches == 0:
            with open(os.path.join(all_strain_analysis_dir, f"{phenotype}_aa_no_perfect.csv"), "a") as no_perfect:
                no_perfect.write(f"{res_gene_object.name},{res_gene_object.function_data}\n")
        else:
            with open(os.path.join(all_strain_analysis_dir, f"{phenotype}_aa_culprit.csv"), "a") as culprit:
                culprit.write(f"{res_gene_object.name},{res_gene_object.function_data},{perfect_matches}")
                for culprits in culprit_organisms_by_gene:
                    culprit.write(f",{'. '.join(culprits)}")
                culprit.write("\n")
        print(f"{str(complete_list_progress)} / {str(len(complete_gene_dict))}")


PHENOTYPES = ["res"]
DRUGS = ["CIPRO", "SULF"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_uniques(DRUG, PHENOTYPE)
