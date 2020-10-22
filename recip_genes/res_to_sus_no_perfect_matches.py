import os

import pandas as pd

from models import blast as b
from models.gene import Gene
from utils.output_util import DataOutput
from utils.dir_utils import OrganismDirs, DrugDirs
import utils.gen_utils as gen_utils


def get_uniques(drug, phenotype):
    # create output file and write the headers
    output_file = DataOutput(f"{phenotype}_aa_UniqueBlastInfo.csv", drug, phenotype)
    output_file.add_headers("organism, feature_id, recip_organism, recip_feature_id, bitscore, match_ratio, qcov, "
                            "contig_id, type_data, location, start, stop, strand, function_data,"
                            "aliases, figfam, evidence_codes \n")
    drug_dirs = DrugDirs(drug, phenotype)
    op_phenotype = gen_utils.get_op_phenotype(phenotype)
    drug_dirs.set_opposite_phenotype_file(op_phenotype)

    # Main Organism
    organism_file = pd.read_csv(drug_dirs.target_phenotype_file, header=None)
    all_orgs = list(organism_file.iloc[:, 0])
    target_organism_dirs = OrganismDirs(all_orgs[0], converted_ncbi_data=True)
    op_organism_file = pd.read_csv(drug_dirs.op_phenotype_file, header=None)
    op_orgs = list(op_organism_file.iloc[:, 0])

    target_res_orgs_hit_count = len(all_orgs) - 1
    target_sus_orgs_hit_count = len(op_orgs)

    # Get file with all the Genes
    gene_file = pd.read_csv(os.path.join(os.getcwd(), "..", "sorted_data", drug,
                                         f"{phenotype}_aa_recip_detailed_info.csv"))
    res_genes = list(gene_file["gene"])
    res_genes_hit_counts = list(gene_file["hit_count"])
    res_genes_length = len(res_genes)

    print("Number of all Genes: ", res_genes_length)

    qcov_genes = {}
    count = 0
    for i in range(0, res_genes_length):
        current_gene = res_genes[i]
        current_gene_hit_count = res_genes_hit_counts[i]
        if current_gene_hit_count != target_res_orgs_hit_count:
            continue

        for op_organism in op_orgs:
            # For each other genome we want to compare to, Recip BLAST those genes and only keep the ones are the same
            op_organism = str(op_organism)
            # print(f"{str(count)} / {str(len(op_orgs))}")
            op_organism_dirs = OrganismDirs(op_organism, converted_ncbi_data=True)

            # First blast the first organism gene against the database of the gene in the list.
            res_gene = Gene(target_organism_dirs.organism, current_gene, get_info=False)
            blast_data = b.blast(res_gene, op_organism_dirs.database_dir, op_organism)

            # if there is a hit
            if not blast_data:
                # all_genes.remove(gene)
                continue

            if current_gene in qcov_genes:
                old_data = qcov_genes[current_gene]
                new_count = old_data[0] + 1
                new_qcov = (old_data[1] + blast_data.qcov) / 2
                new_pident = (old_data[2] + blast_data.pident) / 2
                if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                    new_perfect_matches = old_data[3] + 1
                else:
                    new_perfect_matches = old_data[3]
                qcov_genes[current_gene] = [new_count, new_qcov, new_pident, new_perfect_matches, old_data[4],
                                            old_data[5]]
            else:
                hit_count = 0
                if int(blast_data.qcov) == 100 and int(blast_data.pident) == 100:
                    perfect_matches = 1
                else:
                    perfect_matches = 0
                qcov_genes[current_gene] = [hit_count, blast_data.qcov, blast_data.pident, perfect_matches,
                                            blast_data.blast_gene.length, blast_data.blast_gene.description]

        print(f"{i} / {res_genes_length}")

    with open(drug_dirs.unique_genes_file, "w") as unique_genes:
        for gene in qcov_genes.keys():
            unique_genes.write(f"{gene}\n")

    with open(os.path.join(drug_dirs.drug_dir, f"{phenotype}_aa_sus_no_perfect_matches.csv"), "w") as qcov_gene_files:
        qcov_gene_files.write("gene,hit_count,qcov_average,pident_average,perfect_matches,gene_length,function\n")
        for gene, info in qcov_genes.items():
            qcov_gene_files.write(f"{gene},{str(info[0])},{str(info[1])},{str(info[2])},{str(info[3])},"
                                  f"{str(info[4])},{str(info[5])}\n")


PHENOTYPES = ["res"]
DRUGS = ["CIPRO"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_uniques(DRUG, PHENOTYPE)
