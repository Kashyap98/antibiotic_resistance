import os

from models import blast as b
from models.gene import Gene
from utils.dir_utils import OrganismDirs, DrugDirs
from utils import dir_utils, output_util
import utils.gen_utils as gen_utils


def get_homologs(drug, phenotype):
    # create output file and write the headers
    drug_dirs = DrugDirs(drug, phenotype)
    dir_utils.generate_dir(drug_dirs.homolog_dir)

    res_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.res_file)
    sus_organisms = gen_utils.get_organisms_by_phenotype(drug_dirs.sus_file)

    first_organism = res_organisms.pop(0)

    gene_labels = ["hydrofolate", "hydropteroate"]
    genes = ["WP_001389366.1.fasta", "WP_001043265.1.fasta"]
    final_genes_dict = {}
    gene_objects = []

    for gene in genes:
        new_gene = Gene(first_organism, gene, get_info=False)
        gene_objects.append(new_gene)
        final_genes_dict[gene] = [new_gene]

    count = 0
    for organism in res_organisms:
        count += 1
        print(f"count {count}")
        organism_dirs = OrganismDirs(organism, converted_ncbi_data=True)

        for i in range(len(gene_objects)):
            gene_name = genes[i]
            compare_gene = gene_objects[i]

            blast_data = b.blast(compare_gene, organism_dirs.database_dir, organism)

            for sus_organism in sus_organisms:

                sus_organism_dir = OrganismDirs(sus_organism, converted_ncbi_data=True)

                sus_blast_data = b.blast(blast_data.blast_gene, sus_organism_dir.database_dir, sus_organism)

                if sus_blast_data:
                    if sus_blast_data.qcov == 100 and sus_blast_data.pident == 100:
                        blast_data.blast_gene.match_list.append(sus_blast_data.blast_gene.organism)

            temp_list = final_genes_dict[gene_name]
            temp_list.append(blast_data.blast_gene)
            final_genes_dict[gene_name] = temp_list

    print(len(final_genes_dict))

    for i in range(len(gene_objects)):
        gene_label = gene_labels[i]
        gene_name = genes[i]

        gene_objects_for_gene_name = final_genes_dict[gene_name]

        output_file = output_util.OutputFile(os.path.join(drug_dirs.homolog_dir, f"{drug}_{gene_label}_homologs.csv"),
                                             ["organism", "gene_description", "perfect_matches_count",
                                              "perfect_matches", "sequence_length", "sequence"])

        for gene_info in gene_objects_for_gene_name:
            output_file.write_data_list_to_output_file([gene_info.organism, gene_info.description,
                                                        len(gene_info.match_list), gene_info.match_list,
                                                        len(gene_info.nuc_sequence), gene_info.nuc_sequence])


PHENOTYPES = ["res"]
DRUGS = ["SULF"]

for DRUG in DRUGS:
    for PHENOTYPE in PHENOTYPES:
        get_homologs(DRUG, PHENOTYPE)
