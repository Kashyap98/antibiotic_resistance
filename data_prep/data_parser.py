import csv
import glob
import os
import shutil

from Bio import SeqIO

import utils.dir_utils as dir_utils
import pandas as pd


def parse_tab_deliminated_info(write_nuc=None, write_prot=None, write_info=None):
    # Used to create the fasta files and gene information files for future use.
    # for each organism file
    for file in glob.glob(os.path.join(dir_utils.ANNOTATION_DIR, '*.txt')):
        name = os.path.basename(file).strip(".txt")
        file_name = os.path.join(dir_utils.ANNOTATION_DIR, os.path.basename(file))

        data = list(csv.reader(open(file_name, 'r'), delimiter='\t'))
        df = pd.DataFrame.from_records(data)
        # 0 = contig_id, 1 = feature_id, 2 = type, 3 = location, 4 = start, 5 = stop, 6 = strand, 7 = function
        # 8 = aliases, 9 = figfam, 10 = evidence_codes, 11 = nucleotide_sequence, 12 = aa_sequence

        # Make a new set of folders for each organism.
        organism_folder = os.path.join(dir_utils.CONVERTED_DATA_DIR, name)
        organism_gene_folder = os.path.join(organism_folder, "genes")
        organism_gene_info_folder = os.path.join(organism_folder, "gene_info")

        dir_utils.generate_dir(organism_folder)
        dir_utils.generate_dir(organism_gene_folder)
        dir_utils.generate_dir(organism_gene_info_folder)

        # for each gene row, make fasta for gene file.
        gene = 1
        while gene < (len(df)):
            # Parse out information for each gene.
            gene_data = df.iloc[gene]
            contig_id = str(gene_data[0]).replace(",", "_").replace(".", "_").replace("|", "_")
            feature_id = str(gene_data[1]).replace(",", "_").replace(".", "_").replace("|", "_")
            type_data = str(gene_data[2])
            location = str(gene_data[3]).replace(",", "_").replace(".", "_").replace("|", "_")
            start = str(gene_data[4])
            stop = str(gene_data[5])
            strand = str(gene_data[6])
            function_data = str(gene_data[7]).replace(",", "_").replace(".", "_").replace("|", "_")
            aliases = str(gene_data[8]).replace(",", "_").replace(".", "_").replace("|", "_")
            figfam = str(gene_data[9])
            evidence_codes = str(gene_data[10]).replace(",", "_").replace(".", "_").replace("|", "_")
            nucleotide_sequence = str(gene_data[11])
            aa_sequence = str(gene_data[12])

            if write_nuc:
                with open(os.path.join(organism_gene_folder, f"{feature_id}.fasta"), "w") as fasta_file:
                    fasta_file.write(">" + feature_id + "\n" + nucleotide_sequence + "\n")
                    fasta_file.write(f"> {feature_id} \n{nucleotide_sequence}\n ")

            if write_prot:
                with open(os.path.join(organism_gene_folder, f"{feature_id}-protein.fasta"), "w") as fasta_file:
                    fasta_file.write(f"> {feature_id} \n{aa_sequence}\n ")

            if write_info:
                # for each gene row, make text document with info.
                with open(os.path.join(organism_gene_info_folder, f"{feature_id}.txt"), "w") as info_file:
                    info_file.write(f"contig_id, {contig_id}\n")
                    info_file.write(f"feature_id, {feature_id}\n")
                    info_file.write(f"type_data, {type_data}\n")
                    info_file.write(f"location, {location}\n")
                    info_file.write(f"start, {start}\n")
                    info_file.write(f"stop, {stop}\n")
                    info_file.write(f"strand, {strand}\n")
                    info_file.write(f"function_data, {function_data}\n")
                    info_file.write(f"aliases, {aliases}\n")
                    info_file.write(f"figfam, {figfam}\n")
                    info_file.write(f"evidence_codes, {evidence_codes}\n")
                    info_file.write(f"nucleotide_sequence, {nucleotide_sequence}\n")
                    info_file.write(f"aa_sequence, {aa_sequence}\n")

            gene += 1


def parse_fasta_file_data():
    organism_labels = pd.read_csv(dir_utils.ORGANISM_LABELS_FILE, header=0)
    CONVERSIONS = dict(zip(list(organism_labels.iloc[:, 2]), list(organism_labels.iloc[:, 0])))
    dir_utils.generate_dir(dir_utils.CONVERTED_DATA_DIR)

    for org_file in os.listdir(dir_utils.NCBI_DATA_DIR):
        org_file_path = os.path.join(dir_utils.NCBI_DATA_DIR, org_file)
        ncbi_label = org_file.split(".")[0]
        if ncbi_label == "converted_ncbi_data":
            continue

        print(CONVERSIONS[ncbi_label])
        if "UMB" in CONVERSIONS[ncbi_label]:
            converted_name = CONVERSIONS[ncbi_label]
        else:
            converted_name = f"CP-{CONVERSIONS[ncbi_label]}"

        converted_org_path = os.path.join(dir_utils.CONVERTED_DATA_DIR, converted_name)
        gene_folder_path = os.path.join(converted_org_path, "genes")
        dir_utils.generate_dir(converted_org_path)
        dir_utils.generate_dir(gene_folder_path)

        os.chdir(dir_utils.NCBI_DATA_DIR)
        shutil.copy(org_file_path, os.path.join(converted_org_path, f"{converted_name}.faa"))
        os.chdir(converted_org_path)
        os.system(f'makeblastdb -in {os.path.join(converted_org_path, f"{converted_name}.faa")} '
                  f'-title {converted_name} -out {converted_name} -dbtype prot')

        for record in SeqIO.parse(org_file_path, "fasta"):
            with open(os.path.join(gene_folder_path, f"{record.id}.fasta"), "w") as new_record_file:
                new_record_file.write(f"> {record.description.replace(',', '_')}\n")
                new_record_file.write(f"{record.seq}")
