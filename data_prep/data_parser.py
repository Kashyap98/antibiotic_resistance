import csv
import glob
import os
import utils.dir_utils as dir_utils
import pandas as pd

# Used to create the fasta files and gene information files for future use.
# for each organism file
for file in glob.glob(os.path.join(dir_utils.ANNOTATION_DIR, '*.txt')):
    name = os.path.basename(file).strip(".txt")
    fileName = os.path.join(dir_utils.ANNOTATION_DIR, os.path.basename(file))

    data = list(csv.reader(open(fileName, 'r'), delimiter='\t'))
    df = pd.DataFrame.from_records(data)
    # 0 = contig_id, 1 = feature_id, 2 = type, 3 = location, 4 = start, 5 = stop, 6 = strand, 7 = function
    # 8 = aliases, 9 = figfam, 10 = evidence_codes, 11 = nucleotide_sequence, 12 = aa_sequence

    # Make a new set of folders for each organism.
    organismFolder = os.path.join(dir_utils.CONVERTED_DATA_DIR, name)
    organismGeneFolder = os.path.join(organismFolder, "genes")
    organismGeneInfoFolder = os.path.join(organismFolder, "gene_info")
    
    dir_utils.generate_dir(organismFolder)
    dir_utils.generate_dir(organismGeneFolder)
    dir_utils.generate_dir(organismGeneInfoFolder)

    # for each gene row, make fasta for gene file.
    gene = 1
    while gene < (len(df)):
        # Parse out information for each gene.
        geneData = df.iloc[gene]
        contig_id = str(geneData[0]).replace(",", "_").replace(".", "_").replace("|", "_")
        feature_id = str(geneData[1]).replace(",", "_").replace(".", "_").replace("|", "_")
        type_data = str(geneData[2])
        location = str(geneData[3]).replace(",", "_").replace(".", "_").replace("|", "_")
        start = str(geneData[4])
        stop = str(geneData[5])
        strand = str(geneData[6])
        function_data = str(geneData[7]).replace(",", "_").replace(".", "_").replace("|", "_")
        aliases = str(geneData[8]).replace(",", "_").replace(".", "_").replace("|", "_")
        figfam = str(geneData[9])
        evidence_codes = str(geneData[10]).replace(",", "_").replace(".", "_").replace("|", "_")
        nucleotide_sequence = str(geneData[11])
        aa_sequence = str(geneData[12])

        # Uncomment to create fasta files at the same time.
        with open(os.path.join(organismGeneFolder, feature_id + ".fasta"), "w") as fastaFile:
            fastaFile.write(">" + feature_id + "\n" + nucleotide_sequence + "\n")

        # with open(os.path.join(organismGeneFolder, feature_id + "-protein.fasta"), "w") as fastaFile:
        #     fastaFile.write(">" + feature_id + "\n" + aa_sequence + "\n")

        # # for each gene row, make text document with info.
        # with open(os.path.join(organismGeneInfoFolder, feature_id + ".txt"), "w") as infoFile:
        #     infoFile.write(
        #         "contig_id, " + contig_id + "\n" + "feature_id, " + feature_id + "\n" + "type_data, " +
        #         type_data + "\n" + "location, " + location + "\n" + "start, " + start + "\n" + "stop, " + stop + "\n" +
        #         "strand, " + strand + "\n" + "function_data, " + function_data + "\n" + "aliases, " + aliases + "\n" +
        #         "figfam, " + figfam + "\n" + "evidence_codes, " + evidence_codes + "\n"
        #         + "nucleotide_sequence, " + nucleotide_sequence + "\n" + "aa_sequence, " + aa_sequence)
        print(f"{gene}")
        gene += 1
