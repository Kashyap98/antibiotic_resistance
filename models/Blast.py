from Bio.Blast.Applications import NcbiblastnCommandline, NcbiblastpCommandline
from models.BlastResult import BlastResult


# Perform the blast
def blast(target_gene, database_path, blast_organism):
    blastn = NcbiblastpCommandline(query=target_gene.fasta_file, db=database_path,
                                   outfmt='"10 sseqid length bitscore qcovs pident"', max_target_seqs=1)
    result = list(blastn())
    result_list = result[0].split(",")

    if len(result_list) > 1:
        blast_result = BlastResult(result_list, target_gene, blast_organism)
    else:
        blast_result = None

    return blast_result
