from Bio.Blast.Applications import NcbiblastpCommandline


# organize the blast result in an easily accessible class
from models.gene import Gene


class BlastResult:

    def __init__(self, blast_data, target_gene, blast_organism, bitscore_threshold=1000, qcov_threshold=90):
        self.gene_name = blast_data[0]
        self.bitscore_threshold = bitscore_threshold
        self.qcov_threshold = qcov_threshold
        self.target_gene = target_gene
        self.blast_gene = Gene(blast_organism, self.gene_name)

        self.match_length = int(blast_data[1])
        self.bitscore = int(float(bitscore_fixer(blast_data, blast_data[2])))
        self.qcov = int(bitscore_fixer(blast_data, blast_data[3]))
        self.pident = float(bitscore_fixer(blast_data, blast_data[4]))
        self.match_ratio = float(round(int(self.match_length) / int(self.target_gene.length), 2))

        self.blast_in_threshold = self.above_bitscore_threshold_check()
        self.qcov_above_threshold = self.above_qcov_threshold_check()
        self.within_size_threshold = self.within_size_threshold()

    def above_bitscore_threshold_check(self):
        if self.bitscore >= self.bitscore_threshold:
            return True
        return False

    def within_size_constraints_check(self):
        if self.blast_gene.upper_length >= self.target_gene.length >= self.blast_gene.lower_length:
            return True
        return False

    def above_qcov_threshold_check(self):
        if self.qcov > self.qcov_threshold:
            return True
        return False


def bitscore_fixer(result, bitscore):
    if len(result) > 3:
        bitscore = bitscore.split("\n")[0]
    else:
        bitscore = bitscore.strip("\n")

    return bitscore


# Perform the blast
def blast(target_gene, blast_organism_database_path, blast_organism):
    blast = NcbiblastpCommandline(query=target_gene.fasta_file, db=blast_organism_database_path,
                                  outfmt='"10 sseqid length bitscore qcovs pident"', max_target_seqs=1)
    result = list(blast())
    result_list = result[0].split(",")

    if len(result_list) > 1:
        blast_result = BlastResult(result_list, target_gene, blast_organism)
    else:
        blast_result = None

    return blast_result


