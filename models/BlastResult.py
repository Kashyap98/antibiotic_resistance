from models.Gene import Gene


def bitscore_fixer(result, bitscore):
    if len(result) > 3:
        bitscore = bitscore.split("\n")[0]
    else:
        bitscore = bitscore.strip("\n")

    return bitscore


# organize the blast result in an easily accessible class
class BlastResult:

    def __init__(self, blast_data, target_gene, blast_organism):
        self.gene_name = blast_data[0]
        self.bitscore_threshold = 1000
        self.qcov_threshold = 90
        self.match_length = int(blast_data[1])
        self.bitscore = int(float(bitscore_fixer(blast_data, blast_data[2])))
        self.qcov = int(bitscore_fixer(blast_data, blast_data[3]))
        self.pident = float(bitscore_fixer(blast_data, blast_data[4]))
        self.target_gene = target_gene
        self.blast_gene = Gene(blast_organism, self.gene_name)
        self.blast_in_threshold = self.threshold_check()
        self.qcov_above_threshold = self.qcov_check()
        self.match_ratio = float(round(int(self.match_length) / int(self.target_gene.length), 2))

    def threshold_check(self):
        if self.bitscore >= self.bitscore_threshold and \
                self.blast_gene.upper_length >= self.target_gene.length >= self.blast_gene.lower_length:
            return True
        else:
            return False

    def qcov_check(self):
        if self.qcov > self.qcov_threshold:
            return True
        else:
            return False




