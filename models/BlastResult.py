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
        self.match_length = int(blast_data[1])
        self.bitscore = int(float(bitscore_fixer(blast_data, blast_data[2])))
        self.target_gene = target_gene
        self.blast_gene = Gene(blast_organism, self.gene_name, get_info=True)
        self.blast_in_threshold = self.threshold_check()

    def threshold_check(self):
        if self.bitscore >= 1000 and self.blast_gene.upper_length <= self.target_gene.length and self.blast_gene.lower_length <= self.target_gene.length:
            return True
        else:
            return False




