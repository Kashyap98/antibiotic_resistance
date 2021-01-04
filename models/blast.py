from Bio.Blast.Applications import NcbiblastpCommandline


# organize the blast result in an easily accessible class
from models.gene import Gene
from utils import dir_utils


class BlastResult:

    def __init__(self, blast_data, blasted_gene, blast_organism, bitscore_threshold=1000, qcov_threshold=80):
        self.gene_name = blast_data[0]
        self.bitscore_threshold = bitscore_threshold
        self.qcov_threshold = qcov_threshold
        self.target_gene = blasted_gene
        self.blast_gene = Gene(blast_organism, self.gene_name)

        self.match_length = int(blast_data[1])
        self.bitscore = int(float(bitscore_fixer(blast_data, blast_data[2])))
        self.qcov = int(bitscore_fixer(blast_data, blast_data[3]))
        self.pident = float(bitscore_fixer(blast_data, blast_data[4]))
        self.match_ratio = float(round(int(self.match_length) / int(self.target_gene.length), 2))

        self.is_perfect_match = bool(int(self.qcov == 100) and int(self.pident == 100))
        self.is_homolog = bool(self.above_bitscore_threshold_check() and self.within_size_constraints_check()
                               and self.above_qcov_threshold_check())

    def above_bitscore_threshold_check(self):
        if self.bitscore >= self.bitscore_threshold:
            return True
        return False

    def within_size_constraints_check(self):
        if self.blast_gene.upper_length >= self.target_gene.length >= self.blast_gene.lower_length:
            return True
        return False

    def above_qcov_threshold_check(self):
        if self.qcov >= self.qcov_threshold:
            return True
        return False


class CombinedResult:

    def __init__(self, gene_name: str, gene_description: str):
        self.gene_name = gene_name
        self.gene_description = gene_description
        self.results = 0

        self.bitscore_average = None
        self.qcov_average = None
        self.pident_average = None

        self.match_count = 0
        self.perfect_match_count = 0

        self.match_organisms = []
        self.perfect_match_organisms = []

    def add_new_result(self, result: BlastResult):
        if self.results == 0:
            self.bitscore_average = result.bitscore
            self.qcov_average = result.qcov
            self.pident_average = result.pident
        else:
            self.bitscore_average = round((self.bitscore_average + result.bitscore) / 2, 2)
            self.qcov_average = round((self.qcov_average + result.qcov) / 2, 2)
            self.pident_average = round((self.pident_average + result.pident) / 2, 2)
            
        self.results += 1
        if result.is_perfect_match:
            self.perfect_match_count += 1
            self.perfect_match_organisms.append(result.blast_gene.organism)
        else:
            self.match_count += 1
            self.match_organisms.append(result.blast_gene.organism)
            
    def header(self):
        return ["gene_name", "gene_description", "match_count", "perfect_match_count", "bitscore_average",
                "qcov_average", "pident_average", "match_organisms", "perfect_match_organisms"]
            
    def data(self):
        return [self.gene_name, self.gene_description, self.match_count, self.perfect_match_count,
                self.bitscore_average, self.qcov_average, self.pident_average, self.match_organisms,
                self.perfect_match_organisms]
            

def bitscore_fixer(result, bitscore):
    if len(result) > 3:
        bitscore = bitscore.split("\n")[0]
    else:
        bitscore = bitscore.strip("\n")

    return bitscore


# Perform the blast
def blast(gene_to_blast: Gene, blast_organism: str):
    blast_organism_dirs = dir_utils.OrganismDirs(blast_organism)
    blast = NcbiblastpCommandline(query=gene_to_blast.fasta_file, db=blast_organism_dirs.database_dir,
                                  outfmt='"10 sseqid length bitscore qcovs pident"', max_target_seqs=1)
    result = list(blast())
    result_list = result[0].split(",")

    if len(result_list) > 1:
        blast_result = BlastResult(result_list, gene_to_blast, blast_organism_dirs.organism)
    else:
        blast_result = None

    return blast_result


