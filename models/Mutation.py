from Bio import pairwise2


NONPOLAR = ["G", "A", "V", "L", "I", "M", "F", "W", "P"]
POLAR = ["S", "T", "C", "Y", "N", "Q"]
BASIC = ["K", "R", "H"]
ACIDIC = ["D", "E"]


class Mutation:
    def __init__(self, seq_1, seq_2):
        # 2, -1, -11, -1
        alignment = pairwise2.align.localms(seq_1, seq_2, 2, -1, -11, -1, one_alignment_only=True)[0]
        self.seq_1 = alignment[0]
        self.seq_2 = alignment[1]
        self.score = int(alignment[2])

        self.gap_reference = 0
        self.gap_matched = 0

        self.nonpolar_polar = 0
        self.nonpolar_acidic = 0
        self.nonpolar_basic = 0

        self.polar_nonpolar = 0
        self.polar_acidic = 0
        self.polar_basic = 0

        self.acidic_nonpolar = 0
        self.acidic_polar = 0
        self.acidic_basic = 0

        self.basic_nonpolar = 0
        self.basic_polar = 0
        self.basic_acidic = 0

        self.nonpolar_reference = 0
        self.polar_reference = 0
        self.acidic_reference = 0
        self.basic_reference = 0

        self.nonpolar_matched = 0
        self.polar_matched = 0
        self.acidic_matched = 0
        self.basic_matched = 0

        self.percentage = float(round(self.score / len(seq_1), 2))

        self.order_mutations()

    def order_mutations(self):
        for pos in range(0, len(self.seq_1)):
            x = self.seq_1[pos]
            y = self.seq_2[pos]
            if x is not y:
                # print(x)
                # print(y)
                # print("-----------------------")

                if x in POLAR:
                    if y in NONPOLAR:
                        self.polar_nonpolar += 1
                    elif y in ACIDIC:
                        self.polar_acidic += 1
                    elif y in BASIC:
                        self.polar_basic += 1

                if x in NONPOLAR:
                    if y in POLAR:
                        self.nonpolar_polar += 1
                    elif y in ACIDIC:
                        self.nonpolar_acidic += 1
                    elif y in BASIC:
                        self.nonpolar_basic += 1

                if x in BASIC:
                    if y in POLAR:
                        self.basic_polar += 1
                    elif y in NONPOLAR:
                        self.basic_nonpolar += 1
                    elif y in ACIDIC:
                        self.basic_acidic += 1

                if x in ACIDIC:
                    if y in POLAR:
                        self.acidic_polar += 1
                    elif y in NONPOLAR:
                        self.acidic_nonpolar += 1
                    elif y in BASIC:
                        self.acidic_basic += 1

                if x is "-":
                    self.order_gaps(y, "reference")
                    self.gap_reference += 1

                if y is "-":
                    self.order_gaps(x, "matched")
                    self.gap_matched += 1

    def order_gaps(self, oppositeAcid, currentStrand):
        if currentStrand is "reference":
            if oppositeAcid in POLAR:
                self.polar_matched += 1
            if oppositeAcid in NONPOLAR:
                self.nonpolar_matched += 1
            if oppositeAcid in BASIC:
                self.basic_matched += 1
            if oppositeAcid in ACIDIC:
                self.acidic_matched += 1
        else:
            if oppositeAcid in POLAR:
                self.polar_reference += 1
            if oppositeAcid in NONPOLAR:
                self.nonpolar_reference += 1
            if oppositeAcid in BASIC:
                self.basic_reference += 1
            if oppositeAcid in ACIDIC:
                self.acidic_reference += 1

    def get_gaps(self):
        return [self.nonpolar_reference, self.polar_reference, self.acidic_reference, self.basic_reference,
                self.nonpolar_matched, self.polar_matched, self.acidic_matched, self.basic_matched]

    def get_mutations(self):
        return [self.nonpolar_polar, self.nonpolar_acidic, self.nonpolar_basic, self.polar_nonpolar,
                self.polar_acidic, self.polar_basic, self.acidic_nonpolar, self.acidic_polar, self.acidic_basic,
                self.basic_nonpolar, self.basic_polar, self.basic_acidic, self.gap_reference, self.gap_matched]
