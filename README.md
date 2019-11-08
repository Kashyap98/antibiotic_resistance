Fellowship Project

-Files & Explanations.

-Data Preparation (in order).

data_parser - Used to create the fasta files and the gene information files based on the gene annotation files that were given.

organism_sorter - Used to sort through each organism and placing them in the correct phenotype file so that they may be compared later.

database_gen - Used to compile all the fasta files generated above into 1 fasta file. This file is then made into a Blastn database for each organism.

-BLASTing (in order).

GeneMatcher - Blasts all the genes of 1 organism against all the organisms of the same phenotype. This gives us a list of all recip genes.

uniqueChecker - This blasts all the recip genes from above against all the organisms of the opposite phenotypes. If there is any match the gene is removed as it is not unique to the phenotype.

commonGenes - This is used to create a csv with the information of the genes that are all found in the same phenotype but are not found in any organisms of the opposite phenotype.

-Misc.

Blast - Used to return the result of a BLAST against a target organism using a specific gene.

BlastResult - Contains commonly used information for each BLAST result

Gene - Parsed out gene information into a class

Mutation - Generates the list of mutations and gaps between each gene in a pairwise fashion

