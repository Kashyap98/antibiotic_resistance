Fellowship Project

-Files & Explanations.

-Data Preparation (in order).

data_parser - Used to create the fasta files and the gene information files based on the gene annotation files that were given.

organism_sorter - Used to sort through each organism and placing them in the correct phenotype file so that they may be compared later.

database_gen - Used to compile all the fasta files generated above into 1 fasta file. This file is then made into a Blastn database for each organism.

-BLASTing (in order).

get_recip_genes - Blasts all the genes of 1 organism against all the organisms of the same phenotype. This gives us a list of all recip genes.

unique_checker - This blasts all the recip genes from above against all the organisms of the opposite phenotypes. If there is any match the gene is removed as it is not unique to the phenotype.

- Mutations
compare_recip_mutations - Compare all the mutations between two genes found through reciprocal blast.

filter_mutations - Remove uninteresting mutation counts

get_mutation_info - Find all the information associated with a gene with mutations

get_mutations - Needs refactoring, gets the mutations in a gene compared to it's reciprocal.

mutations_cleanser - Similar to filter_mutations

- Variations

find_variations - Find variations in trains in opposite phenotypes

get_trees - Generate trees and msa from each set of genes found in all genes of both organisms

variations_filter - Filter all the genes into separate groups for analysis

-Misc.

Blast - Used to return the result of a BLAST against a target organism using a specific gene.

BlastResult - Contains commonly used information for each BLAST result

Gene - Parsed out gene information into a class

Mutation - Generates the list of mutations and gaps between each gene in a pairwise fashion

