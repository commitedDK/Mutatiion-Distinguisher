# Mutatiion-Distinguisher
Distinguishes between asynonymous and synonymous mutations

The code takes 3 input files -a) genome fasta file, b) gff file, c) vcf file. 
Using information from vcf file about the position of mutations and the altered nucleotide, the mutated genome is created. The mutations which occur only in CDS are considered. The code saves sequence of every CDS into value of a dictionary with position and phase in key and do the same for mutated sequence. Next step checks every codon with reference codon if it cause any amino acid change. If change in codon result in new amino acid then it saves the mutation as non-synonymous and if the there is no change in amino acid then it notes into synonymous. Finally the code extract the position of the mutation in the genome.
