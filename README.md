# dnds_for_individual_genes
We will do the following here:

1. Get about ten long genes, gene names, the aa sequence, and the nt sequence
2. Use simple mutation model to mutate the nt sequence, but making sure that the resulting codon is not a stop codon
3. Get the mutated aa sequence as well
4. Get the kmers from the original and the mutated nt and aa sequences
5. Using the kmers, estimate dnds using our estimator
6. Using the entire sequences, get the true dnds
7. Record the true and the estimated values, and do this for many simulation runs, and for many mutation rates
8. Get the results for many individual genes as well
9. Plot



## Downloaded from
1. FOXP2 downloaded from https://rest.uniprot.org/uniprotkb/O15409.fasta
2. FOXP2 nucleotide sequence obtained from https://useast.ensembl.org