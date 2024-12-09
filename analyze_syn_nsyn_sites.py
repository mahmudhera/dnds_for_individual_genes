from itertools import product

# Genetic code dictionary
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def count_synonymous_nonsynonymous_sites(codon):
    """
    Calculate the number of synonymous and nonsynonymous sites for a given codon.

    Args:
        codon (str): A three-letter string representing a codon (e.g., 'ATG').

    Returns:
        tuple: (synonymous_sites, nonsynonymous_sites)
    """
    if len(codon) != 3 or not all(base in 'ATGC' for base in codon):
        raise ValueError("Input must be a valid codon (3 letters, A/T/G/C).")

    synonymous_count = 0
    nonsynonymous_count = 0

    for i in range(3):  # Iterate over each nucleotide position
        for nucleotide in 'ATGC':
            if nucleotide == codon[i]:
                continue  # Skip the original nucleotide
            mutated_codon = list(codon)
            mutated_codon[i] = nucleotide
            mutated_codon = ''.join(mutated_codon)
            
            # Determine the type of substitution
            if GENETIC_CODE[mutated_codon] == GENETIC_CODE[codon]:
                synonymous_count += 1
            else:
                nonsynonymous_count += 1

    # Normalize by 3 positions
    synonymous_sites = synonymous_count / 3.0
    nonsynonymous_sites = nonsynonymous_count / 3.0

    return synonymous_sites, nonsynonymous_sites


if __name__ == '__main__':
    # get all the amino acids
    amino_acids = set(GENETIC_CODE.values())

    # go through all the amino acids
    for amino_acid in amino_acids:
        print('-----')
        print('Amino Acid:', amino_acid)
        print('Codons:')
        codons = [codon for codon, aa in GENETIC_CODE.items() if aa == amino_acid]
        for codon in codons:
            syn_sites, nonsyn_sites = count_synonymous_nonsynonymous_sites(codon)
            print(codon, syn_sites, nonsyn_sites)
            
