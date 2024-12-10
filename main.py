from Bio import SeqIO
from Bio.Seq import Seq
import random
from matplotlib import pyplot as plt
from itertools import product
import functools


from calculate_kaks_by_kakscalculator2 import compute_kaks_ng



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

@functools.cache
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


def count_total_syn_nsyn_sites(nt_seq):
    if len(nt_seq) % 3 != 0:
        raise ValueError("Input must be a valid nucleotide sequence (multiple of 3 letters).")
    
    total_syn_sites = 0
    total_nsyn_sites = 0
    for i in range(0, len(nt_seq), 3):
        codon = nt_seq[i:i+3]
        if codon not in GENETIC_CODE:
            raise ValueError(f"Invalid codon found: {codon}")
        
        syn_sites, nsyn_sites = count_synonymous_nonsynonymous_sites(codon)
        total_syn_sites += syn_sites
        total_nsyn_sites += nsyn_sites

    return total_syn_sites, total_nsyn_sites



"""
This fn will randomly mutate the input string with a given mutation rate.
It will not touch the last stop codon in the sequence.
It will also make sure that there are no stop codons in the mutated sequence.
Resulting sequence may have a higher mutation rate than intended due to the above constraints.
"""
def mutate_string(in_str, mut_rate):
    
    out_str = list(in_str)
    other_bases = {
        'A': ['T', 'C', 'G'],
        'C': ['T', 'A', 'G'],
        'G': ['T', 'C', 'A'],
        'T': ['A', 'C', 'G']
    }
    
    for i in range(len(out_str)-3):
        if random.uniform(0,1) < mut_rate:
            out_str[i] = random.choice( other_bases[out_str[i]] )

    return_str = ''.join(out_str)
    stop_codons = ['TAA', 'TAG', 'TGA']
    for stop_codon in stop_codons:
        # find all occurences of stop codon
        stop_codon_indices = [i for i in range(len(return_str)-2) if return_str[i:i+3] == stop_codon and i%3==0]
        
        # if stop codon found, replace with random codon
        for stop_codon_index in stop_codon_indices:
            if stop_codon_index == len(return_str)-3:
                continue
            while True:
                random_codon = ''.join([random.choice(['A', 'C', 'G', 'T']) for _ in range(3)])
                if random_codon not in stop_codons:
                    break
            return_str = return_str[:stop_codon_index] + random_codon + return_str[stop_codon_index+3:]


    return return_str



def get_kmer_set(seq, k):
    kmer_set = set()
    for i in range(len(seq)-k+1):
        kmer_set.add(seq[i:i+k])
    return kmer_set



def get_containment_index_from_kmer_sets(kmer_set1, kmer_set2):
    return len(kmer_set1.intersection(kmer_set2))/len(kmer_set1)


def get_mutation_rate_from_containment(containment_index, k):
    return 1.0 - (containment_index)**(1/k)


def count_syn_nsyn_mutations(dna_seq1, dna_seq2, aa_seq1, aa_seq2):
    num_total_mutations = 0
    num_syn_mutations = 0
    num_non_syn_mutations = 0
    for i in range(len(dna_seq1)):
        if dna_seq1[i] != dna_seq2[i]:
            num_total_mutations += 1
            codon_index = i//3
            if aa_seq1[codon_index] == aa_seq2[codon_index]:
                num_syn_mutations += 1
            else:
                num_non_syn_mutations += 1

    return num_syn_mutations, num_non_syn_mutations


def determine_correct_dnds(dna_seq1, dna_seq2, aa_seq1, aa_seq2):
    num_syn_mutations, num_non_syn_mutations = count_syn_nsyn_mutations(dna_seq1, dna_seq2, aa_seq1, aa_seq2)
    num_syn_sites, num_non_syn_sites = count_total_syn_nsyn_sites(dna_seq1)
    try:
        ret_val = (1.0*num_non_syn_mutations/num_non_syn_sites) / (1.0*num_syn_mutations/num_syn_sites)
    except ZeroDivisionError:
        ret_val = 200

    return ret_val




def simulate_and_plot(gene_name, nt_seq, aa_seq, ksize_nt, ksize_aa, mutation_rate, num_simulations):
    
    orig_nt_seq = str(nt_seq)
    orig_aa_seq = str(aa_seq)

    orig_nt_kmer_set = get_kmer_set(orig_nt_seq, ksize_nt)
    orig_aa_kmer_set = get_kmer_set(orig_aa_seq, ksize_aa)

    dnds_est_list = []
    dnds_correct_list = []
    dnds_using_tool_list = []

    for i in range(num_simulations):
        mutated_nt_seq = Seq(mutate_string(str(nt_seq), mutation_rate))
        mutated_aa_seq = str(mutated_nt_seq.translate())[:-1]
        mutated_nt_seq = str(mutated_nt_seq)

        mutated_nt_kmer_set = get_kmer_set(mutated_nt_seq, ksize_nt)
        mutated_aa_kmer_set = get_kmer_set(mutated_aa_seq, ksize_aa)

        nt_containment_index = get_containment_index_from_kmer_sets(mutated_nt_kmer_set, orig_nt_kmer_set)
        aa_containment_index = get_containment_index_from_kmer_sets(mutated_aa_kmer_set, orig_aa_kmer_set)

        nt_mutation_rate = get_mutation_rate_from_containment(nt_containment_index, ksize_nt)
        aa_mutation_rate = get_mutation_rate_from_containment(aa_containment_index, ksize_aa)

        try:
            dnds_est = aa_mutation_rate/( 1.0 - aa_mutation_rate - (1-nt_mutation_rate)**3 )
        except ZeroDivisionError:
            dnds_est = 200

        num_syn_sites, num_nsyn_sites = count_total_syn_nsyn_sites(orig_nt_seq)
        dnds_est = dnds_est * num_syn_sites / num_nsyn_sites
        
        dnds_correct = determine_correct_dnds(orig_nt_seq, mutated_nt_seq, orig_aa_seq, mutated_aa_seq)
        dnds_using_tool = compute_kaks_ng(orig_nt_seq, mutated_nt_seq)[2]
        
        dnds_est_list.append(dnds_est)
        dnds_correct_list.append(dnds_correct)
        dnds_using_tool_list.append(dnds_using_tool)

    low, high = 0, 5

    plt.cla()
    plt.xlim(low, high)
    plt.ylim(low, high)
    plt.plot([low, high], [low, high], color='red', alpha=0.5, linestyle='--')
    plt.scatter(dnds_using_tool_list, dnds_est_list, alpha=0.5)
    plt.xlabel('Correct dN/dS by KaKs_calculator')
    #plt.ylabel('Correct dN/dS by my calculation')
    plt.ylabel('Estimated dN/dS by our estimator')
    # title using gene name, ksize, and mutation rate
    plt.title(gene_name + ' ksize_nt=' + str(ksize_nt) + ' mutation_rate=' + str(mutation_rate))
    
    # plot a line y=x
    
    plt.legend()
    plt.savefig(gene_name + '_ksize_nt=' + str(ksize_nt) + '_mutation_rate=' + str(mutation_rate) + '.pdf')




def main():

    # set random seed
    random.seed(0)

    genes_filename = "genes.fasta"
    num_simulations = 1000

    genes_to_nt_seqs =  {}
    genes_to_aa_seqs = {}

    for record in SeqIO.parse(genes_filename, "fasta"):
        gene_name = record.id.split("|")[0]
        seq_type = record.id.split("|")[1]
        if seq_type == "nt_seq":
            genes_to_nt_seqs[gene_name] = record.seq
        elif seq_type == "aa_seq":
            genes_to_aa_seqs[gene_name] = record.seq

    for gene_name, nt_seq in genes_to_nt_seqs.items():
        print('Working on gene:', gene_name)
        print(count_total_syn_nsyn_sites(nt_seq))
        for ksize_nt in [39, 45]:
            for mutation_rate in [0.005, 0.01]:
                simulate_and_plot(gene_name, nt_seq, genes_to_aa_seqs[gene_name], ksize_nt, ksize_nt//3, mutation_rate, num_simulations)
                print('Done with gene:', gene_name, 'ksize:', ksize_nt, 'mutation_rate:', mutation_rate)
        

            

if __name__ == "__main__":
    main()