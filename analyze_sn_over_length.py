from main import count_total_syn_nsyn_sites
import random
from Bio import SeqIO
from matplotlib import pyplot as plt

if __name__ == '__main__':
    length = 3**12
    genome_filename = 'ecoli_k12.fna'
    
    # Read the genome sequence
    dna_seq = ''
    for record in SeqIO.parse(genome_filename, 'fasta'):
        dna_seq += str(record.seq)

    s_n_ratios = []  
    for i in range(4,14):
        syn_sites, nonsyn_sites = count_total_syn_nsyn_sites(dna_seq[:3**i])
        s_n_ratios.append(syn_sites/nonsyn_sites)

    plt.plot(range(4,14), s_n_ratios)
    plt.xlabel('Length of DNA sequence')
    plt.ylabel('Ratio of synonymous to nonsynonymous sites')
    plt.title('$S/N$ ratio over DNA sequence length (E. coli K12)')
    # set x tick labels
    plt.xticks(range(4,14), ['$3^{%s}$' % str(i) for i in range(4,14)])
    plt.savefig('syn_nsyn_sites_ratio.pdf')


    genes_filename = "genes.fasta"
    for record in SeqIO.parse("genes.fasta", "fasta"):
        if 'nt_seq' in record.id:
            nt_seq = str(record.seq)
            syn_sites, nonsyn_sites = count_total_syn_nsyn_sites(nt_seq)
            print('Gene:', record.id)
            print('Synonymous sites:', syn_sites)
            print('Nonsynonymous sites:', nonsyn_sites)
            print('Ratio:', syn_sites/nonsyn_sites)
