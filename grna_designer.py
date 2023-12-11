"""
Author: Yixin Zhang
Last updated: Dec. 11, 2023

This file contains functions for designing gRNAs for various Cas systems.
"""

'''
1. SpCas9 is targeted to genomic loci matching a 20-nt spacer sequence within the crRNA,
   immediately upstream of a required 5'-NGG PAM. sgRNA scaffold sequence from Hsu NBT
   is used:
   GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT
2. SaCas9 achieves the highest editing efficiency in mammalian cells with spacers
   between 21 and 23 nucleotides long. This algorithm will choose the average 22 as
   spacer length. SaCas9 has a PAM sequence of 5'-NNGRR immediately downstream of the
   protospacer. R stands for a purine (A or G). sgRNA scaffold sequence from Ran Nature
   2015 is used:
   GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTT
3. The PAM for FnCpf1 (FnCas12a) is located upstream of the 5' end of the displaced
   strand of the protospacer and has the sequence 5'-TTN.  FnCpf1 requires a minimum of
   18 nt of spacer sequence to achieve efficient DNA cleavage in vitro. sgRNA scaffold
   sequence from Zetsche Cell 2015 is used: AATTTCTACTGTTGTAGAT; termination: TTTTTT
4. LbCas12a has the PAM TTTV, which is upstream of the protospacer, as the optimal. V
   stands for A, C, or G. 23-nt spacer LbCas12a mediates gene targeting effectively.
   sgRNA scaffold sequence from Vu Front Plant Sci 2021 is used:
   TAATTTCTACTAAGTGTAGAT; termination: TTTTTT
5. LshCas13a has the best efficiency with a spacer length of 20-28 nucleotides long.
   This algorithm will choose the average 24 as spacer length. Unlike Cas9 and Cas12
   orthologs, Cas13 cleaves RNAs and does not require a PAM sequence. Studies have shown
   targeting regions with few secondary structures will improve efficiency. 5' sgRNA
   scaffold sequence from Bandaru Scientific Reports 2020 is used:
   GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTG
6. LwCas13a is similar to LshCas13a, but displays higher efficiency. Gootenberg Science
   2017 uses a 28-nt long spacer. Their 5'scaffold sequence is used:
   GGGGAUUUAGACUACCCCAAAAACGAAGGGGACUAAAAC
'''

import requests
from utils import reverse_complement, hairpin_counter

SpCas9 = {'pam': ['GG'],
          'n_length': 1,
          'pam_length': 2,
          'spacer_length': 20,
          'scaffold': 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT'}
SaCas9 = {'pam': ['GAA', 'GAG', 'GGA', 'GGG'],
          'n_length': 2,
          'pam_length': 3,
          'spacer_length': 22,
          'scaffold': 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTT'}
FnCas12a = {'pam': ['TT'],
            'n_length': 1,
            'pam_length': 3,
            'spacer_length': 18,
            'scaffold': 'AATTTCTACTGTTGTAGAT'}
LbCas12a = {'pam': ['TTTA', 'TTTC', 'TTTG'],
            'n_length': 0,
            'pam_length': 4,
            'spacer_length': 23,
            'scaffold': 'TAATTTCTACTAAGTGTAGAT'}
LshCas13a = {'spacer_length': 24,
            'scaffold': 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTG'}
LwCas13a = {'spacer_length': 28,
            'scaffold': 'GGGGATTTAGACTACCCCAAAAACGAAGGGGACTAAAAC'}
SYSTEM_TO_CAS = {'1': SpCas9, '2': SaCas9, '3': FnCas12a, '4': LbCas12a, '5': LshCas13a, '6': LwCas13a}

def generate_gRNA_from_seq(gene, system, num):
    """
    Generate gRNAs for a gene with scaffold attached.

    Args:
    gene (str): The gene sequence.
    system (str): The chosen Cas system.
    num (int): The number of gRNAs for each gene.

    Returns:
    list: A list of lists, where each sublist contains gRNAs.
    """
    if system == '5' or system == '6':
        return generate_PAMless_gRNA(gene, system, num)
    
    pams = SYSTEM_TO_CAS[system]['pam']
    n_length = SYSTEM_TO_CAS[system]['n_length']
    pam_length = SYSTEM_TO_CAS[system]['pam_length']
    spacer_length = SYSTEM_TO_CAS[system]['spacer_length']
    scaffold = SYSTEM_TO_CAS[system]['scaffold']
    grnas = []
    prev_start = 0
    current_start = 0
    i = 0
    if system == '1' or system == '2':  # Cas9
        while i < len(gene)-pam_length:
            if len(grnas) == num:  # Return if there are enough gRNAs
                return grnas
            potential_pam = gene[i:i+pam_length]
            spacer = None
            if i >= spacer_length+n_length and potential_pam in pams:  # Check for upstream spacer
                spacer = gene[i-n_length-spacer_length:i-n_length]
                current_start = i-n_length-spacer_length
            if not spacer and i <= len(gene)-spacer_length-n_length-pam_length:  # Check for downstream spacer
                if potential_pam in [reverse_complement(pam) for pam in pams]:
                    spacer = gene[i+pam_length+n_length:i+pam_length+n_length+spacer_length]
                    current_start = i+pam_length+n_length
            i += 1
            if not spacer or current_start - prev_start <= spacer_length:
                continue
            grnas.append(spacer + scaffold)
            prev_start = current_start
            i += spacer_length  # Guarantees distinct gRNAs
    else:  # Cas12a
        while i < len(gene)-pam_length:
            if len(grnas) == num:  # Return if there are enough gRNAs
                return grnas
            potential_pam = gene[i:i+pam_length]
            spacer = None
            if i <= len(gene)-spacer_length-n_length-pam_length and potential_pam in pams:  # Check for downstream spacer
                spacer = gene[i+pam_length+n_length:i+pam_length+n_length+spacer_length]
                current_start = i-n_length-spacer_length
            if not spacer and i >= spacer_length+n_length:  # Check for upstream spacer
                if potential_pam in [reverse_complement(pam) for pam in pams]:
                    spacer = gene[i-n_length-spacer_length:i-n_length]
                    current_start = i+pam_length+n_length
            i += 1
            if not spacer or current_start - prev_start <= spacer_length:
                continue
            prev_start = current_start
            grnas.append(scaffold + spacer + 'TTTTTT')
            i += spacer_length  # Guarantees distinct gRNAs
    return grnas
            
def generate_gRNA_from_coordinate(system, chromosome, begin, end, num):
    """
    Generate gRNAs for a gene at the specified location with scaffold attached.

    Args:
    chromosome_genome (str): The genomic sequence where the gene belongs.
    system (str): The chosen Cas system.
    begin (int): The gene's starting location in the chromosome.
    end (int): The gene's ending location in the chromosome.
    num (int): The number of gRNAs for each gene.

    Returns:
    list: A list of lists, where each sublist contains gRNAs.
    """
    data = requests.get(f"https://api.genome.ucsc.edu/getData/sequence?genome=mm10;chrom={chromosome};start={begin};end={end}", timeout=30)
    gene = data.json()['dna'].upper()
    return generate_gRNA_from_seq(gene.upper(), system, num)

def generate_PAMless_gRNA(gene, system, num):
    """
    Generate gRNAs for a gene with scaffold attached for an RNA-cleaving Cas system.

    Args:
    gene (str): The gene sequence.
    system (str): The chosen Cas system.
    num (int): The number of gRNAs for each gene.

    Returns:
    list: A list of lists, where each sublist contains gRNAs.
    """
    spacer_length = SYSTEM_TO_CAS[system]['spacer_length']
    scaffold = SYSTEM_TO_CAS[system]['scaffold']
    spacers = []
    for i in range(len(gene)-spacer_length+1):
        window = gene[i:i+spacer_length]
        g = hairpin_counter(window)
        spacers.append((window, g))
    spacers = sorted(spacers, key=lambda x: x[1])[:num]
    spacers = [reverse_complement(spacer[0]) for spacer in spacers]
    return [scaffold + spacer for spacer in spacers]

def tiled_gRNAs_from_seq(gene, system, spacing):
    """
    Generate tiled gRNAs with specified spacing for a gene with scaffold attached.

    Args:
    gene (str): The gene sequence.
    system (str): The chosen Cas system.
    spacing (int): The tile size in terms of nucleotides between each gRNA

    Returns:
    list: A list of lists, where each sublist contains gRNAs.
    """
    spacer_length = SYSTEM_TO_CAS[system]['spacer_length']
    scaffold = SYSTEM_TO_CAS[system]['scaffold']
    grnas = []
    i = 0
    while i <= len(gene)-spacer_length:
        if system == '1' or system == '2':  # Cas9
            grnas.append(gene[i:i+spacer_length] + scaffold)
        else:  # Cas12a
            grnas.append(scaffold + gene[i:i+spacer_length] + 'TTTTTT')
        i += spacing
    return grnas

def tiled_gRNAs_from_coordinate(system, chromosome, begin, end, spacing):
    """
    Generate tiled gRNAs with specified spacing for a gene at the specified location with scaffold attached.

    Args:
    chromosome_genome (str): The genomic sequence where the gene belongs.
    system (str): The chosen Cas system.
    begin (int): The gene's starting location in the chromosome.
    end (int): The gene's ending location in the chromosome.
    spacing (int): The tile size in terms of nucleotides between each gRNA

    Returns:
    list: A list of lists, where each sublist contains gRNAs.
    """
    data = requests.get(f"https://api.genome.ucsc.edu/getData/sequence?genome=mm10;chrom={chromosome};start={begin};end={end}", timeout=30)
    gene = data.json()['dna'].upper()
    return tiled_gRNAs_from_seq(gene.upper(), system, spacing)
