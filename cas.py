"""
Author: Yixin Zhang
Last updated: Dec. 11, 2023

This file contains Cas system definitions.
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

from dataclasses import dataclass

@dataclass(frozen=True)
class SpCas9:
    name = 'SpCas9'
    pams = ['GG']
    pam_length = len(pams[0])
    n_length = 1
    spacer_length = 20
    scaffold = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT'

@dataclass(frozen=True)
class SaCas9:
    name = 'SaCas9'
    pams = ['GAA', 'GAG', 'GGA', 'GGG']
    pam_length = len(pams[0])
    n_length = 2
    spacer_length = 22
    scaffold = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTT'

@dataclass(frozen=True)
class FnCas12a:
    name = 'FnCas12a'
    pams = ['TT']
    pam_length = len(pams[0])
    n_length = 1
    spacer_length = 18
    scaffold = 'AATTTCTACTGTTGTAGAT'

@dataclass(frozen=True)
class LbCas12a:
    name = 'LbCas12a'
    pams = ['TTTA', 'TTTC', 'TTTG']
    pam_length = len(pams[0])
    n_length = 0
    spacer_length = 23
    scaffold = 'TAATTTCTACTAAGTGTAGAT'

@dataclass(frozen=True)
class LshCas13a:
    name = 'LshCas13a'
    spacer_length = 24
    scaffold = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTG'

@dataclass(frozen=True)
class LwCas13a:
    name = 'LwCas13a'
    spacer_length = 28
    scaffold = 'GGGGATTTAGACTACCCCAAAAACGAAGGGGACTAAAAC'
