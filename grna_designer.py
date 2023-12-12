"""
Author: Yixin Zhang
Last updated: Dec. 11, 2023

This file contains functions for designing gRNAs for various Cas systems.
"""

import requests
from cas import SpCas9, SaCas9, FnCas12a, LbCas12a, LshCas13a, LwCas13a
from utils import DNATools

class gRNADesigner:
    """
    A class for designing gRNAs for various CRISPR systems.

    Args:
    system (str): The selected CRISPR system.
    num (int): The number of gRNAs should be generated for gene knockout,
               or the tile size between each gRNA for CRISPR screen.
    """
    
    def __init__(self, system, num) -> None:
        if system == '1':
            self.cas = SpCas9()
        if system == '2':
            self.cas = SaCas9()
        if system == '3':
            self.cas = FnCas12a()
        if system == '4':
            self.cas = LbCas12a()
        if system == '5':
            self.cas = LshCas13a()
        if system == '6':
            self.cas = LwCas13a()
        self.num = num

    def generate_gRNA_from_seq(self, gene):
        """
        Generate gRNAs for a gene with scaffold attached.

        Args:
        gene (str): The gene sequence.

        Returns:
        list: A list of lists, where each sublist contains gRNAs.
        """
        if 'Cas13a' in self.cas.name:
            return self.generate_PAMless_gRNA(gene)
        
        grnas = []
        prev_start = 0
        current_start = 0
        i = 0
        if 'Cas9' in self.cas.name:
            while i < len(gene)-self.cas.pam_length:
                if len(grnas) == self.num:  # Return if there are enough gRNAs
                    return grnas
                potential_pam = gene[i:i+self.cas.pam_length]
                spacer = None
                if i >= self.cas.spacer_length+self.cas.n_length and potential_pam in self.cas.pams:  # Check for upstream spacer
                    spacer = gene[i-self.cas.n_length-self.cas.spacer_length:i-self.cas.n_length]
                    current_start = i-self.cas.n_length-self.cas.spacer_length
                if not spacer and i <= len(gene)-self.cas.spacer_length-self.cas.n_length-self.cas.pam_length:  # Check for downstream spacer
                    rc_pams = []
                    for pam in self.cas.pams:
                        rc_pams.append(DNATools(pam).reverse_complement())
                    if potential_pam in rc_pams:
                        spacer = gene[i+self.cas.pam_length+self.cas.n_length:i+self.cas.pam_length+self.cas.n_length+self.cas.spacer_length]
                        current_start = i+self.cas.pam_length+self.cas.n_length
                i += 1
                if not spacer or current_start - prev_start <= self.cas.spacer_length:
                    continue
                grnas.append(spacer + self.cas.scaffold)
                prev_start = current_start
                i += self.cas.spacer_length  # Guarantees distinct gRNAs
        else:  # Cas12a
            while i < len(gene)-self.cas.pam_length:
                if len(grnas) == self.num:  # Return if there are enough gRNAs
                    return grnas
                potential_pam = gene[i:i+self.cas.pam_length]
                spacer = None
                if i <= len(gene)-self.cas.spacer_length-self.cas.n_length-self.cas.pam_length and potential_pam in self.cas.pams:  # Check for downstream spacer
                    spacer = gene[i+self.cas.pam_length+self.cas.n_length:i+self.cas.pam_length+self.cas.n_length+self.cas.spacer_length]
                    current_start = i-self.cas.n_length-self.cas.spacer_length
                if not spacer and i >= self.cas.spacer_length+self.cas.n_length:  # Check for upstream spacer
                    rc_pams = []
                    for pam in self.cas.pams:
                        rc_pams.append(DNATools(pam).reverse_complement())
                    if potential_pam in rc_pams:
                        spacer = gene[i-self.cas.n_length-self.cas.spacer_length:i-self.cas.n_length]
                        current_start = i+self.cas.pam_length+self.cas.n_length
                i += 1
                if not spacer or current_start - prev_start <= self.cas.spacer_length:
                    continue
                prev_start = current_start
                grnas.append(self.cas.scaffold + spacer + 'TTTTTT')
                i += self.cas.spacer_length  # Guarantees distinct gRNAs
        return grnas
                
    def generate_gRNA_from_coordinate(self, chromosome, begin, end):
        """
        Generate gRNAs for a gene at the specified location with scaffold attached.

        Args:
        chromosome (str): The chromosome where the gene belongs.
        begin (int): The gene's starting location in the chromosome.
        end (int): The gene's ending location in the chromosome.

        Returns:
        list: A list of lists, where each sublist contains gRNAs.
        """
        data = requests.get(f"https://api.genome.ucsc.edu/getData/sequence?genome=mm10;chrom={chromosome};start={begin};end={end}", timeout=30)
        gene = data.json()['dna'].upper()
        return self.generate_gRNA_from_seq(gene.upper())

    def generate_PAMless_gRNA(self, gene):
        """
        Generate gRNAs for a gene with scaffold attached for an RNA-cleaving Cas system.

        Args:
        gene (str): The gene sequence.

        Returns:
        list: A list of lists, where each sublist contains gRNAs.
        """
        spacers = []
        for i in range(len(gene)-self.cas.spacer_length+1):
            window = gene[i:i+self.cas.spacer_length]
            g = DNATools(window).hairpin_counter()
            spacers.append((window, g))
        spacers = sorted(spacers, key=lambda x: x[1])[:self.num]
        rc_spacers = []
        for spacer in spacers:
            rc_spacers.append(DNATools(spacer[0]).reverse_complement())
        return [self.cas.scaffold + spacer for spacer in rc_spacers]

    def tiled_gRNAs_from_seq(self, gene):
        """
        Generate tiled gRNAs with specified spacing for a gene with scaffold attached.

        Args:
        gene (str): The gene sequence.

        Returns:
        list: A list of lists, where each sublist contains gRNAs.
        """
        grnas = []
        i = 0
        while i <= len(gene)-self.cas.spacer_length:
            if 'Cas9' in self.cas.name:
                grnas.append(gene[i:i+self.cas.spacer_length] + self.cas.scaffold)
            else:  # Cas12a and Cas13a
                grnas.append(self.cas.scaffold + gene[i:i+self.cas.spacer_length] + 'TTTTTT')
            i += self.num
        return grnas

    def tiled_gRNAs_from_coordinate(self, chromosome, begin, end):
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
        return self.tiled_gRNAs_from_seq(gene.upper())
