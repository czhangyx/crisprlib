"""
Author: Yixin Zhang
Last updated: Dec. 11, 2023

This file contains utility functions for CrisprLib.
"""

import os, requests

LETTER_CHROMOSOMES = 'XYM'
STOP_CODONS = ['TAA', 'TAG', 'TGA']
FWD_TAIL = 'CCATAACTAGT'  # 5' hang + SpeI cut site
REV_TAIL = 'CTCAGGAATTC'  # 5' hang + EcoRI cut site

class CheckFileType:
    """
    A class for checking input file's extension.

    Args:
    filename (str): The name of the file.
    """
    def __init__(self, filename):
        _, self.file_extension = os.path.splitext(filename)
        self.file_extension = self.file_extension.lower()

    def is_csv(self):
        """
        Check if a file has a CSV format based on its file extension.

        Returns:
        bool: True if the file has a CSV format, False otherwise.
        """
        return self.file_extension == '.csv'

    def is_txt(self):
        """
        Check if a file has a TXT format based on its file extension.

        Returns:
        bool: True if the file has a TXT format, False otherwise.
        """
        return self.file_extension == '.txt'

class DNATools:
    """
    A class for performing calculations on a DNA sequence.

    Args:
    sequence (str): The input DNA sequence (e.g., "ATGC").
    """
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    def __init__(self, sequence) -> None:
        self.sequence = sequence.upper()

    def reverse_complement(self):
        """
        Calculate the reverse complement of a DNA sequence.

        Returns:
        str: The reverse complement DNA sequence.
        """
        reverse_sequence = self.sequence[::-1]
        reverse_complement_sequence = ''.join(DNATools.complement_dict[base] for base in reverse_sequence)
        return reverse_complement_sequence

    def hairpin_counter(self):
        """
        Calculates a score for folding for all combined hairpins in a sequence.

        Returns:
        float: a deltaG of stabilization of haipins

        Adapted from Project 3 code.
        """
        rcSeq = self.reverse_complement()
        charSeq = list(self.sequence)
        len_seq = len(self.sequence)
        out = 0.0

        for spaces in range(4, 10):  # For each spacer length between 4 and 9
            for i in range(len_seq - spaces - 12):
                hbonds = self._count_hbonds(charSeq, i, i + spaces + 12, rcSeq, len_seq)
                out += 2 ** hbonds - 1

        return out

    def _count_hbonds(self, charSeq, start_inc, end_excl, revC, len_seq):
        """
        Counts the number of hydrogen bonds participating in hairpins.

        Returns:
        hbonds (int): The number of hydrogen bonds.
        
        Adapted from Project 3 code.
        """
        suffix_start = end_excl - 6
        match_length = 0

        for dist_from_end in range(6):
            rev_char = revC[len_seq - 1 - start_inc - 5 + dist_from_end]
            for_char = charSeq[suffix_start + dist_from_end]
            if for_char == rev_char:
                match_length = dist_from_end
            else:
                break

        if match_length < 2:
            return 0

        hbonds = 0
        for i in range(match_length):
            achar = revC[len_seq - 1 - start_inc - 5 + i]
            if achar == 'C' or achar == 'G':
                hbonds += 3
            if achar == 'A' or achar == 'T':
                hbonds += 2

        return hbonds

    def generate_primers(self):
        """
        Design forward and reverse primers to amplify the input sequence.

        Returns:
        fwd (str): Forward primer of the input sequence.
        rev (str): Reverse primer of the input sequence.
        """
        return FWD_TAIL+self.sequence[:15], REV_TAIL+self.reverse_complement()[:15]

class CheckOperations:
    """
    A class for checking basic operation conditions.

    Args:
    data (str): The input DNA sequence (e.g., "ATGC") or genomic coordinate.
    genome (str): The selected genome.
    """
    def __init__(self, data, genome) -> None:
        self.data = data
        self.genome = genome

    def invalid_gene(self):
        """
        Check if the input gene is valid.

        Returns:
        str: The error message if the input sequence is invalid,
            or
        None: If the input sequence is valid.
        """
        if any([char.upper() not in 'ATCG' for char in self.data]):
            return "Gene has invalid nucleotides"
        if len(self.data) < 30:
            return "Gene cannot be shorter than 30 nucleotides"

    def invalid_coordinate(self):
        """
        Check if a genomic coordinate is valid.

        Returns:
        str: The error message if the input coordinate is invalid,
            or
        str: The chromosome label if the input sequence is valid.
        """
        chromosome = ""
        start = ""
        end = ""
        chr_num = False
        start_done = False
        for char in self.data[3:]:
            # Parse coordinate
            if char == ":":
                chr_num = True
                continue
            if char == "-":
                start_done = True
                continue
            if not chr_num:
                chromosome += char
                continue
            if not start_done:
                start += char
                continue
            end += char
        
        # Check coordinate position
        try:
            start = int(start)
            end = int(end)
        except:
            return "Coordinate positions should be integers", None, None

        # Start position should be before end position
        if start >= end:
            return "Coordinate has invalid range", None, None
        
        # Check if chromosome number is valid
        try:
            if int(chromosome) < 1:
                return "Coordinate has invalid chromosome", None, None
            if self.genome == "1" or self.genome == "2":  # Mouse
                if int(chromosome) > 19:
                    return "Coordinate has invalid chromosome", None, None
            else:
                if int(chromosome) > 22:
                    return "Coordinate has invalid chromosome", None, None
        except:
            chromosome = chromosome.upper()
            if chromosome not in LETTER_CHROMOSOMES:
                return "Coordinate has invalid chromosome", None, None

        return 'chr'+chromosome, start, end

    def check_internet_connection(self):
        """
        Check for an internet connection.

        Return:
        bool: True if the connection is successful, False otherwise.
        """
        try:
            # Attempt to make a request to the specified URL
            response = requests.get('http://www.google.com/', timeout=5)
            # Return True if the request was successful
            return response.status_code == 200
        except requests.ConnectionError:
            # Return False if a connection error occurs
            return False
