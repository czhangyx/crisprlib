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

def is_csv(filename):
    """
    Check if a file has a CSV format based on its file extension.

    Args:
    filename (str): The name of the file.

    Returns:
    bool: True if the file has a CSV format, False otherwise.
    """
    _, file_extension = os.path.splitext(filename)
    return file_extension.lower() == '.csv'

def is_txt(filename):
    """
    Check if a file has a TXT format based on its file extension.

    Args:
    filename (str): The name of the file.

    Returns:
    bool: True if the file has a TXT format, False otherwise.
    """
    _, file_extension = os.path.splitext(filename)
    return file_extension.lower() == '.txt'

def reverse_complement(sequence):
    """
    Calculate the reverse complement of a DNA sequence.

    Args:
    sequence (str): The input DNA sequence (e.g., "ATGC").

    Returns:
    str: The reverse complement DNA sequence.
    """
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_sequence = sequence[::-1]
    reverse_complement_sequence = ''.join(complement_dict[base] for base in reverse_sequence)
    return reverse_complement_sequence

def hairpin_counter(seq):
    """
    Calculates a score for folding for all combined hairpins in a sequence.

    Args:
    sequence (str): A gene sequence.

    Returns:
    float: a deltaG of stabilization of haipins

    Adapted from Project 3 code.
    """
    seq = seq.upper()
    rcSeq = reverse_complement(seq)
    charSeq = list(seq)
    len_seq = len(seq)
    out = 0.0

    for spaces in range(4, 10):  # For each spacer length between 4 and 9
        for i in range(len_seq - spaces - 12):
            hbonds = count_hbonds(charSeq, i, i + spaces + 12, rcSeq, len_seq)
            out += 2 ** hbonds - 1

    return out

def count_hbonds(seq, start_inc, end_excl, revC, len_seq):
    """
    Counts the number of hydrogen bonds participating in hairpins.
    
    Adapted from Project 3 code.
    """
    suffix_start = end_excl - 6
    match_length = 0

    for dist_from_end in range(6):
        rev_char = revC[len_seq - 1 - start_inc - 5 + dist_from_end]
        for_char = seq[suffix_start + dist_from_end]
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

def generate_primers(sequence):
    """
    Design forward and reverse primers to amplify the input sequence.

    Args:
    sequence (str): A gene sequence.

    Returns:
    fwd (str): Forward primer of the input sequence.
    rev (str): Reverse primer of the input sequence.
    """
    return FWD_TAIL+sequence[:15], REV_TAIL+reverse_complement(sequence)[:15]

def invalid_gene(gene):
    """
    Check if a gene is valid.

    Args:
    gene (str): The input gene sequence.

    Returns:
    str: The error message if the input sequence is invalid,
        or
    None: If the input sequence is valid.
    """
    if any([char.upper() not in 'ATCG' for char in gene]):
        return "Gene has invalid nucleotides"
    if len(gene) < 30:
        return "Gene cannot be shorter than 30 nucleotides"

def invalid_coordinate(coordinate, genome):
    """
    Check if a genomic coordinate is valid.

    Args:
    coordinate (str): The input genomic coordinate.
    genome (str): The selected genome where the coordinate belongs.

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
    for char in coordinate[3:]:
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
        if genome == "1" or genome == "2":  # Mouse
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

def check_internet_connection(url='http://www.google.com/', timeout=5):
    """
    Check for an internet connection.

    Args:
    url (str): URL to test the connection.
    timeout (int): Time in seconds to wait for a response.

    Return:
    bool: True if the connection is successful, False otherwise.
    """
    try:
        # Attempt to make a request to the specified URL
        response = requests.get(url, timeout=timeout)
        # Return True if the request was successful
        return response.status_code == 200
    except requests.ConnectionError:
        # Return False if a connection error occurs
        return False
