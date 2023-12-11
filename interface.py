"""
Author: Yixin Zhang
Last updated: Dec. 11, 2023

This file codes for the user interface that is displayed on terminal.
"""

import tkinter as tk
from tkinter import filedialog

VALID_SYSTEMS = {'1': 'SpCas9', '2': 'SaCas9', '3': 'SpCas12a', '4': 'LbCas12a', '5': 'LshCas13a', '6': 'LwCas13a'}
VALID_GENOME = {'1': 'mm10', '2': 'mm39', '3': 'hg19', '4': 'hg38'}
VALID_FUNCTIONS = {'1': 'gene knockout', '2': 'tiled screen'}
VALID_FILETYPES = [("Excel-1", ".xlsx"), ("Excel-2", ".xls"), ("Excel-3", ".xlsm"), ("Excel-4", ".xlsb"),
                   ("Excel-5", ".odf"), ("Excel-6", ".ods"), ("Excel-7", ".odt"),
                   ("CSV", ".csv"), ("TXT", ".txt")]

def start():
    """
    Displays the main menu for user interface. Collects preliminary information.

    Returns:
    chosen_system (str): The Cas system user would like to use.
    chosen_genome (str): The genome user would like to work with.
    chosen_function (str): The experiment user would like to perform.
    export_file_path (str): The destination directory where files will be exported to.
    export_file_name (str): The file name for exported construction file.
    """
    # Internet connection warning
    input("""\n-----------------------------------------------------------
Welcome to CrisprLib, a comprehensive tool for performing CRISPR-related tasks.
    This program may require internet connection.
    Please make sure you are connected to the Internet before proceeding.
    If you wish to quit the program, please press control+D simultaneously.
    Otherwise, press enter to continue.""")

    # Choose CRISPR system
    chosen_system = input("""
    Please choose a supported CRISPR system you would like to use:
    1) SpCas9 (Streptococcus pyogenes)
    2) SaCas9 (Streptococcus aureus)
    3) FnCas12a (Francisella novicida, previously named Cpf1)
    4) LbCas12a (Lachnospiraceae bacterium, previously named Cpf1)
    5) LshCas13a (Leptotrichia shahii)
    6) LwCas13a (Leptotrichia wadei)\n""")
    while chosen_system not in VALID_SYSTEMS:
        chosen_system = input("""Please select a supported system by typing a number.
    If you wish to quit the program, please press control+D simultaneously.
    """)
    print(f'You have chosen to use {VALID_SYSTEMS[chosen_system]}.')
    
    # Choose genome
    chosen_genome = input("""\n-----------------------------------------------------------
    Please choose an organism's genome for perturbation:
    1) Mouse GRCm38/mm10
    2) Mouse GCRm39/mm39
    3) Human GRCh37/hg19
    4) Human GRCh38/hg38\n""")
    while chosen_genome not in VALID_GENOME:
        chosen_genome = input("""Please select a supported organism by typing a number.
    If you wish to quit the program, please press control+D simultaneously.
    """)
    print(f'You have chosen to work with {VALID_GENOME[chosen_genome]}.')

    # Choose experiment
    chosen_function = input("""\n-----------------------------------------------------------
    Please choose what experiment you would like to perform:
    1) Gene knockout for specific genes
    2) Tiled CRISPR screen over sequences\n""")
    while chosen_function not in VALID_FUNCTIONS:
        chosen_function = input("""Please select a supported experiment by typing a number.
    If you wish to quit the program, please press control+D simultaneously.
    """)
    print(f'You have chosen to perform a {VALID_FUNCTIONS[chosen_function]}.')

    # Choose construction file output directory
    input("""\nYou will be asked to select a directory for file output.
    All files generated in the session will go into this directory.
    Press enter to continue.\n""")
    root = tk.Tk()
    root.withdraw()
    export_file_path = filedialog.askdirectory(title="Please select an output directory")
    while not export_file_path:
        input("""\nYou did not select a directory. If you wish to abort, please press control+D simultaneously.
    Otherwise, please press enter to select a directory.\n""")
        export_file_path = filedialog.askdirectory(title="Please select an output directory")

    return chosen_system, chosen_genome, chosen_function, export_file_path

def import_data(function):
    """
    User imports data in this step. The data should follow the structures detailed in README.md.

    Args:
    function (str): The experiment user chooses to perform in the main menu.

    Returns:
    import_file_path (str): The file path to the imported data.
    grna_num (int) or spacing (int): The number of gRNAs should be generated; The spacing between each tiled gRNA.
    """
    input("""\nYou will be asked to select a file to import genes you would like to target. Please review README.md for detailed requirements of file format before proceeding.
    Press enter to continue.\n""")
    root = tk.Tk()
    root.withdraw()
    import_file_path = filedialog.askopenfilename(title="Please select a file", filetypes=VALID_FILETYPES)
    while not import_file_path:
        input("""\nYou did not select a file. If you wish to abort, please press control+D simultaneously.
    Otherwise, please press enter to select a file.\n""")
        import_file_path = filedialog.askopenfilename(title="Please select a file", filetypes=VALID_FILETYPES)
    if function == '1':  # Gene knockout
        grna_num = input('\nPlease enter the number of gRNAs you want to generate for each gene: ')
        while True:
            try:
                grna_num = int(grna_num)
                if grna_num:
                    break
                grna_num = input('\nPlease enter a valid number of base pairs you want as spacing between gRNAs: ')
                continue
            except:
                grna_num = input('\nPlease enter a valid number of base pairs you want as spacing between gRNAs: ')
        return import_file_path, grna_num
    else:  # Tiled screen
        spacing = input('\nPlease enter a numeric tile size: ')
        while True:
            try:
                spacing = int(spacing)
                if spacing:
                    break
                spacing = input('\nPlease enter a valid tile size: ')
                continue
            except:
                spacing = input('\nPlease enter a valid tile size: ')
        return import_file_path, spacing

def export_grnas():
    """
    Ask user if they want a copy of gRNA sequences.

    Returns:
    response (str): Whether the user would like a copy of gRNAs.
    grna_name (str): The file name for exported gRNAs.
    """
    response = input('\nWould you like to receive a copy of gRNAs generated for you? Please enter yes or no.\n')
    while True:
        try:
            response = response.upper()
            if response == 'YES' or response == 'NO':
                break
            response = input('\nPlease enter yes if you want to receive a copy of gRNAs or enter no otherwise.\n')
        except:
            response = input('\nPlease enter yes if you want to receive a copy of gRNAs or enter no otherwise.\n')
    
    if response == 'NO':
        return response, None
    # Choose name of output file
    print('\nThe file will be exported to your previously selected folder.')
    grna_name = input("Please enter a name for the output file (without the extension):\n")
    while not grna_name.strip():
        grna_name = input("Name cannot be empty. Please enter a name for the output file (without the extension):\n")
    return response, grna_name

def export_primers():
    """
    Ask user if they want a copy of primer sequences.

    Returns:
    response (str): Whether the user would like a copy of primers.
    primer_name (str): The file name for exported primers.
    """
    response = input('\nWould you like to receive a copy of primers generated for you? Please enter yes or no.\n')
    while True:
        try:
            response = response.upper()
            if response == 'YES' or response == 'NO':
                break
            response = input('\nPlease enter yes if you want to receive a copy of primers or enter no otherwise.\n')
        except:
            response = input('\nPlease enter yes if you want to receive a copy of primers or enter no otherwise.\n')
    
    if response == 'NO':
        return response, None
    # Choose name of output file
    print('\nThe file will be exported to your previously selected folder.')
    primer_name = input("Please enter a name for the output file (without the extension):\n")
    while not primer_name.strip():
        primer_name = input("Name cannot be empty. Please enter a name for the output file (without the extension):\n")
    return response, primer_name
    
def export_cf():
    """
    Get a file name for exported construction files.

    Returns:
    cf_name (str): The file name for exported construction files.
    """
    print('\nConstruction files in readable form and json form will be exported.')
    cf_name = input("Please enter a name for the output file (without the extension):\n")
    while not cf_name.strip():
        cf_name = input("Name cannot be empty. Please enter a name for the output file (without the extension):\n")
    return cf_name
