"""
Author: Yixin Zhang
Last updated: Dec. 11, 2023

This file is the main engine of CrisprLib.
"""

import os, sys
import numpy as np
import pandas as pd
from dna import DNA
from interface import Interface
from grna_designer import gRNADesigner
from construction_file_builder import CFBuilder
from utils import CheckFileType, CheckOperations, DNATools

VALID_GENOME = {'1': 'mm10', '2': 'mm39', '3': 'hg19', '4': 'hg38'}

started = False
internet_checked = False
while True:
    try: 
        # Collect preliminary information and initialize required data
        if not started:
            user_interface = Interface()
            chosen_system, chosen_genome, chosen_function, export_file_path = user_interface.start()
            started = True

        # Process imported data
        import_file_path, num = user_interface.import_data(chosen_function)
        grna_designer = gRNADesigner(chosen_system, num)
        file_checker = CheckFileType(os.path.basename(import_file_path))
        if file_checker.is_csv():
            df = pd.read_csv(import_file_path)
        elif file_checker.is_txt():
            df = pd.read_csv(import_file_path, delimiter='\t')
        else:
            df = pd.read_excel(import_file_path)

        # Populate primers

        # Generate gRNAs
        array = np.array(df)
        print('\nGenerating gRNAs and primers...')
        grnas = []
        primers = [DNA('backbone_fwd', 'TTTTTGAATTCTCTAGAGTCGACCTGCAGA', 'primer'),
                   DNA('backbone_rev', 'CGATGACTAGTATTATACCTAGGACT', 'primer')]
        pcr_products = []
        for row in array:
            gene_data = row[1]
            sanity_checker = CheckOperations(gene_data, chosen_genome)
            if gene_data[:3] == 'chr':  # Genomic coordinate
                if not internet_checked:
                    if sanity_checker.check_internet_connection():
                        internet_checked = True
                    else:
                        print('Failed to connect to the Internet. Please make sure you are connected before proceeding.')
                        sys.exit(0)
                chromosome, begin, end = sanity_checker.invalid_coordinate()
                if chromosome[:3] != 'chr':  # Coordinate is invalid
                    print(f'\nYour coordinate for gene named "{row[0]}" is invalid: {chromosome}')
                    raise ValueError()
                if chosen_function == '1':  # Specific gene knockout
                    grna = grna_designer.generate_gRNA_from_coordinate(chromosome, begin, end)
                else:  # Tiled screen
                    grna = grna_designer.tiled_gRNAs_from_coordinate(chromosome, begin, end)
            else:  # Sequence
                gene_data = gene_data.upper()
                message = sanity_checker.invalid_gene()
                if message:  # Gene sequence is invalid
                    print(f'\nYour coordinate for gene named "{row[0]}" is invalid: {message}')
                    raise ValueError()
                if chosen_function == '1':  # Specific gene knockout
                    grna = grna_designer.generate_gRNA_from_seq(gene_data)
                else:  # Tiled screen
                    grna = grna_designer.tiled_gRNAs_from_seq(gene_data)
            grnas.append(grna)
            
            # Generate primers
            x = 1
            for guide in grna:
                fwd, rev = DNATools(guide).generate_primers()
                primers.append(DNA(f'{row[0]}_gRNA{x}_fwd', fwd, 'primer'))
                primers.append(DNA(f'{row[0]}_gRNA{x}_rev', rev, 'primer'))
                pcr_products.append(DNA(f'{row[0]}_gRNA{x}', fwd + guide[15:-14] + rev, 'pcr'))
                x += 1
        print('Generation complete!')
        
        # Output gRNA list upon request
        grna_response, grna_name = user_interface.export_grnas()
        if grna_response == 'YES':  # Generate a list of gRNAs
            exported_grnas = [' '.join(str(g) for g in gene_gRNAs) for gene_gRNAs in grnas]
            num_columns = len(max(exported_grnas, key=len).split())  # Get the maximum number of gRNAs in a gene
            new_columns = {}
            for i in range(num_columns):
                new_columns[f'gRNA {i + 1}'] = [gRNAs.split()[i] if len(gRNAs.split()) > i else '' for gRNAs in exported_grnas]
            new_df = pd.DataFrame(new_columns)
            df = pd.concat([df, new_df], axis=1)
            df.to_csv(f'{export_file_path}/{grna_name}.csv', index=False)
            print('\ngRNAs have been saved to your specified destination.')
        
        # Output primer list upon request
        primer_response, primer_name = user_interface.export_primers()
        primer_names = [primer.name for primer in primers]
        primer_sequences = [primer.sequence for primer in primers]
        if primer_response == 'YES':  # Generate a list of primers
            primer_df = pd.DataFrame({'Primer Name': primer_names, 'Primer Sequence': primer_sequences})
            primer_df.to_csv(f'{export_file_path}/{primer_name}.csv', index=False)
            print('\nPrimers have been saved to your specified destination.')
        
        # Generate construction files
        cf_name = user_interface.export_cf()
        print('\nGenerating construction files...')
        CFBuilder(f'{export_file_path}/{cf_name}').output_CF(primers, pcr_products)
        print('Your files have been successfully exported. Thank you for using the program.\n')
        break

    except (KeyboardInterrupt, EOFError):
        # Handle Ctrl+C and Ctrl+D
        print("\nYou have successfully quitted. Thank you for using the program.\n")
        sys.exit(0)

    except ValueError:
        # Handle bad data structure
        try:
            input("""\nThe uploaded file has incorrect structure or information. Please review README.md for detailed requirements of file format.
If you wish to upload a different file, please press enter to continue.
Otherwise, press control+D simultaneously to quit the program.""")
        except (KeyboardInterrupt, EOFError):
            # Handle Ctrl+C and Ctrl+D
            print("\nYou have successfully quitted. Thank you for using the program.\n")
            sys.exit(0)

'''
    except:
        # Handle crashes
        print("""\n************************************************************************
The program crashed due to an unexpected error. Apologies for your inconvenience.
If you wish to report the bug, please contact czhangyx@berkeley.edu.
************************************************************************n""")
        sys.exit(0)
'''