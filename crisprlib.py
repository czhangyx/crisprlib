"""
Author: Yixin Zhang
Last updated: Dec. 11, 2023

This file is the main engine of CrisprLib.
"""

import os, sys
import numpy as np
import pandas as pd
from interface import start, import_data, export_grnas, export_primers, export_cf
from utils import is_csv, is_txt, invalid_gene, invalid_coordinate, generate_primers, check_internet_connection
from grna_designer import generate_gRNA_from_seq, generate_gRNA_from_coordinate
from grna_designer import tiled_gRNAs_from_seq, tiled_gRNAs_from_coordinate
from construction_file_builder import output_CF

VALID_GENOME = {'1': 'mm10', '2': 'mm39', '3': 'hg19', '4': 'hg38'}

started = False
internet_checked = False
while True:
    try: 
        # Collect preliminary information and initialize required data
        if not started:
            chosen_system, chosen_genome, chosen_function, export_file_path = start()
            started = True

        # Process imported data
        import_file_path, num = import_data(chosen_function)
        if is_csv(os.path.basename(import_file_path)):
            df = pd.read_csv(import_file_path)
        elif is_txt(os.path.basename(import_file_path)):
            df = pd.read_csv(import_file_path, delimiter='\t')
        else:
            df = pd.read_excel(import_file_path)

        # Generate gRNAs
        array = np.array(df)
        print('\nGenerating gRNAs and primers...')
        grnas = []
        primer_names = ['backbone_fwd', 'backbone_rev']
        primer_sequences = ['TTTTTGAATTCTCTAGAGTCGACCTGCAGA', 'CGATGACTAGTATTATACCTAGGACT']
        pcr_product = []
        for row in array:
            gene_data = row[1]
            if gene_data[:3] == 'chr':  # Genomic coordinate
                if not internet_checked:
                    if check_internet_connection():
                        internet_checked = True
                    else:
                        print('Failed to connect to the Internet. Please make sure you are connected before proceeding.')
                        sys.exit(0)
                chromosome, begin, end = invalid_coordinate(gene_data, chosen_genome)
                if chromosome[:3] != 'chr':  # Coordinate is invalid
                    print(f'\nYour coordinate for gene named "{row[0]}" is invalid: {chromosome}')
                    raise ValueError()
                if chosen_function == '1':  # Specific gene knockout
                    grna = generate_gRNA_from_coordinate(chosen_system, chromosome, begin, end, num)
                else:  # Tiled screen
                    grna = tiled_gRNAs_from_coordinate(chosen_system, chromosome, begin, end, num)
            else:  # Sequence
                gene_data = gene_data.upper()
                message = invalid_gene(gene_data)
                if message:  # Gene sequence is invalid
                    print(f'\nYour coordinate for gene named "{row[0]}" is invalid: {message}')
                    raise ValueError()
                if chosen_function == '1':  # Specific gene knockout
                    grna = generate_gRNA_from_seq(gene_data, chosen_system, num)
                else:  # Tiled screen
                    grna = tiled_gRNAs_from_seq(gene_data, chosen_system, num)
            grnas.append(grna)
            
            # Generate primers
            x = 1
            for guide in grna:
                primer_names.extend([f'{row[0]}_gRNA{x}_fwd', f'{row[0]}_gRNA{x}_rev'])
                fwd, rev = generate_primers(guide)
                primer_sequences.extend([fwd, rev])
                pcr_product.append(fwd + guide[15:-14] + rev)
                x += 1
        print('Generation complete!')
        
        # Output gRNA list upon request
        grna_response, grna_name = export_grnas()
        if grna_response == 'YES':  # Generate a list of gRNAs
            exported_grnas = [' '.join(str(g) for g in gene_gRNAs) for gene_gRNAs in grnas]
            num_columns = len(max(exported_grnas, key=len).split())  # Get the maximum number of gRNAs in a gene

            # Prepare a dictionary to hold the new columns
            new_columns = {}
            for i in range(num_columns):
                new_columns[f'gRNA {i + 1}'] = [gRNAs.split()[i] if len(gRNAs.split()) > i else '' for gRNAs in exported_grnas]

            # Create a new DataFrame from the new columns
            new_df = pd.DataFrame(new_columns)

            # Concatenate the new DataFrame with the original one
            df = pd.concat([df, new_df], axis=1)

            df.to_csv(f'{export_file_path}/{grna_name}.csv', index=False)
            print('\ngRNAs have been saved to your specified destination.')
        
        # Output primer list upon request
        primer_response, primer_name = export_primers()
        if primer_response == 'YES':  # Generate a list of primers
            primer_df = pd.DataFrame({'Primer Name': primer_names, 'Primer Sequence': primer_sequences})
            primer_df.to_csv(f'{export_file_path}/{primer_name}.csv', index=False)
            print('\nPrimers have been saved to your specified destination.')
        
        # Generate construction files
        cf_name = export_cf()
        print('\nGenerating construction files...')
        output_CF(primer_names, primer_sequences, pcr_product, f'{export_file_path}/{cf_name}')
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
