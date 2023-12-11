# CrisprLib, a comprehensive tool for generating lab guidelines for CRISPR projects.

## Usage
This program generates gRNAs, appropriate primers, and experiment procedures for CRISPR-related projects. The program prompts user for the CRISPR/Cas systems, genome, and experiment. Upon request, the program outputs a table of gRNA sequences in the selected destination folder. The program then outputs a table of primers required to execute experiments, as well as human-readable and json [construction files](https://doi.org/10.1101/2023.06.28.546630).

To execute the program, navigate to the crisprlib directory on your terminal and enter the following command:
```python crisprlib.py```
You can then follow the instructions in the program.

### Supported CRISPR/Cas systems:
- SpCas9 (Streptococcus pyogenes)
- SaCas9 (Streptococcus aureus)
- FnCas12a (Francisella novicida, previously named Cpf1)
- LbCas12a (Lachnospiraceae bacterium, previously named Cpf1)
- LshCas13a (Leptotrichia shahii)
- LwCas13a (Leptotrichia wadei)
### Supported genomes:
- Mouse GRCm38/mm10
- Mouse GCRm39/mm39
- Human GRCh37/hg19
- Human GRCh38/hg38

## Dependencies
This program require NumPy and Pandas. Both packages can be installed via pip.
The program may require Internet connection if the input data contains genomic coordinates.
The construction file outlines a set of experiments based on [this plasmid](https://www.addgene.org/62226/). The gRNAs should be ordered as ssDNA and amplified with PCR. They should then be digested with SpeI and EcoRI and ligated with SpeI/EcoRI digested backbone and transformed.

## Input file structure
Accepted file formats: .csv, .txt, .xlsx, .xls, .xlsm, .xlsb, .odf, .ods, .odt
The file should only have two columns. The first column contains distinct names of the genes, and the second column contains gene sequences or genomic coordinates. Gene sequences should be at least 30 nucleotides long and should only include A/T/C/G. Your data should be free of space, punctuations, and other irrelevant characters. There should be a header row in the file.

### References:
Genome data source: [UCSC Genome Browser](https://genome.ucsc.edu)  
[Enabling AI in Synthetic Biology through Construction File Specification](https://doi.org/10.1101/2023.06.28.546630)  
[DNA targeting specificity of RNA-guided Cas9 nucleases](https://doi.org/10.1038/nbt.2647)  
[In vivo genome editing using Staphylococcus aureus Cas9](https://doi.org/10.1038/nature14299)  
[Cpf1 Is a Single RNA-Guided Endonuclease of a Class 2 CRISPR-Cas System](https://doi.org/10.1016/j.cell.2015.09.038)  
[A Cas12a ortholog with stringent PAM recognition followed by low off-target editing rates for genome editing](https://doi.org/10.1186/s13059-020-01989-2)  
[Improvement of the LbCas12a-crRNA System for Efficient Gene Targeting in Tomato](https://doi.org/10.3389/fpls.2021.722552)  
[Structure-based design of gRNA for Cas13](https://doi.org/10.1038/s41598-020-68459-4)  
[Nucleic acid detection with CRISPR-Cas13a/C2c2](https://doi.org/10.1126/science.aam9321)  

Author: Yixin Zhang
To report bugs, please send an email to czhangyx@berkeley.edu.