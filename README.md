# RATT-with-GFF
Pipeline for Rapid Annotation Transfer Tool (RATT) when your starting annotation file is a GFF. 

Utlilizes GFF to EMBL conversion tools and corrects for errors that occur during conversion.

Designed and tested with Linux operating systems.

## PREREQUISITE:


NBISweden/EMBLmyGFF3 (and all prerequisites)
  - https://github.com/NBISweden/EMBLmyGFF3
  
Rapid Annotation Transfer Tool(RATT) (and all prerequisites)
  - http://ratt.sourceforge.net/download.html
  
The European Molecular Biology Open Software Suite (EMBOSS)
  - http://emboss.open-bio.org/html/adm/ch01s01.html
  
Quality Assesment Tools for Genome Assemblies (QUAST)
 - http://bioinf.spbau.ru/en/quast
  
## USAGE:
Requires all files be in the same directory that you are running the script from

> ./RATTwithGFF.py  [reference-GFF] [reference-FASTA] [query-FASTA] [run-ID]  [RATT-Transfer-Type]

**reference-GFF:** Genomic GFF containing the annotations you wish to transfer

**reference-FASTA:** Genomic FASTA containing the sequences corresponding to the reference-GFF

**query-FASTA:** Genomic FASTA that you wish to transfer the annotations to

**run-ID:** A name added as a prefix to each RATT result file

**RATT-Transfer-Type:** The set of parameters you wish RATT to use when transferring annotations
- valid options: Assembly, Assembly.Repetitive, Strain, Strain.Repetitive, Species, Species.Repetitive, Multiple
- See for more info on RATT parameters: http://ratt.sourceforge.net/documentation.html
  
### EXAMPLE
> ./RATTwithGFF.py ref.gff ref.fasta query.fasta refToQuery Strain

## OUTPUT:
**contig_fasta:** contains the reference sequence split into separate fasta files for each contig/chromosome

**contig_gff:** contains the reference annotations split into seperate gff files for each contig/chromosome

**contig_embl:** contains the annotations for each reference contig/chromosome as an EMBL file

**[run-ID]_RATT:** contains the RATT output files organized into subdirectories
- **final_embl** contains the RATT transferred annotations in EMBL format for each query contig/chromosome
- **final_gff** contains the RATT transferred annotation in GFF format for each query contig/chromosome as well as a genomic gff
- **transferStats.csv** comma-delimited file containing the unique feature counts for the reference GFF and newly generated GFF. Also contains statistics related to to the quality of the assembly.
