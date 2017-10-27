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
  
## USAGE:

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
