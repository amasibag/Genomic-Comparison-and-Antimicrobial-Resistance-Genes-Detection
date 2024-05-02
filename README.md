# Genomic-Comparison-and-Antimicrobial-Resistance-Genes-Detection
The primary aim of this project is to perform a comparative genomic analysis between two bacterial strains to identify and characterize antimicrobial resistance genes. By leveraging genome annotation data from RAST, the project seeks to pinpoint differences and commonalities in ARG profiles.

This project performs gene function analysis for comparative genomics studies between two strains of Vibrio coralliilyticus, namely RE22 and ATCC BAA-450. V. coralliilyticus is a Gram-negative, rod-shaped bacterium, that is known to infect coral causing  tissue lysis, bleaching, and drastic losses in coral reefs. It is considered to be one of the major pathogens due to its ability to cause major economic damage by infecting various aquatic organisms including oysters. V. coralliilyticus ATCC BAA-450 is the type-strain while RE22 specifically infects larval bivalve shellfish and causes high larval mortality in oyster hatcheries. Analysis of the genome of these two strains and determination of antimicrobial resistance genes present may shed light on the mechanism as to why these strains are difficult to eradicate in oyster aquacultures and aid in the choice of proper antimicrobials to be used as well as in probiotic development. This script encompasses various stages of data manipulation, statistical summaries, and visualizations, making extensive use of R libraries such as dplyr, ggplot2, and igraph. The analysis includes loading data, identifying common gene functions, examining gene structure, assessing gene orientations and locations, identifying antimicrobial resistance genes (ARG), and constructing a network of shared and unique ARGs.  

The genome of the bacterial strains used in this project was downloaded from NCBI GenBank. ATCC BAA-450 is the type-strain and the RE22 is the strain commonly associated with oysters.

Rapid Annotation using Subsystem Technology (RAST) was also used to annotate the bacterial genomes. These annotations include information such as gene functions, locations, and other metadata essential for downstream analysis.
    Parameters used for RAST:
    Genetic code:	11
    Annotation scheme:	RASTtk
    Preserve gene calls:	no
    Automatically fix errors:	yes
    Fix frameshifts:	yes
    Backfill gaps:	yes

Annotated CSV files were loaded to R studio.

Gene Size was performed by getting the difference between the start and stop of each gene.

Most common Protein-encoding gene (PEGS) functions were identified by filtering the data using "PEGS" and counting the occurence of each gene function.

Gene orientation and cluster analysis were done using simple count function and ggplot.

Identification and analysis of antimicrobial resistance genes (ARGs) present in one or both strain were done using keyword search, the script filters for genes related to "resistance". *Other ARG detection packages were not use on this project due to the limitation set by the R studio version used. Other current packages tried for ARG detection does not function well with R studio version 4.1.

Network of ARGs is still a work in progress. Other gene functions can also be mapped using different visualizations and parameters.

R script were converted to a RMD file. 


