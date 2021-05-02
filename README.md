# SeniorProject
## Identification and Characterization of Butyrate-Producing Species in the Human Gut Microbiome

## Introduction
In this project we explore genomes from HMP 1 to attempt to identify species that may be capable of butyrate production based on the presence of enzyme-encoding genes involved in the pathway(s) of butyrate production. In this repository, you will find the code and methods used for data formatting, protein alignments, and data analysis.

## DATA:
Species = http://downloads.ihmpdcc.org/data/reference_genomes/body_sites/Gastrointestinal_tract.cds.fsa ("Gastrointestinal_tract.cds.fsa")
Pathway Proteins = https://www.ncbi.nlm.nih.gov/nuccore/DQ987697.1?report=fasta** (Save as "butyrate_genes_BTC&BK.fasta"), 
https://www.ncbi.nlm.nih.gov/nuccore/AY796317.2?report=fasta** (Save as "butyryl_COA_transferase_BCT.fasta"), https://www.ncbi.nlm.nih.gov/protein/WP_003420701.1?report=fasta (Save as "phosphate_butyryltransferase_BK.fasta", https://www.ncbi.nlm.nih.gov/protein/WP_003722496.1?report=fasta (Save as "butyrate_kinase_BK.fasta")
** partial CDS that will be translated to protein sequence during formatting.

## HOW TO USE:
### Step 1 - Data formatting
1. Collect data from HMP
2. Collect data from NCBI
3. run DataFormatting.py --> Creates .fasta files for all protein sequences and all species

### Step 2 - Alignments
1. run RunAlignments.py --> creates xml files for each alignment, csv files for unprocessed results, and csv files for scored results

### Step 3 - Analyses
1. pull 16s rRNA data and taxon info from Silva Database (batch download)
2. run Analyze_Data.py --> Creates output files for labeling taxon info, for formating input ANOVA files, and for running the multiple sequence alignment on Clustal Omega
3. Run ANOVA.R with the output files from step 2
4. Run Clustal Omega multiple sequence alignment with output files from step 2
5. Interpret results with output taxon files from step 2, from the phylogenetic analysis, and from the ANOVA.
