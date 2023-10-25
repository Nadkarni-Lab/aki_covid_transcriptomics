# aki_covid_transcriptomics
This Repository contains code for the paper:

TITLE: **Peripheral Transcriptomics in Acute and Long-Term Kidney Dysfunction in SARS-CoV2 Infection**
AUTHORS: Pushkala Jayaraman , Madhumitha Rajagopal , Ishan Paranjpe , Lora Liharska, Mayte Suarez-Farinas , Ryan Thompson , Diane Marie Del Valle , Noam Beckmann , Wonsuk Oh , Ankit Sakhuja , Annie Chen, Steven Chen, Ephraim Kenigsberg , Akhil Vaid , Edgar Gonzalez-Kozlova , Justin Kauffman , Faris F Gulamali , Sergio Dellepiane , George Vasquez-Rios , John Cijiang He, Steven G Coca , Lili Chan , Eric Schadt, Sacha Gnjatic , Miram Merad, Seunghee Kim-Schulze, Ephraim Tsalik , Raymond Langley , Alexander W Charney , and Girish N Nadkarni 

The scripts are run with R v4.3

Structure:
All scripts are run on our HPC cluster using an LSF job scheduler. 

Main file:
RNASEQQC_newAKICaseCtrl_withcellfrac.R - The script develops a linear model factoring in clinical,  phenotype, and transcriptomic (gene expression) information for COVID-19 patients and finds significant associations with AKI. 
It will begin with creating a results directory with an automated timestamp(start of run time) and will publish all results within the folder. 
It takes in a metadata file and the raw gene expressions to identify differentially expressed genes associated with AKI in COVID-19 patients. 
The output is an rData file in tabular format with the fold change and adjusted pValue significance for each gene. 

associate_protein_with_egfr.R - This script performs a mixed model(lmm) analysis using temporal eGFR values calculated from multiple post-discharge follow-up for patients from the COVID-AKi cohort. The results from the lmm are then analyzed further to find associations between overall eGFR values using a negative beta estimate of the model against tertiles of gene expressions (high-midle-low) that showed significant associations in the lmm.  
