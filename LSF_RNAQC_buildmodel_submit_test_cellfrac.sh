#!/bin/bash

#BSUB -J AKICVD_RNAQC
#BSUB -P "acc_EHR_ML"
#BSUB -q express 
#BSUB -n 12 
#BSUB -R span[hosts=1]
#BSUB -W 12:00 
#BSUB -o %J.stdout 
#BSUB -eo %J.stderr 
#BSUB -oo RNAQC_%J.out

module load R/4.0.3 

#Rscript RNASEQQC_lbp_allBatches_QCwithDream_forPushkala_wfactors_buildmodel_test.R

Rscript RNASEQQC_lbp_allBatches_QCwithDream_forPushkala_wfactors_buildmodel_test_withcellfrac.R 
