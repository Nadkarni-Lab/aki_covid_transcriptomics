#.libPaths( c( .libPaths(), "/hpc/users/jayarp02/R/lib/4.0.3/") )
.libPaths( c( .libPaths(), "/Users/jayaramanp/workspace/R/lib/4.04/") )
.libPaths()

library("dplyr")
library("stringr")
library("lubridate")

setwd(getwd())

### Load somalogic datas
#source("/sc/arion/projects/EHR_ML/pushkala_work/AKI_COVID_WGCNA/tag_aki_covid.R")
#write.csv(file = "somalogic_input_df.csv", somalogic)
somalogic = readRDS("somalogic_inputdf.rds")
lnames = load("somalogic_after_aki_encounter_clinical_data_with_protein.RData")
lnames
dim(somalogic)
dim(final_df)
colnames(final_df)
somalogic_protein = final_df[,1:30]

#data.frame(somalogic %>% filter(Uncorrected_Blood_Sample=="PICR7118T1") %>% select(Subject_ID, Uncorrected_Blood_Sample, aki, aki_stage, COVID19_Within_Encounter, Days_Since_Admission, RNA_Corrected_Sample))

#rnaseq = somalogic[somalogic$Uncorrected_Blood_Sample %in% somalogic_protein$Uncorrected_Blood_Sample ,]
#rnaseq = somalogic %>% 
#  filter(Uncorrected_Blood_Sample %in% somalogic_protein$Uncorrected_Blood_Sample)

rnaseq = somalogic %>% 
  filter(Uncorrected_Blood_Sample %in% somalogic_protein$Uncorrected_Blood_Sample & 
           between(Event_Date, somalogic_protein$somalogic_date + 2, somalogic_protein$somalogic_date - 2))
data.frame(rnaseq$Uncorrected_Blood_Sample)

##remove not duplicated samples - PICR7118T1, PICR7392T4 - they both have multiple rnaseq IDs but neither have aki and they all had covid, so just pick one. 
rnaseq_not_dup = data.frame(rnaseq[!duplicated(rnaseq$Uncorrected_Blood_Sample),])
#dim(rnaseq_not_dup)
#dim(somalogic_protein)
identical(sort(rnaseq_not_dup$Uncorrected_Blood_Sample), sort(somalogic_protein$Uncorrected_Blood_Sample))
rnaseq_samples = rnaseq_not_dup %>% filter(!is.na(rnaseq_not_dup$RNA_Corrected_Sample) & !grepl('Pediatric', rnaseq_not_dup$RNA_Corrected_Sample))
#table(rnaseq_samples$aki, rnaseq_samples$COVID19_Most_Recent_Order_Result)
somalogic_protein_ordered = somalogic_protein[order(somalogic_protein$Uncorrected_Blood_Sample),]
rnaseq_samples_ordered = rnaseq_samples[order(rnaseq_samples$Uncorrected_Blood_Sample),] 
rownames(rnaseq_samples_ordered) <- rnaseq_samples_ordered$Uncorrected_Blood_Sample
#rnaseq_merged = merge(somalogic_protein_ordered, rnaseq_samples_ordered, on = "Uncorrected_Blood_Sample")
rnaseq_merged = merge(somalogic_protein_ordered, rnaseq_samples_ordered, by = "Uncorrected_Blood_Sample")
save(rnaseq_samples_ordered, file = "rnaseq_sample_for_analysis.rData")
save(rnaseq_merged, file = "rnaseq_final_merged.rData")

