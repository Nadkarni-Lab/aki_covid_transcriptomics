#.libPaths( c( .libPaths(), "/Users/jayaramanp/workspace/R/lib/4.04/") )
.libPaths( c( .libPaths(), "/hpc/users/jayarp02/R/lib/4.0.3/") ) 

library("tidyverse")
library("stringr")
library("lubridate")
library("transplantr")
library("ggplot2")
library("plotly")
#install.packages("lcsm")
#devtools::install_github("milanwiedemann/lcsm")
library("BiocManager")
#BiocManager::install("lcsm", lib="/Users/jayaramanp/workspace/R/lib/4.04/")
#library("lcsm")    ### does not work with R/4.0.3 not sure why.. however we do not need it right now. 
library("lme4")
library("lmerTest")
library("foreach")
library("doParallel")
#setwd("/Users/iparanjpe/Documents/mount_sinai_research/covid19/biobank_somalogic/")
setwd("./")

# load("data/clarity/pheno_df_egfr_merged_final_with_discharge_date.RData")
#load("/Users/jayaramanp/Dropbox/MSSM/LAB/AKI_COVID_WGCNA/AKI_COVID_ATTEMPTS/Manuscript/comparison_with_proteomics/somalogic_after_aki_encounter_clinical_data_with_protein (1).RData")


load("pheno_df_egfr_merged_final (7).RData")
#pheno_orig = data.frame(get(load("pheno_df_egfr_merged_final (7).RData")))
final_df_rnaseq = readRDS("geneexpr.rds")
final_df_transposed_rnaseq = data.frame(t(final_df_rnaseq))
final_df_transposed_rnaseq$Subject_ID = data.frame(t(data.frame(strsplit(rownames(final_df_transposed_rnaseq), "T"))))$X1

#load("data/preprocessed_data/somalogic_after_aki_encounter_clinical_data_with_protein.RData")
metadata = readRDS("./metadata.rds")
metadata$Subject_ID = data.frame(t(data.frame(strsplit(metadata$RNA_Corrected_Sample, "T"))))$X1

table(pheno$max_aki_stage)
pheno = pheno[, -c(41:4536)]

common_pheno_finaldf_cols = intersect(colnames(pheno), colnames(final_df_transposed_rnaseq))
pheno_data = merge(pheno, final_df_transposed_rnaseq, by = "Subject_ID")
pheno = pheno_data

#### add baseline eGFR
pheno$sex_char <- ifelse(pheno$sex=="Male", "M","F")
pheno$sex_val = ifelse(pheno$sex=="Male", 1.000, 0.742)

table(pheno$age)
table(pheno$baseline_cr)
pheno$eth = "non-black"
pheno$units = "US"
print("starting ckd-epi")

### use non-black for everyone - workaround for race-free calculations
pheno <- pheno %>%
  mutate(
    baseline_eGFR = ckd_epi_US(creat = baseline_cr, age = age, sex = sex_char, eth = eth),
    baseline_eGFR2 = 175 * ((baseline_cr)**(-1.154)) * ((age)**(-0.203)) * (sex_val)
  )

pheno$baseline_eGFR
### re-compute eGFR - because some values come up with Inf
#pheno_df_2$cr <- ifelse(pheno_df_2$cr==0, pheno_df_2$ValueNumeric, pheno_df_2$cr)
#pheno_df_2$eGFR <- ckd_epi(creat = pheno_df_2$cr, age = pheno_df_2$age, sex = pheno_df_2$sex_char, eth = pheno_df_2$eth, units = pheno_df_2$units)

print("done with ckd-epi")
length(unique(pheno$Subject_ID))

####

### Add aki stage to data frame
#pheno$max_aki_stage <- final_df$max_aki_stage[match(pheno$Subject_ID, final_df$Subject_ID)]
pheno$max_aki_stage <- metadata$max_aki_stage[match(pheno$Subject_ID, metadata$Subject_ID)]
table(pheno$max_aki_stage)
pheno$max_aki_outcome <- metadata$max_aki_outcome[match(pheno$Subject_ID, metadata$Subject_ID)]
table(pheno$max_aki_outcome)
pheno$RNA_Corrected_Sample <- metadata$RNA_Corrected_Sample[match(pheno$Subject_ID, metadata$Subject_ID)]
#pheno_orig = pheno
pheno_uq = pheno %>% drop_na(RNA_Corrected_Sample)
pheno = pheno_uq
#source("code/preprocess_somalogic_data.R")
table(pheno$max_aki_stage)

### Subset to egfr after discharge
pheno <- pheno[pheno$date > (pheno$chosen_somalogic_date + (pheno$hospital_stay- pheno$somalogic_hospital_day)),]
# pheno <- pheno[pheno$Encounter_Discharge_Date]
# subset3 <- test[test$Subject_ID=="PICR3008",]

pheno<- pheno %>% group_by( Subject_ID, date) %>% 
  slice(1) %>%
  as.data.frame()

pheno <- pheno[pheno$eGFR<150,]
# pheno <- pheno[pheno$ckd==0,]

pheno$days_since_somalogic <- as.numeric(pheno$date-  pheno$chosen_somalogic_date)

pheno$days_since_discharge <- as.numeric(pheno$date- (pheno$chosen_somalogic_date + (pheno$hospital_stay- pheno$somalogic_hospital_day)))

summary(pheno$days_since_discharge)
plot(density(pheno$days_since_discharge))

pheno_summary = pheno %>% group_by(Subject_ID) %>% 
  summarize(
    Subject_ID = unique(Subject_ID),
    num_egfr = n(),
    age = unique(age),
    sex = unique(sex),
    
    median_egfr = median(eGFR, na.rm = TRUE),
    max_egfr = max(eGFR, na.rm=TRUE),
    min_egfr = min(eGFR, na.rm = TRUE),
    mean_egfr = mean(eGFR, na.rm = TRUE),
    var_egfr = var(eGFR, na.rm = TRUE),
    max_egfr_date= max(days_since_discharge, na.rm=T),
    min_egfr_date= min(days_since_discharge, na.rm=T),
    median_egfr_date= median(days_since_discharge, na.rm=T),
    history_ckd = unique(ckd),
    history_aki = unique(aki)
)


head(pheno_summary)
dim(pheno_summary)
summary(pheno_summary$max_egfr_date)
table(pheno_summary$num_egfr>1)
table(pheno_summary$max_egfr_date>100)

summary(pheno_summary$num_egfr)
plot(density(pheno_summary$num_egfr))
table(pheno_summary$num_egfr >1000)

ggplot(data= pheno_summary, aes(x= min_egfr_date))+
  geom_density(size=3)+
  
  theme_minimal() +
  ylab("eGFR") +
  xlab("Days since protein measurement")+
  scale_x_continuous(breaks = round(seq(0,550, by = 100),1))+
  # scale_color_manual(values=c( "black","blue", "red")) +
  # labs(color=paste0(protein_name, " level"))+
  scale_color_brewer(palette="Set1")+
  theme(axis.title = element_text(size = 0),
        plot.subtitle=element_text(size=18, hjust=0.5, face="italic", color="black"),
        axis.text = element_text(size = 25, color="black", face="bold"),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 15),
        legend.position = "none",
        title =element_text(size=20, face='bold'),
        # plot.title = element_text(size = 30,hjust = 0.5), 
  )

print("saving pheno summary file")
save(pheno_summary, file="pheno_summary.rData")


#########  ANALYZE  #######################

### Load proteins to plot
#load(file="results/merged_proteomic_mcgill_sinai/all_sinai_proteins_severe_aki_no_interaction_term.RData")
#sinai_results = get(load("/Users/jayaramanp/Dropbox/MSSM/LAB/AKI_COVID_WGCNA/AKI_COVID_ATTEMPTS/AKI_COVID_newAKI123_vs_0_attempt3/RNA_Library_Prep_Plate/Neutrophils/RNASEQ_FIRST_OF_PAIR.PF_INDEL_RATE/Plasma_cells/MECHANICAL_VENTILATION/Sex/Macrophages_M0/Age/T_cells_CD4_memory_activated/ckd/RNASEQ_UNMAPPED_READS/RNASEQ_DEG (56).rData"))
sinai_results = get(load("RNASEQ_DEG (56).rData"))

#sinai_results2 = readxl::read_excel("/Users/jayaramanp/Dropbox/MSSM/LAB/AKI_COVID_WGCNA/AKI_COVID_ATTEMPTS/AKI_COVID_newAKI123_vs_0_attempt3/Analysis_results_DEG56_829pm/Molecules_Significant.xlsx")
sinai_results2 = readxl::read_excel("Molecules_Significant.xlsx")

sinai_results <- sinai_results[sinai_results$adj.P.Val<=0.05,]
replicated_results <- data.frame(sinai_results2[sinai_results2$`Expr False Discovery Rate (q-value)`<=0.05,])
replicated_results  <- replicated_results[!is.na(replicated_results$Symbol),]
rownames(replicated_results) = replicated_results$Ensembl
sinai_results  <- sinai_results[!is.na(sinai_results$gene),]
rownames(sinai_results) = sinai_results$gene
dim(replicated_results)

### Filter to only proteins with effect in same direction between mcgill and sinai
#replicated_results <- replicated_results[(replicated_results$logFC>0 &replicated_results$mcgill_logFC>0) | (replicated_results$logFC<0 &replicated_results$mcgill_logFC<0),]


### Add number of egfr measurements to data frame
pheno$num_egfr <- pheno_summary$num_egfr[match(pheno$Subject_ID, pheno_summary$Subject_ID)]
table(pheno$num_egfr)
#### Associate AKI with egfr
table(pheno$max_aki_outcome)
table(pheno$composite_outcome)


### get protein names
#pheno = pheno[,1:40]
pheno$Subject_ID
#vars <- colnames(proteomic)[-1]
#vars <-colnames(pheno)[colnames(pheno) %in% vars]
rownames(sinai_results)
vars <- rownames(sinai_results)
vars <- vars[vars %in% rownames(sinai_results)]
vars

# pheno <- pheno %>%
#   mutate(
#     baseline_eGFR = ckd_epi(creat = baseline_cr, age = age, sex = sex_char, eth = eth, units = units),
#     baseline_eGFR2 = 175 * ((baseline_cr)**(-1.154)) * ((age)**(-0.203)) * (sex_val)
#   )

common_columns = intersect(colnames(pheno), colnames(metadata))
common_columns

pheno_meta = merge(pheno, metadata, by = common_columns)
table(pheno_meta$RNASEQ_UNMAPPED_READS)

print("saving pheno_meta rData")
save(pheno_meta, file="pheno_meta.rData")


quit(status=1)





### Loop over genes:
print("starting parallel process")

#run_model_get_lmme <- function(pheno_meta, gene){
#  gene = eval(parse(text = paste("log2(as.numeric(", str(gene), ")+1)")))
#  print(gene)
#  #model <- lmer(eGFR-baseline_eGFR ~ num_egfr+days_since_somalogic+baseline_cr+max_aki_stage+age+sex+log2(as.numeric(pheno[[gene_expr]])+1) + (1|Subject_ID), data = pheno)
#  #model <- lmer(eGFR - baseline_eGFR ~ age + baseline_cr + days_since_discharge + num_egfr + max_aki_outcome + ckd + icu+ RIN + (1|Subject_ID) + log2(as.numeric(gene_expr)+1), data = pheno_meta)
#  #model <- lmer(eGFR - baseline_eGFR ~ age + baseline_cr + days_since_discharge + num_egfr + max_aki_outcome + ckd + icu + log2(as.numeric(gene)+1) + (1|Subject_ID), data = pheno_meta)
#  model <- lmer(eGFR - baseline_eGFR ~ age + baseline_cr + days_since_discharge + num_egfr + max_aki_outcome + ckd + icu + gene + (1|Subject_ID), data = pheno_meta)
#  print(summary(model))
#  return(model)
#}

### loop over proteins
#results_df2 <- data.frame(var = vars, p =NA, coef= NA) 
#vars[2635]
intermediate_directory <- './tempdir_072723'

# #vars = vars[1:10]
# table(pheno_meta$RNA_DV200_Percent)
# #model <- run_model_get_lmme(pheno_meta=pheno_meta, gene=str(vars[1]))
# sampleindex = sample(1:2635, 5, replace=FALSE)
# sampleindex
# vars2 = vars[sampleindex]
# vars2
# for (i in 1:length(vars2)) {
#   gene = vars2[i]
#   print(gene)
#   #lmer_func = eval(parse(text = paste("eGFR - baseline_eGFR ~ age + baseline_cr + days_since_discharge + num_egfr + max_aki_outcome + ckd + icu + log2(as.numeric(", gene, ")+1) + (1|Subject_ID)")))
#   lmer_func = eval(parse(text = paste("eGFR - baseline_eGFR ~ log2(as.numeric(", gene, ")+1) + (1|Subject_ID)")))
#   lmer_func
#   lmme_model = lmer(lmer_func, data = pheno_meta)
#   summary = summary(lmme_model)
#   print(summary)
#   ano = anova(lmme_model)
#   print(ano)
# }
#print(summary$coefficients[9,])
#print(summary$coefficients[9,5])
#pval = summary$coefficients[nrow(summary$coefficients),5]
#print(pval)
#padj = p.adjust(pval, method = "bonferroni")
#coef(model)
#getME(model, "fixef")
#confint(model)

print("setting cores")

#cores_2=4
cores_2 = 20
cores_2
cl <- makeCluster(cores_2[1]-1) #not to overload your computer
clusterExport(cl, c("vars", "pheno_meta"))
clusterEvalQ(cl, require(c(lme4, lmerTest)))
registerDoParallel(cl)

print("registered and running clusterCall")

#clusterCall(cl, function(run_model_get_lmme) .libPaths("/Users/jayaramanp/workspace/R/lib/4.04/"), .libPaths())
clusterCall(cl, function(run_model_get_lmme) .libPaths("/hpc/users/jayarp02/R/lib/4.0.3/"), .libPaths())

print("running foreach!")
#i=1:length(vars)
final_results_df2 <- data.frame(foreach(i=1:length(vars), .packages=c("lme4", "lmerTest"), .combine=rbind) %dopar% {
  each_filename <- paste0('DEGMATRIX_TEST_truth_', as.character(vars[i]),'_' , as.character(i), '.rdata') 
  each_filepath <- file.path(intermediate_directory, each_filename)
  #lmme_model <- run_model_get_lmme(pheno_meta = pheno_meta, gene_expr = gene)
  #lmme_model <- lmer(eGFR - baseline_eGFR ~ age + baseline_cr + days_since_discharge + num_egfr + max_aki_outcome + ckd + icu + eval(parse(text = paste("log2(as.numeric(", vars[i], ")+1)")))  + (1|Subject_ID), data = pheno_meta)
  #lmer_func = eval(parse(text = paste("eGFR - baseline_eGFR ~ age + baseline_cr + days_since_discharge + num_egfr + max_aki_outcome + ckd + icu + log2(as.numeric(", vars[i], ")+1) + (1|Subject_ID)")))
  lmer_func = eval(parse(text = paste("eGFR ~ log2(as.numeric(", vars[i], ")+1)*(as.numeric(days_since_somalogic)) + (1|Subject_ID) ")))
  #lmer_func
  lmme_model = lmer(lmer_func, data = pheno_meta)
  summary = summary(lmme_model)
  summary$coefficients
  p_val <- summary$coefficients[nrow(summary$coefficients), 5]
  #coef_val <- summary$coefficients[nrow(summary$coefficients), 1]
  #results_subset_df = data.frame(var = gene, p = p_val, coef = coeff)
  #results_subset_df = c(vars[i], p_val, coef_val)
  
  ##### save model if the interaction term is significant 
  #if (p_val <= 0.05){
      #save(lmme_model, file = each_filepath)
  #}
  
  #return(results_subset_df)
  results_subset_df = summary$coefficients
  return(results_subset_df)
})

print("parallel rbind completed!")

save(final_results_df2, file = "./tempdir_072723/Final_results_summaryCoeff_F1_SomalogicOnly_072623.rData")

stopCluster(cl)
print("done with cluster")
head(final_results_df2)
#stopCluster(cl_2)

quit(status=1)

# ############################# TESTING WITHOUT PARALLELIZATION##################
# for(i in 1:nrow(results_df2)){
#   # model <- lmer(eGFR ~ log2(as.numeric(pheno$`Tankyrase-1`)) +age+sex+ckd+ (1|Subject_ID), data = pheno)
#   print(i)
# 
# ### with random intercept
# #  model <- lmer(eGFR-baseline_egfr ~num_egfr+days_since_somalogic+baseline_cr+max_aki_stage+ age+sex+days_since_somalogic*x) + (1|Subject_ID), data = pheno)
# ##  #model <- lmer(eGFR-baseline_egfr ~num_egfr+days_since_somalogic+baseline_cr+max_aki_stage+ age+sex+log2(as.numeric(pheno[[results_df$var[i]]])) + (1|Subject_ID), data = pheno)
# ## model <- lmer(log2(as.numeric(merged[merged$aki==1,results_df$var[i]])) ~ modified_sofa+age+sex+ckd+type+ (1|Subject_ID), data = merged[merged$aki==1,])
# # model <- lm(log2(as.numeric(merged[[results_df$var[i]]])) ~ somalogic_hospital_day+type*aki, data = merged)
#   #model <- lmer(eGFR-baseline_eGFR ~ num_egfr+days_since_somalogic+baseline_cr+max_aki_outcome+age+sex+log2(as.numeric(pheno[[results_df2$var[i]]])+1) + (1|Subject_ID), data = pheno)
#   #model <- lmer(eGFR-baseline_eGFR ~ num_egfr+days_since_somalogic+baseline_cr+max_aki_stage+age+sex+log2(as.numeric(pheno[[results_df2$var[i]]])+1) + (1|Subject_ID), data = pheno)
#   model <- lmer(eGFR - baseline_eGFR ~ age + baseline_cr + days_since_discharge + num_egfr + max_aki_outcome + ckd + icu + (1|Subject_ID) + log2(as.numeric(pheno[[results_df2$var[i]]])+1), data = pheno_meta)
# 
#   summary <- summary(model)
#   # summary$coefficients
#   results_df2$p[i] <- summary$coefficients[nrow(summary$coefficients), 5]
#   results_df2$coef[i] <- summary$coefficients[nrow(summary$coefficients), 1]
# }
# head(results_df2)
################################ENDOF TESTING ###############################

final_results_df2_genes = final_results_df2 %>% filter(grepl("log2", rownames(final_results_df2)))
final_results_df2_genes$fdr <- p.adjust(final_results_df2_genes$Pr...t.., method = 'BH')
final_results_df2_genes <- final_results_df2_genes[order(final_results_df2_genes$Pr...t.., decreasing = F),]
head(final_results_df2_genes)

results_df2$EntrezGeneSymbol <- sinai_results2$Symbol[match(results_df2$var, sinai_results2$Ensembl)]
table(results_df2$fdr<0.05)
results_df2$adj.P <- p.adjust(results_df2$p, method = "bonferroni")
table(results_df2$adj.P<0.05)
#results_df$gene_name <- protein_metadata$EntrezGeneSymbol[match(results_df$var, protein_metadata$TargetFullName)]
save(results_df2, file="egfr_association_all_proteins.RData")
sig_results <- results_df[results_df$adj.P<0.05,]
sig_results <- sig_results[order(sig_results$coef, decreasing = F),]
write.csv(sig_results, file="egfr_associated.csv")

### Plot
protein_name <- results_df$var[dim(results_df)[1]-1]
protein_name <- results_df$var[1]

cutoff <- quantile(log2(final_df_transposed_rnaseq[[protein_name]]), na.rm = T, probs = c(.33,.67,1))
cutoff <- c(-Inf, cutoff)
# pheno$protein_high <- ifelse(log2(pheno[[protein_name]])>cutoff, "High protein ", "Low protein")
pheno$protein_cat <- cut(log2(pheno[[protein_name]]), cutoff, labels=c("low","intermediate","high"))
summary(pheno$protein_cat)
# p <- ggplot(data = pheno, aes(x = date, y = eGFR, group = Subject_ID, col = Subject_ID))
# p + geom_line()


# Trefoil factor 3

# plot_ly( x = pheno$days_since_somalogic, y= pheno$eGFR,group=pheno$Subject_ID,
#        mode = 'lines', color=~pheno$protein_cat)
# %>%
#   add_surface()

ggplot(data = pheno, aes(x = days_since_somalogic, y =eGFR, group=protein_cat, color= factor(protein_cat)))+
  # geom_line()+
  stat_smooth()+
  # facet_wrap(~protein_high, nrow = 1)
  # scale_color_manual(labels = c(0,1), values = c("blue", "red")) +
  
  theme_minimal() +
  ylab("eGFR") +
  xlab("Days since protein measurement")+
  # scale_color_manual(values=c( "black","blue", "red")) +
  labs(color=paste0(protein_name, " level"))+
  scale_color_brewer(palette="Set1")+
  theme(axis.title = element_text(size = 0),
        plot.subtitle=element_text(size=18, hjust=0.5, face="italic", color="black"),
        axis.text = element_text(size = 35, color="black", face="bold"),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 15),
        legend.position = "none",
        title =element_text(size=20, face='bold'),
        # plot.title = element_text(size = 30,hjust = 0.5), 
  )

# plot_trajectories(data = pheno,
#                   id_var = "Subject_ID", 
#                   var_list =c("eGFR") ,
#                   xlab = "days_since_somalogic", ylab = "eGFR",
#                   connect_missing = F, 
#                   random_sample_frac = 0.018, 
#                   title_n = T) +
#   facet_wrap(~Subject_ID)
#> Warning: Removed 1 row(s) containing missing values

# ggplot(pheno,aes(x=days_since_somalogic,y=Subject_ID,fill=eGFR)) + geom_tile(colour='black') 

######## 
##### Plot egfr vs AKi associations
########

combined_results <- results_df
combined_results$aki_sinai_logFC <- sinai_results$logFC[match(combined_results$var, rownames(sinai_results))]
#combined_results$aki_mcgill_logFC <- replicated_results$mcgill_logFC[match(combined_results$var, rownames(replicated_results))]

combined_results$sig_egfr <- ifelse(combined_results$adj.P<0.05,1,0)

combined_results$to_label <- 0
combined_results$to_label[combined_results$sig_egfr==1 &combined_results$coef< -10] <- 1

combined_results$delabel <- NA
# combined_df$delabel[combined_df$diffexpressed != "NO"] <-combined_df$gene_name[combined_df$diffexpressed != "NO"]
combined_results$delabel[combined_results$to_label ==1] <- combined_results$EntrezGeneSymbol[combined_results$to_label ==1] 



ggplot(data=combined_results , aes(x = aki_sinai_logFC, y =coef, color= factor(sig_egfr), label = delabel))+
  # stat_smooth()+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=4)+
  geom_label_repel( size=7, max.overlaps = 30, nudge_y= .2) +
  # facet_wrap(~protein_high, nrow = 1)
  # scale_color_manual(labels = c(0,1), values = c("blue", "red")) +
  
  theme_minimal() +
  ylab("eGFR") +
  xlab("Days since protein measurement")+
    
  # scale_color_manual(values=c( "black","blue", "red")) +
  # labs(color=paste0(protein_name, " level"))+
  scale_color_brewer(palette="Set1")+
  theme(axis.title = element_text(size = 0),
        plot.subtitle=element_text(size=18, hjust=0.5, face="italic", color="black"),
        axis.text = element_text(size = 35, color= "black", face= "bold"),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 15),
        legend.position = "none",
        title =element_text(size=20, face='bold'),
        # plot.title = element_text(size = 30,hjust = 0.5), 
  )
