#module load R/4.0.3 #very important that you load this R version, dream will not work with the regular version that auto loads on minerva
#R
#############
### final code to run the full analysis from the metadata to the cell type deconvolution to the differential gene expression. Every time the analysis is run, it will create a new folder with all analysis results within the foler. 
##############
rm(list=ls())
options(stringsAsFactors=F)
##################################################
#load libraries
##################################################

.libPaths( c( .libPaths(), "/hpc/users/jayarp02/R/lib/4.0.3/") ) #  need to do this before loading any packages!!!

library(assertive)
#library(dplyr)
#library(tidyr)
library(tidyverse)
library(edgeR)
library(BiocParallel)
library(batchtools)
library(devtools)
library(withr)
library(limma)
library(Glimma)
library(variancePartition)#, lib.loc="~/.Rlib")
library(ggplot2)
library(gridExtra)
library(grid)
library(doParallel)
registerDoParallel(20)
library(sp)
#library(CellMix)
library(biomaRt)
library(gsubfn)
library(data.table)
library(sp)
library(Matrix)

#############################################################################################################
#Run the functions first
#############################################################################################################
###IMPORTANT
multiplot_same_legend <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  #get_legend_info
  mylegend<-g_legend(plots[[1]])

  numPlots = length(plots)

  for(i in 1:numPlots){
    plots[[i]] <- plots[[i]] + theme(legend.position="none")
  }
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
    if(is.null(mylegend)==F){
      grid.draw(mylegend)
    }    
  }
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if(length(leg)>0){
    legend <- tmp$grobs[[leg]]
  }else{
    legend <- c()
  }
  return(legend)}



#options(width=150)

###IMPORTANT - correlation between gene expression and covariates. this is where you choose covariates. you will choose covariates that are highly correlated wth PCs. 
###
canCorAllAgainstAll_Original <- function(X, Y = X,minimum_intersect=0) {
    # Compute canonical correlation of all columns of X against all columns of Y,
    # similar to variancePartition::canCorPairs.
    print("colnames X")
    colnames(X)

    print("colnames Y")
    colnames(Y)

    print("print X_formulas")
    library(stringr)
    X_formulas <- lapply(str_c("~", colnames(X)), as.formula)
    print("print X_varList")
    X_varList <- lapply(X_formulas, function(xf) model.matrix.lm(xf, X, na.action = "na.pass")[,-1, drop = FALSE])
    head(X_varList)
    print("print Y_formulas")
    Y_formulas <- lapply(str_c("~", colnames(Y)), as.formula)
    print(Y_formulas)
    print("Y_varList")
    Y_varList <- lapply(Y_formulas, function(yf) model.matrix.lm(yf, Y, na.action = "na.pass")[,-1, drop = FALSE])
    head(Y_varList)

    print("XY cc")
    XY_cc <- matrix(nrow = ncol(X), ncol = ncol(Y), data = 0,
                    dimnames = list(colnames(X), colnames(Y)))
    head(XY_cc)
    print("X_varList:")
    print(X_varList)
    print("Y_varList:")
    print(Y_varList)
    print("for loop!")
    for (ix in seq_along(X_varList)) {
        print("ix:")
        print(ix)
	print("looping ix in seq_along X_varList")
        keep1 = apply(X_varList[[ix]], 1, function(x) !any(is.na(x)))
	print("keep1")
	head(keep1)
        for (iy in seq_along(Y_varList)) {
	    print("iy:")
            print(iy)
	    print("looping iy in seq_along Y_varList")
            keep2 = apply(Y_varList[[iy]], 1, function(x) !any(is.na(x)))
	    print("getting keep")
            keep = keep1 & keep2
	    head(keep)
            if(sum(keep)>minimum_intersect){
	    print("if sum keep > minimum_intersect")
	    print("fit")
            fit <- cancor(X_varList[[ix]][keep, , drop = FALSE], Y_varList[[iy]][keep, , drop = FALSE])
            head(fit)

            print("XY_cc[ix,iy]")
            # Using root-mean-square to summarize, as discussed with Gabriel Hoffman
	    XY_cc[ix,iy] <- sqrt(mean(fit$cor^2))
            head(XY_cc)
        }else{
	    print("or else NA")
            XY_cc[ix,iy] <- NA
        }

        print("done with iy in seq_along(Y_varList)")
        }
    }

    print("return XY_cc")
    return(XY_cc)
    }

#####################################################. MAIN. ########################################################
print("creating main analysis directory")
results_dir = paste("AAFINAL_AKI123_vs_0", format(Sys.time(), "%d-%b-%Y-%H-%M-%S"), sep = "_")
dir.create(results_dir)
results_dir = paste0(results_dir,"/")
results_dir

cat("README",file=paste(results_dir, "README.txt", sep=""),sep="\n")

lbpcov <- readRDS("rnaseq_final_merged.rds")

##remove samples that do not have cell fraction counts at this step because removing them reduces factor levels!
print("removing samples that arent in cell fraction counts due to bad qc...")
###### samples discarded when cell fracyions counts were adjusted for
outliersamples = c('PICR3019T12_Plate_9','PICR5016T1_Plate_11','PICR5033T1_Plate_9','PICR5047T1_Plate_5','PICR5067T1_Plate_9','PICR5082T1_Plate_11','PICR5090T1_Plate_9','PICR5094T8_Plate_14','PICR7023T12_Plate_9','PICR7059T13_Plate_13','PICR7220T4_Plate_2','PICR7236T8_Plate_4','PICR7263T4_Plate_18','PICR7306T1_Plate_7','PICR7317T12_Plate_9','PICR7366T1_Plate_3','PICR7392T4_Plate_3','PICR7428T12_Plate_17','PICR7433T12_Plate_15','PICR7446T12_Plate_19','PICR7469T1_Plate_13','PICR7472T1_Plate_8','PICR7480T8_Plate_15','PICR7496T8_Plate_14','PICR7501T12_Plate_14','PICR7508T1_Plate_11')

#### remove samples not in cell fraction counts. They have poor quality
print("removing samples not in cell fraction counts")

lbpcov = lbpcov[!lbpcov$RNA_Corrected_Sample %in% outliersamples,]

#### ADDITIONAL STEP add cell fraction file to file!

## read in cell fraction counts
cellfrac <- read.delim("/sc/arion/projects/mscic1/results/Noam/cibersortx/outdir_LM22_MainCovid/CIBERSORTx_Results.txt", sep="\t", header=T, stringsAsFactors=F)

#replace "." in column names with "_"
colnames(cellfrac) <- gsub("\\.", "_", colnames(cellfrac))

print("cellfrac structure")

print("factor cellfrac")
str(cellfrac, list.len=Inf)
colnames(cellfrac)

#'B.cells.naive','B.cells.memory','Plasma.cells','T.cells.CD8','T.cells.CD4.naive','T.cells.CD4.memory.resting','T.cells.CD4.memory.activated','T.cells.follicular.helper','T.cells.gamma.delta','T.cells.regulatory..Tregs.','NK.cells.resting','NK.cells.activated','Monocytes','Macrophages.M0','Macrophages.M1','Macrophages.M2','Dendritic.cells.resting','Dendritic.cells.activated','Mast.cells.resting','Mast.cells.activated','Eosinophils','Neutrophils','P.value','Correlation','RMSE'

cellfrac = cellfrac[,names(cellfrac) %in% c('Mixture','B_cells_naive','B_cells_memory','Plasma_cells','T_cells_CD8','T_cells_CD4_naive','T_cells_CD4_memory_resting','T_cells_CD4_memory_activated','T_cells_follicular_helper','T_cells_gamma_delta','T_cells_regulatory__Tregs_','NK_cells_resting','NK_cells_activated','Monocytes','Macrophages_M0','Macrophages_M1','Macrophages_M2','Dendritic_cells_activated','Dendritic_cells_resting','Mast_cells_resting','Eosinophils','Neutrophils')]
#cellfrac

cellfrac[sapply(cellfrac, is.character)] <- lapply(cellfrac[sapply(cellfrac, is.character)], as.factor)

print("now factor")
str(cellfrac, list.len=Inf)

newlbpcov = merge(x=lbpcov, y=cellfrac, by.x=c('RNA_Corrected_Sample'), by.y=c('Mixture'), sort = TRUE)

head(newlbpcov)
nrow(newlbpcov)


### reassign to lbpcov with cellfraction
lbpcov = newlbpcov

#names_new = colnames(lbpcov)
names_new = c(
  "eskd",
  "Specimen_Stool_Collected",
  "Specimen_Saliva_Collected",
  "Specimen_BAL_Collected",
  "Specimen_SomaLogics_Performed",
  "ECMO_Flag",
  "MED_CYCLOPHOSPHAMIDE",
  "MED_RITUXIMAB",
  "MED_HEPARIN",
  "MED_AZITHROMYCIN",
  "MED_HYDROXYCHLOROQUINE",
  "MED_APIXABAN",
  "MED_NOREPINEPHRINE_LEVOPHED",
  "MED_PHENYLEPHRINE",
  "MED_METHYLPREDNISOLONE",
  "MED_VASOPRESSIN_VASOSTRICT",
  "MED_TISSUE_PLASMINOGEN_ACTIVATOR",
  "MED_RIVAROXABAN",
  "MED_WARFARIN_COUMADIN",
  "MED_TOCILIZUMAB",
  "MED_ECULIZUMAB",
  "MED_DOBUTAMINE",
  "MED_GIMISILUMAB",
  "MED_VEDOLIZUMAB_ENTYVIO",
  "MED_NO_ppm",
  "MED_MILRINONE",
  "MED_EPINEPHRINE",
  "MED_DOPAMINE",
  "MED_CISPLATIN_PLATINOL",
  "MED_PEMBROLIZUMAB_KEYTRUDA",
  "MED_ANTI.THYMOCYTE_GLOBULIN",
  "MED_HYDROCORTISONE",
  "MED_SARILUMAB",
  "MED_CARBOPLATIN",
  "MED_METHOTREXATE",
  "MED_DOXORUBICIN",
  "MED_INFLIXIMAB_REMICADE",
  "MED_PREDNISOLONE",
  "MED_ASRPARGINASE",
  "MED_NIVOLUMAB_OPDIVO",
  "MED_EDOXABAN",
  "MED_USTEKINUMAB_STELARA",
  "MED_ARGATROBAN",
  "COVID19_Prior_Infection_Reported",
  "COVID19_Vaccine_Received",
  "COVID19_Vaccine_Clinical_Trial_Flag",
  "RNA_Mislabel_Flag",
  "RNA_Batch_Control_Flag",
  "Subject_ID_Retired",
  "COVID19_No_Test_Result",
  "aki.y",
  "Consent_Status",
  "Specimen_ELISA_Grand_Serology",
  "Encounter_Type",
  "Encounter_Patient_Class",
  "COVID19_Vaccine_Type",
  "Blood_Mislabel_Error_Type",
  "RNA_Uncorrected_Sample",
  "RNA_Uncorrected_Blood_Sample",
  "RNA_Personal_Stock_and_Submission_Plate_Wells",
  "RNA_Library_Prep_Plate_Well",
  "RNA_Uncorrected_Subject_ID",
  "RNA_Mislabel_Error_Type",
  "Uncorrected_Blood_Sample",
  "Subject_ID.x",
  "max_aki_stage",
  "aki_last_date",
  "aki_hospital_day",
  "somalogic_hospital_day",
  "severity",
  "icu",
  "death",
  "Subject_ID.y",
  "Time_Point",
  "Sex",
  "Race_From_Consent",
  "Ethnicity_From_Consent",
  "Race_From_EHR",
  "Ethnicity_From_EHR",
  "Specimen_SST_Tube_Collected",
  "Specimen_CPT_Tube_Collected",
  "Specimen_Tempus_Tube_Collected",
  "Specimen_CyTOF_Performed",
  "Specimen_CyTOF_Sample_QC",
  "Specimen_Olink_Performed",
  "Specimen_WGS_Performed",
  "Specimen_RNAseq_Performed",
  "Convalescent_Plasma",
  "COVID19_Order_Result",
  "RRT_Flag",
  "COVID19_Antibody_Assay_Detected",
  "COVID19_Quantitative_Antibodies",
  "O2_Therapy",
  "VENTILATION_TYPE",
  "MECHANICAL_VENTILATION",
  "ICU_Status",
  "Viral_Load_Run_Batch",
  "Encounter_Discharge_Location",
  "Encounter_ICU_Flag",
  "ENCOUNTER_COMORBID_ARDS",
  "ENCOUNTER_COMORBID_ACUTE_KIDNEY_INJURY",
  "ENCOUNTER_COMORBID_ACUTE_VENOUS_THROMBOEMBOLISM",
  "ENCOUNTER_COMORBID_CEREBRAL_INFARCTION",
  "ENCOUNTER_COMORBID_INTRACEREBRAL_HEMORRHAGE",
  "ENCOUNTER_COMORBID_ACUTE_MI",
  "NURSING_CAC_IP_R_NON_SURGICAL_WOUND_PRESENT",
  "NURSING_CAC_IP_R_PRESSURE_ULCER_PRESENT_ON_HOSP_ADM",
  "NURSING_CAC_IP_R_PRESSURE_ULCER_UNIT_TRANSFER",
  "NURSING_CAC_R_BEDSIDE_CARDIAC_MONITOR_ON",
  "NURSING_CAC_R_LEVEL_OF_CONSCIOUSNESS_MULTISELECT",
  "NURSING_CAC_R_PVS_EDEMA",
  "NURSING_CAC_R_ABDOMEN_TENDERNESS",
  "NURSING_CAC_R_BOWEL_SOUNDS_GENERAL_ALL_4_QUADRANTS",
  "NURSING_CAC_R_RESPIRATORY_PATTERN",
  "NURSING_CAC_R_TELEMETRY_CARDIAC_MONITOR",
  "NURSING_CAC_R_BREATH_SOUNDS_BILATERAL",
  "MED_REMDESIVIR",
  "MED_DEXAMETHASONE",
  "MED_ENOXAPARIN",
  "MED_FAMOTIDINE",
  "MED_PREDNISONE",
  "SMOKING_STATUS",
  "BLOOD_TYPE",
  "DECEASED_INDICATOR",
  "Cytomegalovirus",
  "Herpes_Simplex",
  "Trial_Tocilizumab",
  "Trial_Remdesivir",
  "Trial_Hydroxychloroquine",
  "Trial_Anakinra",
  "Trial_Azithromycin",
  "COMORBID_ASTHMA",
  "COMORBID_COPD",
  "COMORBID_HTN",
  "COMORBID_OBSTRUCTIVE_SLEEP_APNEA",
  "COMORBID_DIABETES",
  "COMORBID_CHRONIC_KIDNEY_DISEASE",
  "COMORBID_CANCER_FLAG",
  "COMORBID_CANCER_DIAGNOSIS_DESCRIPTION",
  "COMORBID_CORONARY_ARTERY_DISEASE",
  "COMORBID_ATRIAL_FIBRILLATION",
  "COMORBID_HEART_FAILURE",
  "COMORBID_CHRONIC_VIRAL_HEPATITIS",
  "COMORBID_ALCOHOLIC_NONALCOHOLIC_LIVER_DISEASE",
  "COMORBID_CROHNS_DISEASE",
  "COMORBID_ULCERATIVE_COLITIS",
  "Immunization_Hepatitis_A_Adult_Flag",
  "Immunization_Hepatitis_B_Adult_Flag",
  "Immunization_Influenza_Flag",
  "Immunization_Meningococcal_Adult_Flag",
  "Immunization_Meningococcal_Pediatric_Flag",
  "Immunization_Pneumococcal_Adult_Flag",
  "Immunization_Pneumococcal_Pediatric_Flag",
  "Immunization_Shingles_Adult_Flag",
  "Immunization_Tetanus_Adult_Flag",
  "Immunization_Tetanus_Pediatric_Flag",
  "COVID19_Vaccine_Doses",
  "Corrected_Blood_Sample",
  "Uncorrected_Subject_ID",
  "Uncorrected_Time_Point",
  "Blood_Mislabel_Flag",
  "Serology_Any_Ig_Detected",
  "Serology_Any_IgAGM_Detected",
  "COVID19_Within_Encounter",
  "COVID19_Days_Since_Encounter_Start",
  "COVID19_Days_Until_Encounter_End",
  "RNA_Technician",
  "RNA_Library_Prep_Plate",
  "RNA_Extraction_Batch",
  "RNA_Uncorrected_Time_Point",
  "RNASEQ_Run_Sample",
  "Post_COVID19_Symptom_Anxiety_Depression_Ever",
  "Post_COVID19_Symptom_Cavities_Teeth_Problems_Ever",
  "Post_COVID19_Symptom_Chest_Pain_Cardiac_Issues_Ever",
  "Post_COVID19_Symptom_Eating_More_Less_Ever",
  "Post_COVID19_Symptom_Hair_Loss_Ever",
  "Post_COVID19_Symptom_Headaches_Ever",
  "Post_COVID19_Symptom_Increased_Mucus_Ever",
  "Post_COVID19_Symptom_Joint_Pain_Ever",
  "Post_COVID19_Symptom_Liver_Issues_Ever",
  "Post_COVID19_Symptom_Lung_Problems_Ever",
  "Post_COVID19_Symptom_Memory_Thought_Problems_Ever",
  "Post_COVID19_Symptom_Muscle_Pain_Ever",
  "Post_COVID19_Symptom_Nausea_Diarrhea_Vomiting_Ever",
  "Post_COVID19_Symptom_Need_Supplemental_O2_Ever",
  "Post_COVID19_Symptom_Pneumonia_Ever",
  "Post_COVID19_Symptom_Shortness_Of_Breath_Ever",
  "Post_COVID19_Symptom_Skin_Rash_Ever",
  "Post_COVID19_Symptom_Sleep_Problems_Ever",
  "Post_COVID19_Symptom_Smell_Taste_Problems_Ever",
  "Post_COVID19_Symptom_Sore_Throat_Ever",
  "Post_COVID19_Symptom_Toe_Problems_Ever",
  "Post_COVID19_Symptom_Weakness_Or_Fatigue_Ever",
  "Post_COVID19_Symptom_Other_Ever",
  "Post_COVID19_Symptom_None_Ever",
  "Post_COVID19_Symptom_Medication_Flag_Ever",
  "Post_COVID19_New_Fatigue_Ever",
  "Post_COVID19_New_Fatigue_Headache_Ever",
  "Post_COVID19_New_Fatigue_Joint_Pain_Ever",
  "Post_COVID19_New_Fatigue_Muscle_Pain_Ever",
  "Post_COVID19_New_Fatigue_Problems_Concentration_Remembering_Ever",
  "Post_COVID19_New_Fatigue_Sore_Neck_Ever",
  "Post_COVID19_New_Fatigue_Sore_Throat_Ever",
  "Post_COVID19_New_Fatigue_Trouble_Sleeping_Ever",
  "Post_COVID19_New_Fatigue_Trouble_Walking_Ever",
  "Post_COVID19_New_Fatigue_Other_Ever",
  "Post_COVID19_New_Fatigue_None_Ever",
  "Post_COVID19_Health_Ever",
  "Post_COVID19_Quality_Of_Life_Ever",
  "Post_COVID19_Daily_Mood_Ever",
  "COVID19_Ever_Positive",
  "COVID19_Most_Recent_Order_Result",
  "COVID19_Days_Since_Most_Recent_Order",
  "COVID19_Any_Antibody_Detected",
  "COVID19_Most_Recent_Antibody_Result",
  "COVID19_Days_Since_Most_Recent_Antibody",
  "COVID19_Non_Consecutive_Positive_Tests",
  "COVID19_Antibody_Loss",
  "SARSCoV2_Recovery_Status",
  "COVID19_Symptoms_At_Presentation_Reviewed",
  "COVID19_Symptoms_Reported_Onset_Is_Exact",
  "COVID19_Deceased_Within_Encounter",
  "Patient_Classification_At_First_Sample",
  "Severity_Ignoring_COVID19_Status",
  "Severity",
  "Severity_Ambiguous_EOD",
  "WHO_Ordinal_Ignoring_COVID19_Status",
  "WHO_Ordinal",
  "aki_rule1",
  "aki_rule2",
  "aki_stage_2",
  "aki_stage_3",
  "aki_stage"
)
lbpcov[, names_new] <- lapply(lbpcov[, names_new] , factor)
#str(lbpcov)
lbpcov[, names_new] <- lbpcov[,names_new] %>% mutate_all(funs(factor(replace(as.character(.), is.na(.), "NA"))))


print("categorical covariates - factors with levels >=2")
toremove = NULL

for (i in 1:ncol (lbpcov)) 
	{
	if (class (lbpcov[,i]) == "factor") 
		{
		print(nlevels(lbpcov[,i]))
		if( nlevels(lbpcov[,i]) < 2 )
			{
			print(names(lbpcov)[i])
			toremove = append(toremove, names(lbpcov)[i])
			}
		}
	}

print("removing columns names with zero variance")
toremove = append(toremove, "Specimen_RNAseq_Performed")

print("removing columns starting with ICD that we do not need")
notneededicd = c(
  "delta_2_day",
  "seven_day_minimum",
  "Serology_Average_Date_Tested",
  "Serology_Any_IgAGM_Detected",
  "Serology_Highest_Titer_AGM_Only",
  "RNA_Num_Samples_For_Blood_Sample",
  "Specimen_Tempus_Tube_Collected",
  "Time_Point",
  "Uncorrected_Subject_ID",
  "Event_Date",
  "Corrected_Blood_Sample",
  "somalogic_hospital_day",
  "aki_date",
  "aki_date_minus_somalogic_date",
  "aki_hospital_day",
  "hospital_stay",
  "two_day_minimum",
  "LAB_CREATINE_KINASE_MB_StDev",
  "RNA_Elution_Volume_ul",
  "ICU_Hours_Until_Discharge",
  "ICU_Hours_Until_Entry",
  "PaO2_Mean",
  "MED_DAYS_SINCE_OUTPATIENT_ORDER_TISSUE_PLASMINOGEN_ACTIVATOR",
  "MED_DAYS_SINCE_OUTPATIENT_ORDER_GIMISILUMAB",
  "Viral_Load_Run_Date",
  "Viral_Load_Run_Batch",
  "LAB_BASE_EXCESS_VENOUS_StDev",
  "LAB_BRAIN_NATRIURETIC_PROTEIN_Mean",
  "LAB_C_REACTIVE_PROTEIN_StDev",
  "LAB_D_DIMER_StDev",
  "LAB_ESR_Mean",
  "LAB_FERRITIN_StDev",
  "LAB_FIBRINOGEN_StDev",
  "LAB_INR_StDev",
  "aki_rule1",
  "aki_rule2",
  "aki_stage_2",
  "aki_stage_3",
  "somalogic_date",
  "baseline_cr",
  "Subject_ID.x",
  "Subject_ID.y",
  "Uncorrected_Blood_Sample",
  "LAB_PROCALCITONIN_StDev",
  "LAB_URIC_ACID_Num_Observations",
  "ICD10_C85.10",
  "ICD10_A41.9",
  "ICD10_R65.20",
  "ICD10_C83.38",
  "ICD10_D70.9",
  "ICD10_R50.81",
  "ICD10_C83.39",
  "ICD10_R07.89",
  "ICD10_I48.91",
  "ICD10_R79.89",
  "ICD10_I48.0",
  "ICD10_Z79.4",
  "ICD10_I10",
  "ICD10_I21.4",
  "ICD10_I25.10",
  "ICD10_E78.2",
  "ICD10_E11.9",
  "ICD10_Z98.61",
  "ICD10_R09.02",
  "ICD10_R06.89",
  "ICD10_R06.00",
  "ICD10_R68.89",
  "ICD10_I25.9",
  "ICD10_R07.9",
  "ICD10_K29.00",
  "ICD10_J80",
  "ICD10_J18.9",
  "ICD10_J96.01",
  "ICD10_L03.119",
  "ICD9_250.00",
  "ICD10_L02.419",
  "ICD10_J45.909",
  "ICD10_I48.92",
  "ICD10_I95.9",
  "ICD10_Z00.00",
  "ICD10_E66.9",
  "ICD10_M25.561",
  "ICD10_S01.81XA",
  "ICD10_L03.116",
  "ICD10_R06.02",
  "ICD10_E66.01",
  "ICD10_L03.115",
  "ICD10_L03.90",
  "ICD10_R05",
  "ICD10_Z76.0",
  "ICD10_I50.22",
  "ICD10_J06.9",
  "ICD10_K08.89",
  "ICD10_I48.19",
  "ICD10_M25.531",
  "ICD10_R00.2",
  "ICD10_M25.579",
  "ICD10_T14.8XXA",
  "ICD10_E78.5",
  "ICD10_N20.0",
  "ICD10_R60.9",
  "ICD10_D21.9",
  "ICD10_G89.29",
  "ICD10_Z95.5",
  "ICD10_I20.0",
  "ICD10_R10.2",
  "ICD10_M85.80",
  "ICD10_R09.89",
  "ICD10_R55",
  "ICD10_N18.3",
  "ICD10_E04.1",
  "ICD10_R73.03",
  "ICD10_R51",
  "ICD10_M16.12",
  "ICD10_M79.10",
  "ICD10_R10.12",
  "ICD10_Z12.11",
  "ICD10_J34.89",
  "ICD10_F20.9",
  "ICD10_R42",
  "ICD10_I48.20",
  "ICD10_I50.1",
  "ICD10_N93.9",
  "ICD10_K92.2",
  "ICD10_R10.9",
  "ICD10_J18.1",
  "ICD10_I74.2",
  "ICD10_E11.65",
  "ICD10_J15.9",
  "ICD10_I73.9",
  "ICD10_L97.402",
  "ICD10_N17.9",
  "ICD10_E08.621",
  "ICD10_R73.9",
  "ICD10_E11.621",
  "ICD10_I70.261",
  "ICD10_K57.32",
  "ICD10_K83.1",
  "ICD10_R10.13",
  "ICD10_E80.6",
  "ICD10_K76.89",
  "ICD10_K75.0",
  "ICD10_K81.0",
  "ICD10_S63.259A",
  "ICD10_K43.9",
  "ICD10_E87.5",
  "ICD10_R31.0",
  "ICD10_W19.XXXA",
  "ICD10_R33.9",
  "ICD10_M25.532",
  "ICD10_N30.01",
  "ICD10_N39.0",
  "ICD10_R78.81",
  "ICD10_J96.02",
  "ICD10_M54.9",
  "ICD10_K59.00",
  "ICD10_Z68.43",
  "ICD10_R11.2",
  "ICD10_Z21",
  "ICD10_J91.0",
  "ICD10_J90",
  "ICD10_J45.901",
  "ICD10_T78.40XA",
  "ICD10_H57.9",
  "ICD10_J30.1",
  "ICD10_J45.41",
  "ICD10_J45.20",
  "ICD10_J45.40",
  "ICD10_J45.21",
  "ICD10_E11.51",
  "ICD10_D64.9",
  "ICD10_D86.9",
  "ICD10_N18.5",
  "ICD10_E78.00",
  "ICD10_E87.70",
  "ICD10_E11.36",
  "ICD10_E55.9",
  "ICD10_I50.20",
  "ICD10_D63.1",
  "ICD10_K65.9",
  "ICD10_M54.2",
  "ICD10_R10.84",
  "ICD10_N92.0",
  "ICD10_I63.9",
  "ICD10_R53.1",
  "ICD10_M54.5",
  "ICD10_M54.42",
  "ICD10_M79.602",
  "ICD10_F41.1",
  "ICD10_R06.09",
  "ICD10_N18.9",
  "ICD10_R65.21",
  "ICD10_R41.0",
  "ICD10_F32.9",
  "ICD10_J84.10",
  "ICD10_R06.03",
  "ICD10_R41.82",
  "ICD10_R26.9",
  "ICD10_I63.412",
  "ICD10_J12.9",
  "ICD10_E87.1",
  "ICD10_N13.30",
  "ICD10_R19.7",
  "ICD10_A08.4",
  "ICD10_B34.9",
  "ICD10_I83.893",
  "ICD10_B97.89",
  "ICD10_K21.9",
  "ICD10_E11.40",
  "ICD10_K82.8",
  "ICD10_G47.00",
  "ICD10_E04.2",
  "ICD10_D25.9",
  "ICD10_R10.11",
  "ICD10_K85.90",
  "ICD10_R52",
  "ICD10_E11.10",
  "ICD10_N30.90",
  "ICD10_R56.9",
  "ICD10_K63.1",
  "ICD10_R33.8",
  "ICD10_N40.0",
  "ICD10_F02.81",
  "ICD10_R31.9",
  "ICD10_T83.511D",
  "ICD10_R45.1",
  "ICD10_M79.605",
  "ICD10_K59.09",
  "ICD10_G25.0",
  "ICD10_IMO0001",
  "ICD10_Z91.81",
  "ICD10_G30.9",
  "ICD10_L03.319",
  "ICD10_L02.219",
  "ICD10_R50.9",
  "ICD10_J22",
  "ICD10_B97.29",
  "ICD10_U07.1",
  "ICD10_N18.4",
  "ICD10_N28.9",
  "ICD10_R35.0",
  "ICD10_E03.9",
  "ICD10_N25.81",
  "ICD10_M54.30",
  "ICD10_J11.1",
  "ICD10_R35.8",
  "ICD10_F41.9",
  "ICD10_M79.606",
  "ICD10_R80.9",
  "ICD10_I16.0",
  "ICD10_D50.9",
  "ICD10_H26.9",
  "ICD10_E16.2",
  "ICD10_I11.0",
  "ICD10_E11.3293",
  "ICD10_E11.22",
  "ICD10_M31.6",
  "ICD10_N30.00",
  "ICD10_J40",
  "ICD10_H52.7",
  "ICD10_H81.10",
  "ICD10_R11.0",
  "ICD10_J96.00",
  "ICD10_J81.0",
  "ICD10_J44.1",
  "ICD10_J96.12",
  "ICD10_J44.9",
  "ICD10_I50.32",
  "ICD10_I50.30",
  "ICD10_R53.81",
  "ICD10_R91.8",
  "ICD10_J41.0",
  "ICD10_Z72.0",
  "ICD10_R04.0",
  "ICD10_J45.902",
  "ICD10_R30.0",
  "ICD10_J45.31",
  "ICD10_M54.31",
  "ICD10_F11.20",
  "ICD10_F10.20",
  "ICD10_B37.3",
  "ICD10_G43.909",
  "ICD10_E78.1",
  "ICD10_Z53.21",
  "ICD10_I49.5",
  "ICD10_S09.90XA",
  "ICD10_K42.9",
  "ICD10_K40.91",
  "ICD10_K45.8",
  "ICD10_K52.9",
  "ICD10_I50.9",
  "ICD10_I35.0",
  "ICD10_M25.512",
  "ICD10_R10.10",
  "ICD10_Y92.009",
  "ICD10_R60.0",
  "ICD10_S90.00XA",
  "ICD10_F10.920",
  "ICD10_M25.569",
  "ICD10_M79.671",
  "ICD10_M25.511",
  "ICD10_I50.23",
  "ICD10_G93.41",
  "ICD10_I47.1",
  "ICD10_E86.0",
  "ICD10_K70.31",
  "ICD10_K72.90",
  "ICD10_M25.562",
  "ICD10_S92.213A",
  "ICD10_S02.2XXA",
  "ICD10_H27.119",
  "ICD10_S93.401A",
  "ICD10_N43.3",
  "ICD10_N50.82",
  "ICD10_N43.1",
  "ICD10_N45.4",
  "ICD10_L03.314",
  "ICD10_R06.82",
  "ICD10_I26.99",
  "ICD10_L55.9",
  "ICD10_S22.49XA",
  "ICD10_S30.1XXA",
  "ICD10_M79.672",
  "ICD10_T83.9XXA",
  "ICD10_R10.32",
  "ICD10_G45.9",
  "ICD10_Z86.73",
  "ICD10_M54.16",
  "ICD10_M25.551",
  "ICD10_R79.1",
  "ICD10_R13.10",
  "ICD10_Z51.89",
  "ICD10_K25.3",
  "ICD10_F10.10",
  "ICD10_G93.40",
  "ICD10_E13.10",
  "ICD10_K50.90",
  "ICD10_G89.18",
  "ICD10_E53.8",
  "ICD10_R10.33",
  "ICD10_K92.1",
  "ICD10_K50.919",
  "ICD10_E44.0",
  "ICD10_K50.019",
  "ICD10_E43",
  "ICD10_J02.9",
  "ICD10_R62.7",
  "ICD10_I96",
  "ICD10_M86.9",
  "ICD10_G51.0",
  "ICD10_H10.9",
  "ICD10_I50.89",
  "ICD10_E11.42",
  "ICD10_I34.0",
  "ICD10_I50.33",
  "ICD10_K51.90",
  "ICD10_C18.9",
  "ICD10_Z93.2",
  "ICD10_K51.919",
  "ICD10_K80.50",
  "ICD10_R10.31",
  "ICD10_C18.7",
  "ICD10_C18.5",
  "ICD10_T83.511S",
  "ICD10_K51.911",
  "ICD10_D62",
  "ICD10_K44.9",
  "ICD10_D69.59",
  "ICD10_I82.499",
  "ICD10_J98.11",
  "ICD10_Z98.890",
  "ICD10_K56.609",
  "ICD10_K31.1",
  "ICD10_K56.690",
  "ICD10_K56.699",
  "ICD10_R18.0",
  "ICD10_R58",
  "ICD10_J45.51",
  "ICD10_M32.9",
  "ICD10_K52.89",
  "ICD10_M47.9",
  "ICD10_M75.30",
  "ICD10_N28.1",
  "ICD10_M17.10",
  "ICD10_L30.8",
  "ICD10_M79.7",
  "ICD10_M79.641",
  "ICD10_G40.909",
  "ICD10_G62.9",
  "ICD10_G47.33",
  "ICD10_S93.609A",
  "ICD10_D17.1",
  "ICD10_M70.50",
  "ICD9_428.0",
  "ICD9_428.32",
  "ICD10_D17.9",
  "ICD10_R25.1",
  "ICD10_G40.319",
  "ICD10_M77.02",
  "ICD10_R19.00",
  "ICD10_M17.11",
  "ICD10_M54.6",
  "ICD10_M79.601",
  "ICD10_K21.0",
  "ICD10_E08.9",
  "ICD10_M75.42",
  "ICD10_G40.509",
  "ICD10_M75.32",
  "ICD10_R53.83",
  "ICD10_S80.01XA",
  "ICD10_K25.4",
  "ICD10_J20.9",
  "ICD10_M25.571",
  "ICD10_M00.9",
  "ICD10_R07.1",
  "ICD10_T78.3XXA",
  "ICD10_R07.2",
  "ICD10_M17.12",
  "ICD10_Z87.09",
  "ICD10_F41.0",
  "ICD10_R94.31",
  "ICD10_R18.8",
  "ICD10_B19.20",
  "ICD10_K74.69",
  "ICD10_D69.6",
  "ICD10_K72.00",
  "ICD10_B18.2",
  "ICD10_C22.0",
  "ICD10_Z94.4",
  "ICD10_I51.7",
  "ICD10_L08.9",
  "ICD10_L97.509",
  "ICD10_E11.8",
  "ICD10_I25.83",
  "ICD10_Z79.01",
  "ICD10_R26.2",
  "ICD10_E08.10",
  "ICD10_L03.032",
  "ICD10_M10.9",
  "ICD10_Z45.02",
  "ICD10_B37.0",
  "ICD10_M11.20",
  "ICD10_M62.838",
  "ICD10_B20",
  "ICD10_R39.15",
  "ICD10_S02.609A",
  "ICD10_R31.29",
  "ICD10_N19",
  "ICD10_D68.59",
  "ICD10_Z23",
  "ICD10_I42.9",
  "ICD10_J96.91",
  "ICD10_Z11.1",
  "ICD10_B59",
  "ICD10_N52.9",
  "ICD10_N28.89",
  "ICD10_I15.1",
  "ICD10_R04.2",
  "ICD10_S82.141A",
  "ICD10_T81.31XA",
  "ICD10_H44.23",
  "ICD10_H25.9",
  "ICD10_H18.463",
  "ICD10_Z94.5",
  "ICD10_Z95.810",
  "ICD10_A63.0",
  "ICD10_K60.2",
  "ICD10_K12.2",
  "ICD10_H10.32",
  "ICD10_S01.512A",
  "ICD10_I42.0",
  "ICD10_L02.416",
  "ICD10_J96.11",
  "ICD10_H33.022",
  "ICD10_L02.415",
  "ICD10_D69.3",
  "ICD10_M25.572",
  "ICD10_L03.012",
  "ICD10_S40.019A",
  "ICD10_M79.642",
  "ICD10_M79.18",
  "ICD10_K64.0",
  "ICD10_S40.012A",
  "ICD10_R10.30",
  "ICD10_N32.89",
  "ICD10_A07.4",
  "ICD10_N48.89",
  "ICD10_R63.4",
  "ICD10_G35",
  "ICD10_T84.498A",
  "ICD10_L97.313",
  "ICD10_E87.6",
  "ICD10_M79.604",
  "ICD10_M48.00",
  "ICD10_M54.41",
  "ICD10_M79.609",
  "ICD10_I50.21",
  "ICD10_I87.8",
  "ICD10_R29.6",
  "ICD10_S80.811A",
  "ICD10_K29.80",
  "ICD10_K27.9",
  "ICD10_S01.111A",
  "ICD10_I61.9",
  "ICD10_E11.00",
  "ICD10_D3A.00",
  "ICD10_Z78.9",
  "ICD10_K65.8",
  "ICD10_K80.80",
  "ICD10_R94.5",
  "ICD10_K57.30",
  "ICD10_D72.825",
  "ICD10_R91.1",
  "ICD10_K91.2",
  "ICD10_B49",
  "ICD10_E44.1",
  "ICD10_E03.4",
  "ICD10_T85.9XXA",
  "ICD10_M05.20",
  "ICD10_G47.01",
  "ICD10_R00.0",
  "ICD10_J95.821",
  "ICD10_K80.12",
  "ICD10_T81.89XD",
  "ICD10_M06.9",
  "ICD10_M25.552",
  "ICD10_K62.5",
  "ICD10_M19.90",
  "ICD10_G81.90",
  "ICD10_R32",
  "ICD10_I60.8",
  "ICD10_Z85.3",
  "ICD10_H40.9",
  "ICD10_B37.49",
  "ICD10_J96.22",
  "ICD10_R21",
  "ICD10_I50.31",
  "ICD10_E66.2",
  "ICD10_M15.0",
  "ICD10_G81.14",
  "ICD10_N39.44",
  "ICD10_I63.89",
  "ICD10_G82.20",
  "ICD10_N31.9",
  "ICD10_K59.04",
  "ICD10_N39.46",
  "ICD10_F22",
  "ICD10_F01.51",
  "ICD10_K31.819",
  "ICD10_Z63.6",
  "ICD10_Z74.01",
  "ICD10_N12",
  "ICD10_K59.01",
  "ICD10_I47.9",
  "ICD10_J42",
  "ICD10_M47.12",
  "ICD10_R07.81",
  "ICD10_R29.898",
  "ICD10_R20.0",
  "ICD10_S39.011A",
  "ICD10_I26.94",
  "ICD10_E87.79",
  "ICD10_M86.172",
  "ICD10_T50.902A",
  "ICD10_J98.8",
  "ICD10_K43.2",
  "ICD10_I25.119",
  "ICD10_Z95.2",
  "ICD10_R74.0",
  "ICD10_R09.82",
  "ICD10_I82.402",
  "ICD10_J10.1",
  "ICD10_K46.0",
  "ICD10_IMO0002",
  "ICD10_N61.1",
  "ICD10_L02.91",
  "ICD10_E11.319",
  "ICD10_H53.8",
  "ICD10_R14.0",
  "ICD10_S16.1XXA",
  "ICD10_S10.93XA",
  "ICD10_S00.93XA",
  "ICD10_M25.559",
  "ICD10_D05.92",
  "ICD10_S99.921A",
  "ICD10_D05.12",
  "ICD10_F17.210",
  "ICD10_E08.00",
  "ICD10_N32.0",
  "ICD10_R40.4",
  "ICD10_I27.20",
  "ICD10_G95.9",
  "ICD10_R40.1",
  "ICD10_J02.0",
  "ICD10_L50.0",
  "ICD10_L29.9",
  "ICD10_G56.03",
  "ICD10_S93.409A",
  "ICD10_I83.90",
  "ICD10_F14.11",
  "ICD10_F17.200",
  "ICD10_N92.6",
  "ICD10_Z01.419",
  "ICD10_M79.669",
  "ICD10_F19.20",
  "ICD10_M79.89",
  "ICD10_J45.30",
  "ICD10_S80.12XA",
  "ICD10_M25.541",
  "ICD10_L97.909",
  "ICD10_T45.1X5A",
  "ICD10_G62.0",
  "ICD10_R57.0",
  "ICD10_E87.2",
  "ICD10_T56.894A",
  "ICD10_F31.31",
  "ICD10_F31.32",
  "ICD10_F31.4",
  "ICD10_F33.9",
  "ICD10_A60.04",
  "ICD10_B97.4",
  "ICD10_N23",
  "ICD10_V87.7XXA",
  "ICD10_M25.519",
  "ICD10_K63.5",
  "ICD10_L02.211",
  "ICD10_H60.399",
  "ICD10_Z79.899",
  "ICD10_Z29.8",
  "ICD10_Z94.0",
  "ICD10_R82.90",
  "ICD10_T83.511A",
  "ICD10_N18.2",
  "ICD10_I12.9",
  "ICD10_Q64.2",
  "ICD10_N18.1",
  "ICD10_F81.9",
  "ICD10_B96.20",
  "ICD10_F03.91",
  "ICD10_E87.0",
  "ICD10_I26.09",
  "ICD10_S32.591A",
  "ICD10_S32.591D",
  "ICD10_D49.4",
  "ICD10_H15.019",
  "ICD9_V42.0",
  "ICD10_N05.9",
  "ICD10_K57.92",
  "ICD10_K80.20",
  "ICD10_C64.9",
  "ICD10_Z80.42",
  "ICD10_E83.42",
  "ICD9_278.02",
  "ICD10_E21.3",
  "ICD10_Z94.9",
  "ICD10_M06.00",
  "ICD10_M34.1",
  "ICD10_M54.10",
  "ICD10_M81.0",
  "ICD10_F43.21",
  "ICD10_M46.1",
  "ICD10_M47.817",
  "ICD10_M47.812",
  "ICD10_H66.90",
  "ICD10_M54.12",
  "ICD10_R06.2",
  "ICD10_M47.819",
  "ICD10_K11.20",
  "ICD10_M54.32",
  "ICD10_H60.392",
  "ICD10_G56.00",
  "ICD10_H92.02",
  "ICD10_I51.89",
  "ICD10_S43.006A",
  "ICD10_M19.019",
  "ICD10_S40.012D",
  "ICD10_I24.9",
  "ICD10_Z95.1",
  "ICD10_I20.8",
  "ICD10_E11.3299",
  "ICD10_M79.2",
  "ICD10_R20.2",
  "ICD10_L02.619",
  "ICD10_L97.519",
  "ICD10_R94.39",
  "ICD10_Z93.0",
  "ICD10_R17",
  "ICD10_D57.1",
  "ICD10_I67.5",
  "ICD10_D57.00",
  "ICD9_517.3",
  "ICD9_282.60",
  "ICD9_282.62",
  "ICD10_F43.22",
  "ICD10_Z86.2",
  "ICD10_G40.409",
  "ICD10_G43.009",
  "ICD10_D57.01",
  "ICD10_S42.322A",
  "ICD10_O22.30",
  "ICD10_N05.1",
  "ICD10_E83.51",
  "ICD10_F33.1",
  "ICD10_H90.3",
  "ICD10_M19.049",
  "ICD10_J38.01",
  "ICD10_H91.93",
  "ICD10_J38.00",
  "ICD10_J38.3",
  "ICD10_Z78.0",
  "ICD10_E03.8",
  "ICD10_Z16.12",
  "ICD10_B96.29",
  "ICD10_I21.3",
  "ICD10_I25.110",
  "ICD10_Z46.6",
  "ICD10_E11.69",
  "ICD10_M86.171",
  "ICD10_L02.611",
  "ICD10_L03.031",
  "ICD10_M15.8",
  "ICD10_M79.675",
  "ICD10_I70.245",
  "ICD10_Z96.652",
  "ICD10_I62.00",
  "ICD10_S06.5X9A",
  "ICD10_J01.10",
  "ICD10_N40.1",
  "ICD10_Y82.9",
  "ICD10_D72.821",
  "ICD10_L84",
  "ICD10_M79.673",
  "ICD10_L98.499",
  "ICD10_H91.90",
  "ICD10_L02.31",
  "ICD10_H60.93",
  "ICD10_L20.89",
  "ICD10_T83.010A",
  "ICD10_L51.1",
  "ICD10_H25.11",
  "ICD10_T83.011A",
  "ICD10_T83.9XXS",
  "ICD10_T83.091D",
  "ICD10_I63.312",
  "ICD10_I63.512",
  "ICD10_T79.6XXA",
  "ICD10_R45.851",
  "ICD10_F32.3",
  "ICD10_R59.1",
  "ICD10_S22.41XA",
  "ICD10_K35.80",
  "ICD10_K35.890",
  "ICD10_M50.322",
  "ICD10_J35.8",
  "ICD10_R29.818",
  "ICD10_I60.9",
  "ICD10_R19.5",
  "ICD10_T82.898A",
  "ICD10_I51.9",
  "ICD10_I11.9",
  "ICD10_I25.118",
  "ICD10_R22.0",
  "ICD10_K02.9",
  "ICD10_M51.26",
  "ICD10_I44.7",
  "ICD10_H92.22",
  "ICD10_R07.82",
  "ICD10_I20.9",
  "ICD10_S02.5XXA",
  "ICD10_E11.628",
  "ICD10_G44.89",
  "ICD10_R25.2",
  "ICD10_F03.90",
  "ICD10_Z01.818",
  "ICD10_I82.401",
  "ICD10_S62.009A",
  "ICD10_S62.90XA",
  "ICD10_M25.549",
  "ICD10_D52.9",
  "ICD10_D72.89",
  "ICD10_N64.4",
  "ICD10_F32.1",
  "ICD9_729.5",
  "ICD10_E10.9",
  "ICD10_S69.92XA",
  "ICD10_N39.3",
  "ICD10_R87.613",
  "ICD10_K86.89",
  "ICD10_I89.0",
  "ICD10_N76.0",
  "ICD10_B35.6",
  "ICD10_J69.0",
  "ICD10_C7B.8",
  "ICD10_Z68.41",
  "ICD10_R26.89",
  "ICD10_D3A.8",
  "ICD10_Z98.84",
  "ICD10_W19.XXXD",
  "ICD10_R09.81",
  "ICD10_H57.89",
  "ICD10_N75.1",
  "ICD10_L72.3",
  "ICD10_H10.13",
  "ICD10_S46.912A",
  "ICD10_R69",
  "ICD10_I26.93",
  "ICD10_K74.60",
  "ICD10_R74.8",
  "ICD10_Z00.6",
  "ICD10_Z95.811",
  "ICD9_428.9",
  "ICD9_402.91",
  "ICD10_R20.9",
  "ICD10_T82.897A",
  "ICD9_428.22",
  "ICD9_428.42",
  "ICD10_T82.9XXS",
  "ICD10_T80.219A",
  "ICD10_T84.7XXA",
  "ICD10_T84.7XXD",
  "ICD10_M25.00",
  "ICD10_K66.1",
  "ICD10_H01.003",
  "ICD10_H01.006",
  "ICD10_D64.89",
  "ICD10_A49.8",
  "ICD10_K29.71",
  "ICD10_R65.10",
  "ICD10_Z76.82",
  "ICD10_T82.9XXD",
  "ICD10_T73.2XXD",
  "ICD10_H43.393",
  "ICD10_D63.8",
  "ICD10_I63.40",
  "ICD10_Z94.1",
  "ICD10_E08.3521",
  "ICD10_C61",
  "ICD10_A41.52",
  "ICD10_K65.1",
  "ICD10_C78.00",
  "ICD10_L03.818",
  "ICD10_I63.539",
  "ICD10_N04.1",
  "ICD9_401.9",
  "ICD10_E13.9",
  "ICD10_T86.19",
  "ICD10_I70.1",
  "ICD10_K91.0",
  "ICD10_L89.159",
  "ICD10_T81.49XA",
  "ICD10_R41.3",
  "ICD10_M48.061",
  "ICD10_G82.54",
  "ICD10_E06.3",
  "ICD10_I42.8",
  "ICD10_G47.9",
  "ICD10_K92.0",
  "ICD10_E89.0",
  "ICD10_B30.9",
  "ICD10_Z99.89",
  "ICD10_J31.0",
  "ICD10_K76.0",
  "ICD10_E66.09",
  "ICD10_F17.290",
  "ICD10_G82.22",
  "ICD10_K85.91",
  "ICD10_J96.21",
  "ICD10_S06.0X0A",
  "ICD10_D24.1",
  "ICD10_Q61.3",
  "ICD10_G25.3",
  "ICD10_K86.1",
  "ICD10_Z87.891",
  "ICD10_K06.8",
  "ICD10_I38",
  "ICD10_D50.8",
  "ICD10_E11.59",
  "ICD10_I48.3",
  "ICD10_I62.9",
  "ICD10_G81.94",
  "ICD10_I67.9",
  "ICD10_Z86.69",
  "ICD10_I69.30",
  "ICD10_C50.919",
  "ICD10_C50.912",
  "ICD10_I26.92",
  "ICD10_H20.9",
  "ICD10_I85.10",
  "ICD10_I25.2",
  "ICD10_F05",
  "ICD10_S70.02XA",
  "ICD10_S62.307A",
  "ICD10_R30.9",
  "ICD10_T83.090A",
  "ICD10_G83.9",
  "ICD10_Z43.3",
  "ICD10_I71.9",
  "ICD10_G83.10",
  "ICD10_Z93.3",
  "ICD10_Z86.79",
  "ICD10_M46.20",
  "ICD10_I16.1",
  "ICD10_I82.890",
  "ICD10_B00.9",
  "ICD10_N45.3",
  "ICD10_T16.2XXA",
  "ICD10_M79.644",
  "ICD9_428.43",
  "ICD10_K81.9",
  "ICD10_G73.7",
  "ICD10_Z55.0",
  "ICD10_E03.2",
  "ICD10_E85.0",
  "ICD10_I51.3",
  "ICD10_Z95.820",
  "ICD10_I63.449",
  "ICD10_G47.39",
  "ICD10_I99.8",
  "ICD10_G93.89",
  "ICD10_C71.1",
  "ICD10_B02.9",
  "ICD10_N17.0",
  "ICD10_R00.1",
  "ICD10_R59.9",
  "ICD10_K29.21",
  "ICD10_I71.01",
  "ICD10_I71.00",
  "ICD10_T79.0XXA",
  "ICD10_K40.90",
  "ICD10_I71.03",
  "ICD10_I63.10",
  "ICD10_L97.319",
  "ICD10_I83.013",
  "ICD10_M75.01",
  "ICD10_B37.2",
  "ICD10_M54.40",
  "ICD10_J20.8",
  "ICD10_L89.154",
  "ICD10_I50.43",
  "ICD10_I21.02",
  "ICD10_R11.15",
  "ICD10_K40.30",
  "ICD10_C91.10",
  "ICD10_R79.9",
  "ICD10_M25.861",
  "ICD10_M00.862",
  "ICD10_C55",
  "ICD10_C54.1",
  "ICD10_R22.1",
  "ICD10_J96.90",
  "ICD10_B34.2",
  "ICD10_R40.0",
  "ICD10_F11.10",
  "ICD10_K75.9",
  "ICD10_K85.10",
  "ICD10_K85.20",
  "ICD10_J84.9",
  "ICD10_S40.011A",
  "ICD10_R62.51",
  "ICD10_I95.89",
  "ICD10_J11.00",
  "ICD10_D12.6",
  "ICD10_G91.2",
  "ICD10_Z87.898",
  "ICD10_T46.0X1A",
  "ICD10_D72.829",
  "ICD10_I77.0",
  "ICD10_Z83.71",
  "ICD10_Z09",
  "ICD10_G62.81",
  "ICD10_C54.9",
  "ICD10_C53.0",
  "ICD10_S72.142A",
  "ICD10_S72.142D",
  "ICD10_M00.819",
  "ICD10_I82.452",
  "ICD10_Z93.1",
  "ICD10_L89.150",
  "ICD10_I85.00",
  "ICD10_K76.6",
  "ICD10_R89.9",
  "ICD10_Z59.0",
  "ICD10_Z65.9",
  "ICD10_M86.10",
  "ICD10_L98.429",
  "ICD10_T86.298",
  "ICD10_I47.2",
  "ICD10_E13.11",
  "ICD10_S22.050A",
  "ICD10_K75.81",
  "ICD10_R60.1",
  "ICD10_R11.11",
  "ICD10_H53.9",
  "ICD10_L02.411",
  "ICD10_C90.01",
  "ICD10_C90.00",
  "ICD10_B19.11",
  "ICD10_K83.09",
  "ICD10_K35.30",
  "ICD10_E05.90",
  "ICD10_S30.0XXA",
  "ICD10_J68.3",
  "ICD10_F41.8",
  "ICD10_R03.0",
  "ICD10_J45.32",
  "ICD10_Z63.8",
  "ICD10_N84.0",
  "ICD10_S03.40XA",
  "ICD10_Z17.0",
  "ICD10_C50.811",
  "ICD10_J30.89",
  "ICD10_Z43.1",
  "ICD10_L30.9",
  "ICD10_Y95",
  "ICD10_S72.001A",
  "ICD10_S72.141A",
  "ICD10_T82.9XXA",
  "ICD10_M65.30",
  "ICD10_M20.012",
  "ICD10_S62.101A",
  "ICD10_Z74.09",
  "ICD10_S52.501A",
  "ICD10_S32.029A",
  "ICD10_S32.001S",
  "ICD10_S82.851A",
  "ICD10_F10.929",
  "ICD10_K29.20",
  "ICD10_Z76.5",
  "ICD10_M25.522",
  "ICD10_G44.229",
  "ICD10_K85.00",
  "ICD10_B00.89",
  "ICD10_T78.01XS",
  "ICD10_T78.00XA",
  "ICD10_I61.4",
  "ICD10_S82.899A",
  "ICD10_E83.52",
  "ICD10_T82.838A",
  "ICD10_I65.22",
  "ICD10_L97.511",
  "ICD10_L97.901",
  "ICD10_J44.0",
  "ICD10_G25.5",
  "ICD10_M25.40",
  "ICD10_J30.2",
  "ICD10_E08.22",
  "ICD10_D49.6",
  "ICD10_F32.89",
  "ICD10_B19.10",
  "ICD10_I63.521",
  "ICD10_I61.2",
  "ICD10_B88.8",
  "ICD10_Z59.3",
  "ICD10_F02.80",
  "ICD10_H04.123",
  "ICD10_S62.231A",
  "ICD10_T74.21XA",
  "ICD10_C16.2",
  "ICD10_J47.9",
  "ICD10_D69.1",
  "ICD10_A04.72",
  "ICD10_L91.0",
  "ICD10_J45.42",
  "ICD10_H40.009",
  "ICD10_H35.40",
  "ICD10_L03.039",
  "ICD10_S61.219A",
  "ICD10_Z48.02",
  "ICD10_S98.132A",
  "ICD10_S81.802A",
  "ICD10_S91.302S",
  "ICD9_428.23",
  "ICD10_Z90.10",
  "ICD10_C50.911",
  "ICD10_N83.209",
  "ICD10_Z90.11",
  "ICD10_R87.810",
  "ICD10_R87.610",
  "ICD10_D06.9",
  "ICD10_B97.7",
  "ICD10_N83.201",
  "ICD10_N94.89",
  "ICD10_N85.8",
  "ICD10_I63.411",
  "ICD10_I63.311",
  "ICD10_I63.319",
  "ICD10_E11.29",
  "ICD10_R11.10",
  "ICD10_Z11.3",
  "ICD10_K63.89",
  "ICD10_C20",
  "ICD10_C78.7",
  "ICD10_G89.3",
  "ICD10_C79.51",
  "ICD10_J93.0",
  "ICD10_C79.31",
  "ICD10_C79.49",
  "ICD10_C79.9",
  "ICD10_Z85.6",
  "ICD10_D71",
  "ICD10_J18.0",
  "ICD10_R57.9",
  "ICD10_L97.401",
  "ICD10_W54.0XXA",
  "ICD10_Z82.49",
  "ICD10_I21.9",
  "ICD10_N21.8",
  "ICD10_J00",
  "ICD10_S10.92XA",
  "ICD10_G08",
  "ICD10_S00.82XA",
  "ICD10_S00.02XA",
  "ICD10_J12.89",
  "ICD10_T42.8X1A",
  "ICD10_G82.50",
  "ICD10_M46.46",
  "ICD10_S34.109D",
  "ICD10_K37",
  "ICD10_K64.8",
  "ICD10_I30.9",
  "ICD10_I07.1",
  "ICD10_R53.82",
  "ICD10_D68.32",
  "ICD10_T45.515A",
  "ICD10_I95.0",
  "ICD10_D89.9",
  "ICD10_T86.21",
  "ICD10_I67.1",
  "ICD10_C67.2",
  "ICD10_Z45.09",
  "ICD10_I72.1",
  "ICD10_S40.022A",
  "ICD10_I77.71",
  "ICD10_K83.8",
  "ICD10_I72.9",
  "ICD10_T82.7XXD",
  "ICD10_T82.7XXS",
  "ICD10_I66.02",
  "ICD10_G81.91",
  "ICD10_I71.4",
  "ICD10_N81.4",
  "ICD10_S23.9XXA",
  "ICD10_V89.2XXA",
  "ICD10_S13.9XXA",
  "ICD10_M62.830",
  "ICD10_J93.83",
  "ICD10_T82.590A",
  "ICD10_Z98.62",
  "ICD10_I70.0",
  "ICD10_M43.02",
  "ICD10_J01.00",
  "ICD10_M35.00",
  "ICD10_M17.0",
  "ICD10_R13.12",
  "ICD10_H16.123",
  "ICD10_H52.223",
  "ICD10_H52.213",
  "ICD10_G31.83",
  "ICD10_H40.039",
  "ICD10_H53.40",
  "ICD10_H25.13",
  "ICD10_H40.053",
  "ICD10_H50.10",
  "ICD10_M50.30",
  "ICD10_R41.89",
  "ICD10_H25.813",
  "ICD10_M34.9",
  "ICD10_H66.002",
  "ICD10_G50.1",
  "ICD10_W46.1XXA",
  "ICD10_Z20.6",
  "ICD10_G93.3",
  "ICD10_M79.5",
  "ICD10_J38.5",
  "ICD10_J38.4",
  "ICD10_J95.2",
  "ICD10_Z92.241",
  "ICD10_S61.001A",
  "ICD10_G06.0",
  "ICD10_I49.9",
  "ICD10_M20.022",
  "ICD10_F14.10",
  "ICD10_I82.511",
  "ICD10_L97.919",
  "ICD10_L97.929",
  "ICD10_I83.009",
  "ICD10_L97.921",
  "ICD10_K57.20",
  "ICD10_K57.90",
  "ICD10_B19.9",
  "ICD10_T83.9XXD",
  "ICD10_S62.304A",
  "ICD10_S61.213A",
  "ICD10_G40.901",
  "ICD10_J96.92",
  "ICD10_J21.0",
  "ICD10_I50.40",
  "ICD10_I83.019",
  "ICD10_E08.01",
  "ICD10_K13.79",
  "ICD10_N13.5",
  "ICD10_W46.0XXA",
  "ICD10_S01.00XA",
  "ICD10_M95.2",
  "ICD10_G91.9",
  "ICD10_Z48.89",
  "ICD10_K50.812",
  "ICD10_K50.811",
  "ICD10_R19.07",
  "ICD10_M48.062",
  "ICD10_E61.1",
  "ICD10_I82.90",
  "ICD10_E02",
  "ICD10_E13.65",
  "ICD10_E08.42",
  "ICD10_E13.29",
  "ICD10_J43.2",
  "ICD10_I15.9",
  "ICD10_E87.8",
  "ICD10_F11.23",
  "ICD10_I61.5",
  "ICD10_Q27.30",
  "ICD10_K50.014",
  "ICD10_K50.814",
  "ICD10_Z45.2",
  "ICD10_M62.82",
  "ICD10_J98.2",
  "ICD10_S22.089S",
  "ICD10_K82.A2",
  "ICD10_I70.209",
  "ICD10_I70.229",
  "ICD10_T42.6X5A",
  "ICD10_S91.302D",
  "ICD10_I70.239",
  "ICD10_I25.5"
)
toremove = append(toremove, notneededicd)

print("structure of 850 variables")
str(lbpcov[, toremove])

print("names to remove")
toremove

### flipped cases and controls. Cases is now 0, Controls is now 1
lbpcov = lbpcov %>% mutate(max_aki_outcome = ifelse(max_aki_stage %in% c("1", "2", "3"), "case", "ctrl"))
#lbpcov$max_aki_outcome <- as.numeric(as.character(lbpcov$max_aki_outcome))
lbpcov$max_aki_outcome<-factor(lbpcov$max_aki_outcome, levels=c("ctrl","case"))

print("find levels:")
head(as.factor(lbpcov$max_aki_outcome))
#head(as.factor(lbpcov$COVID19_Days_Since_Most_Recent_Antibody))

#quit(status=1)

##ensure RIN is numeric
lbpcov$RIN <- as.numeric(as.character(lbpcov$RIN))

#head(lbpcov)
print("removing columns where factor levels are 1 or 0")
lbpcov = lbpcov[, !colnames(lbpcov) %in% toremove] 


#print("remove variables except for Mean, Median and StDEV since theyre all correlated.")
#labtablesnotneeded = lbpcov %>% dplyr::select(starts_with(c("RNA_Uncorrected","ICD","NURSING","ECG","Pulse","Immunization","Post_COVID19","Viral_Load")),ends_with(c("StDev","Max","Num_Observations","IQR","per_mL","Date")))
##print(colnames(labtablesnotneeded))

###now remove everything else in the above df from the main DF. 
#print("removing them from main df")
#lbpcov = lbpcov[, !colnames(lbpcov) %in% colnames(labtablesnotneeded)]


##we should now have only ~300 columns or so
print("check to see how many columns have na values")

nrow(lbpcov)
dim(lbpcov)
print("nacount")
na_count_df <-sapply(lbpcov, function(y) sum(length(which(is.na(y)))))
na_count_df = data.frame(na_count_df)
dim(na_count_df)
colnames(na_count_df) = c("na_count")
data.frame(na_count_df)
print("filtering to remove 90% of NA")

#Change it back to 30!!!
#na_count_filter = na_count_df %>% dplyr::filter(na_count <= 30)
na_count_filter = na_count_df %>% dplyr::filter(na_count <= 100)
data.frame(na_count_filter)

final_cols_toselect = rownames(na_count_filter)
final_cols_toselect

####FINAL METADATA
lbpcov = lbpcov[,final_cols_toselect]


####RUNNIG ADDITIONAL STEP 09/03/2022 PJ
print("remove variables except for Mean, Median and StDEV since theyre all correlated.")
labtablesnotneeded = lbpcov %>% dplyr::select(starts_with(c("RNA_Uncorrected","ICD","NURSING","ECG","Pulse","Immunization","Post_COVID19","Viral_Load")),ends_with(c("StDev","Max","Num_Observations","IQR","per_mL","Date")))
#print(colnames(labtablesnotneeded))

##now remove everything else in the above df from the main DF. 
print("removing them from main df")
lbpcov = lbpcov[, !colnames(lbpcov) %in% colnames(labtablesnotneeded)]


#### ADDITIONAL STEP 


print("look for var = 0 in colnames with numeric class")

rnaseq_samplenames = lbpcov$RNA_Corrected_Sample
#rownames(lbpcov) = rnaseq_samplenames
#make the rownames of lbpcov the SAMPLE_ISMMS ids in order for the formula below to work
rownames(lbpcov) <- lbpcov$RNA_Corrected_Sample

##transpose
lbpcov_tr = data.frame(t(lbpcov))
colnames(lbpcov_tr) = rnaseq_samplenames
#head(lbpcov)

dim(lbpcov)



### read rnaseq expression data
covariate_col = colnames(lbpcov_tr)
fc_orig <- read.csv(file = "/sc/arion/projects/mscic1/results/Noam/MainCovid/expression/MainCovid_all_counts_matrix_2020-08-18_mislabelCorrected_renamed.txt", sep='\t', header=TRUE, stringsAsFactors = FALSE)
#colnames(fc_orig)
#rownames(fc_orig)
rownames(fc_orig) = fc_orig$X
fc = data.frame(fc_orig)
dim(fc)

##remove genes in the PAR region
fc_discard = fc[grepl("_PAR_", fc$X),]
fc = fc[!grepl("_PAR_", fc$X),]
dim(fc_discard)

fc = fc %>% separate(X, c("gene", NA))

head(fc_discard[,1:10])

rownames(fc) = fc$gene
dim(fc)

#rownames(fc)

rnaseq_ct_col = colnames(fc)

#fc = fc_orig[,lbpcov$RNA_Corrected_Sample]
#fc$Geneid <- NULL
#lbpcov <- lbpcov[match(colnames(fc), lbpcov$SAMPLE_ISMMS),]
newidx = match(covariate_col, rnaseq_ct_col)
newidx = newidx[!is.na(newidx)]
newidx
print("now getting df to match order")
#head(fc)
dim(fc)
#fc_ordered = fc %>% dplyr::select(newidx)
fc_ordered = fc[, newidx]
print("checking new ordered df")
head(fc_ordered)
fc = fc_ordered

print("getting index of the filtered not NA samples to filter dataframe rows")
newrowidx = match(colnames(fc), rownames(lbpcov))
newrowidx
newrowidx = newrowidx[!is.na(newrowidx)]
print("now getting df to match order of rows")
lbpcov_ordered = lbpcov[newrowidx,]

lbpcov = lbpcov_ordered
#lbpcov[1:20,]

#make the rownames (again just in case) of lbpcov the SAMPLE_ISMMS ids in order for the formula below to work
rownames(lbpcov) <- lbpcov$RNA_Corrected_Sample
print("testing for congruency in list")

rownames(lbpcov)
colnames(fc)
print("checking if rownames lbpcov are same as colnames fc - they should be!!")
identical(rownames(lbpcov), colnames(fc))

print("drop unused levels in metadata")
lbpcov = droplevels(lbpcov)

####make all factors
#lbpcov = lbpcov %>% dplyr::mutate_all(as.factor) %>% purrr::map(levels)
colnames(lbpcov)
dim(lbpcov)


##########################################################################################################BEGIN TO TEST MODEL


metadata <- lbpcov
countMatrix <- fc

saveRDS(metadata, file=paste(results_dir,"metadata.rds", sep=""))
saveRDS(countMatrix, file=paste(results_dir,"geneexpr.rds", sep=""))




#####testing correlations####


form2=as.formula(paste("~",paste(colnames(lbpcov)[colnames(lbpcov) %in% c('ventilation','severity','O2_Therapy','WHO_Ordinal')],collapse="+")))
form2

C2 = canCorPairs(form2,lbpcov[,colnames(lbpcov) %in% c('ventilation','severity','O2_Therapy','WHO_Ordinal')])
C2

#paste(plotpath, "lbp_allBatches_QC_pca-plot.pdf",sep="")
print("plotting CanCorPairs ventilation WHO")
pdf(file = paste(results_dir, "Ventilation_WHO_QC_CanCorsHeatmap.pdf", sep=""), height=50,width=50)
plotCorrMatrix(C2,margins = c(25, 25))
#dev.off()

ggplot(lbpcov, aes(x = lbpcov$WHO_Ordinal, y = lbpcov$ventilation, color = lbpcov$max_aki_outcome)) +
  geom_point(size = 20) +  scale_fill_discrete()

dev.off()


###########################################################################################
#canCorPairs Table -- correlations between all covariates
###########################################################################################
###correlations geneartd by cancorpairs - anonical correlation
print("running canCorrPairs now!")
library(variancePartition)
form=as.formula(paste("~",paste(colnames(lbpcov)[colnames(lbpcov) !="RNA_Corrected_Sample"],collapse="+")))
form

C = canCorPairs(form,lbpcov[,colnames(lbpcov)!= "RNA_Corrected_Sample"])

#paste(plotpath, "lbp_allBatches_QC_pca-plot.pdf",sep="")
print("plotting CanCorPairs")
pdf(file = paste(results_dir, "lbp_allBatches_QC_CanCorsHeatmap.pdf", sep=""), height=50,width=50)
plotCorrMatrix(C,margins = c(25, 25))
dev.off()



###COMMENTING OUT FOR NOW!!!!! - PJ 09142021
#C = canCorPairs(form,lbpcov[,colnames(lbpcov)!= "RNA_Corrected_Sample"])

#paste(plotpath, "lbp_allBatches_QC_pca-plot.pdf",sep="")
#print("plotting CanCorPairs")
#pdf(file = paste(results_dir, "lbp_allBatches_QC_CanCorsHeatmap.pdf", sep=""), height=50,width=50)
#plotCorrMatrix(C,margins = c(25, 25))
#dev.off()


table(lbpcov$min_lymphocyte)
table(lbpcov$LAB_EOSINOPHIL_NUMBER_Mean)

### STOP HERE UNTIL YOU FIND ALL YOUR COVARIATES. 
##########################################################################################################
#Everything -- this is all covariates I selected being added to the model, but you will do this iteratively starting with no covariates which you did in the code above and now you will add your first covariate -- sex, and then you'll add your first and second -- sex and rin, and then you'll add your first, second, and third -- sex, rin, X, and so forth. 
###########################################################################################################
#### ADD YOUR DESIRED COVARIATES AND THEN RUN. 


#form <- ~ RNA_Library_Prep_Plate + Neutrophils + RNASEQ_FIRST_OF_PAIR.PF_INDEL_RATE + Plasma_cells + MECHANICAL_VENTILATION + Sex + Macrophages_M0 + Age + T_cells_CD4_memory_activated + ckd + RNASEQ_UNMAPPED_READS + max_aki_outcome
form <- ~ RNA_Library_Prep_Plate + Neutrophils + RNASEQ_FIRST_OF_PAIR.PF_INDEL_RATE + Plasma_cells + max_aki_outcome + Severity + Sex + RNASEQ_WIDTH_OF_95_PERCENT + Age + ckd + RNASEQ_GC_DROPOUT + RNASEQ_MEAN_INSERT_SIZE + Macrophages_M0 + max_aki_outcome 

#design <- model.matrix(~ RNA_Library_Prep_Plate + Neutrophils + RNASEQ_FIRST_OF_PAIR.PF_INDEL_RATE + Plasma_cells + MECHANICAL_VENTILATION + Sex + Macrophages_M0 + Age + T_cells_CD4_memory_activated + ckd + RNASEQ_UNMAPPED_READS + max_aki_outcome, metadata)
design <- model.matrix(~ RNA_Library_Prep_Plate + Neutrophils + RNASEQ_FIRST_OF_PAIR.PF_INDEL_RATE + Plasma_cells + Severity + Sex + RNASEQ_WIDTH_OF_95_PERCENT + Age + ckd + RNASEQ_GC_DROPOUT + RNASEQ_MEAN_INSERT_SIZE + Macrophages_M0 + max_aki_outcome, metadata)

head(design)

rownames(design)
print("print form")
form

## to compute residuals for PCA correlation!

#resdesign <- model.matrix(~ RNA_Library_Prep_Plate + Neutrophils + RNASEQ_FIRST_OF_PAIR.PF_INDEL_RATE + Plasma_cells + MECHANICAL_VENTILATION + Sex + Macrophages_M0 + Age + T_cells_CD4_memory_activated + ckd + RNASEQ_UNMAPPED_READS , metadata)
resdesign <- model.matrix(~ RNA_Library_Prep_Plate + Neutrophils + RNASEQ_FIRST_OF_PAIR.PF_INDEL_RATE + Plasma_cells + Severity + Sex + RNASEQ_WIDTH_OF_95_PERCENT + Age + ckd + RNASEQ_GC_DROPOUT + RNASEQ_MEAN_INSERT_SIZE + Macrophages_M0 , metadata)


####write form to file
cat(toString(form),file=paste(results_dir,"README.txt",sep=""),append=TRUE, sep="\n")

###RESIDUALS - ALL OTHER COVARIATES WITHOUT COVARIATE OF INTERES
#### no resform with voom no Dream : just res - res = residuals(fitmm, vobjDream)
#resform <- ~(1|mymet_sex) + mymet_rin + neuronal + RNASeqMetrics_MEDIAN_3PRIME_BIAS + RNASeqMetrics_PCT_MRNA_BASES + (1|IID_ISMMS) + (1|mymet_depletionbatch)
####resform <- ~(RNASEQ_Run_Sample + RNA_Extraction_Batch + Age)

isexpr <- rowSums(cpm(countMatrix)>=1) >= 0.1*ncol(countMatrix)
# Standard usage of limma/voom
geneExpr = DGEList( countMatrix[isexpr,] )
geneExpr = calcNormFactors( geneExpr )

dim(geneExpr)

#library(BiocParallel)
#Sys.setenv(OMP_NUM_THREADS = 26)

#per Noam, there should be no big difference between these two vobjects because the formula is not used for the vobject, the vobject normalizes the distribution
#vobjDream = voomWithDreamWeights( geneExpr, form, metadata, BPPARAM = MulticoreParam(5))
print("running vOBJDream with VOOM")

print("geneExpr")
print(geneExpr)
dim(geneExpr)

print("design")
print(design)
dim(design)

vobjDream = voom( geneExpr, design, plot=TRUE)
print("head vobj")
head(vobjDream)

print("print coeff of fit")
coef(vobjDream)

print("print vobjdream")
print(vobjDream)

##### running fitmm for object with outcome
fitmm = lmFit(vobjDream, design)
print("get colnames of fitmm - not fitDupCor")
names(fitmm)

print("ebayes")
# Fit Empirical Bayes for moderated t-statistics
fitDupCor <- eBayes( fitmm )

print("names fitDupCor")
names(fitDupCor)

pdf("ebayes.pdf")
plotSA(fitDupCor)
dev.off()



##### running resfitmm for object without outcome

fitmm_res = lmFit(vobjDream, resdesign)
print("get colnames of fitmm - not fitDupCor")
names(fitmm_res)

print("ebayes")
# Fit Empirical Bayes for moderated t-statistics
fitDupCor_res <- eBayes( fitmm_res )

print("names fitDupCor")
names(fitDupCor_res)

pdf("ebayes_res.pdf")
plotSA(fitDupCor_res)
dev.off()


print("fitting residual")
res = residuals(fitDupCor_res, vobjDream)
saveRDS(res, file=paste(results_dir,"res_obj.rds", sep=""))
#res[1:100,]

### Run pca
library(PCAtools)
#pdf(paste(results_dir,"res_PCA.pdf", sep=""))
#print("running res pca")
#pca <- PCAtools::pca(res, metadata = metadata)
#print("done with pca")
#colkey = c('0'='gold', '1'='blue')
#print("creating pairs plot")
#pairsplot(pca,colby = "max_aki_outcome" , legendPosition = "none",colkey=colkey,
#          components = c("PC1","PC2","PC3","PC4"))
#print("done with pairs plot")
#dev.off()


lmgroup_DE <- fitDupCor #PJ NOT ignoring the eBayes step
saveRDS(lmgroup_DE, file=paste(results_dir,"voom_obj.rds", sep=""))

print("printing coefcol")
coefcol <- 'max_aki_outcomecase' #the name of your case/control status column
print(coefcol)


print(lmgroup_DE)

print("Group DE Tab")
Group_DE_tab <- topTable(lmgroup_DE, coef=coefcol, number=nrow(vobjDream))
print(Group_DE_tab)

print("getting de")
de <- data.table( gene = rownames(Group_DE_tab), Group_DE_tab)
de

print("getting exprs data and TopTable")
exprs_data = vobjDream$E
print(exprs_data)

print("toptable values****")
topTable(lmgroup_DE, coef=coefcol, number=1)

Group_de_tab_n1 <- rownames(topTable(lmgroup_DE, coef=coefcol, number=1))
typeMean <- tapply(exprs_data[Group_de_tab_n1,], metadata$max_aki_outcome, mean)
typeMean

#save(de, file = paste(results_dir, "RNASEQ_DEG.rData", sep=""))

de <- de[order(logFC)]
de[adj.P.Val<0.05, DEG:="DEG"]
de[adj.P.Val>0.05, DEG:="NOTDEG"]

### This LFC is flipped because 
#[1] 0 1 0 0 0 1
#Levels: 0 1
#we see the reference level is 0 so positive values means it is down-regulated in ER-. 

de[logFC<0, LFC:="NEGLFC"]
de[logFC>0, LFC:="POSLFC"]

print("print the DEGs -  which de$DEG == DEG")
de[which(de$DEG == "DEG"),]
#[1] 17128 vs 17211 without deletion batch added as cov 

print("calculating mean expression for each gene in case control")
#liv <- metadata[mymet_postmortem==0]$RNA_Corrected_Sample
#pmt <- metadata[mymet_postmortem==1]$RNA_Corrected_Sample
print("cases:")
case <- metadata[metadata$max_aki_outcome=="case",]$RNA_Corrected_Sample ### this means aki is 1, 2 or 3
head(case)

print("controls:")
ctrl <- metadata[metadata$max_aki_outcome=="ctrl",]$RNA_Corrected_Sample #### this means aki is 0
head(ctrl)

print("tmp")
tmp <- data.table(gene=names(rowMeans(res)), computedAvgExp=rowMeans(res))
head(tmp)


print("saving res to rdata")
save(res, file = paste(results_dir, "residuals.rData", sep=""))

print("de")
de <- merge(de, tmp)

print("control_mean")
ctrl_mn <- data.table(gene=names(rowMeans(res[,ctrl])), mean_ctrl=rowMeans(res[,ctrl]))
head(ctrl_mn)
str(ctrl_mn)

print("case_mean")
case_mn <- data.table(gene=names(rowMeans(res[,case])), mean_case=rowMeans(res[,case]))
head(case_mn)
str(case_mn)

print("de again")
de <- merge(merge(de, ctrl_mn), case_mn)
head(de)

print("for each gene, the difference between its mean expression in the max_aki_outcome case or ctrl")
de[,DIFF:=mean_case-mean_ctrl] #for each gene, the difference between its mean expression in the living and postmortem sample

#de$DIFF <- ( as.numeric(as.character(case_mn)) - as.numeric(as.character(ctrl_mn)) )
head(de)

save(de, file = paste(results_dir, "RNASEQ_DEG.rData", sep=""))


print("commented out - plot diff gene expression")
pdf(paste(results_dir, "lbp_allBatches_QC_AllCovs_DIFFPlot.pdf", sep=""))
    # ggplot(de, aes(mean_liv, mean_pm)) + geom_point(size=3, pch=21) + facet_wrap(~DEG+LFC) 
  ggplot(de, aes(DIFF, colour=DEG, fill=DEG)) + geom_freqpoly(binwidth=0.1) + facet_wrap(~DEG, scales="free")
  ggplot(de, aes(logFC, colour=DEG, fill=DEG)) + geom_freqpoly(binwidth=0.1) + facet_wrap(~DEG, scales="free")
  ggplot(de, aes(adj.P.Val, colour=DEG, fill=DEG)) + geom_freqpoly(binwidth=0.1) + facet_wrap(~DEG, scales="free")
dev.off()



vp = fitExtractVarPartModel( vobjDream, form, metadata, BPPARAM = MulticoreParam(5))

pdf(file = paste(results_dir,"lbp_allBatches_vobjDream_variancePartitionPlot.pdf", sep=""), width=25, height=25)
plotVarPart( sortCols(vp))
dev.off()


###############################################
print("starting PCA analysis")
DATA <- metadata


###convert all factors back to to numeric
#lbpcov = mutate_if(lbpcov, is.factor, ~ as.numeric(levels(.x))[.x])


plotpath = results_dir
vobj <- vobjDream

level3=pnorm(3,mean=0,sd=1,lower.tail=T) - pnorm(3,lower.tail=F)
level2=pnorm(2,mean=0,sd=1,lower.tail=T) - pnorm(2,lower.tail=F)
level1=pnorm(1,mean=0,sd=1,lower.tail=T) - pnorm(1,lower.tail=F)
type="norm"


is.empty.vector=function(x) return(length(x)==0)

#clonename<-rownames(SampleByVariable) #if you uncomment this, you will get the sample names next to the dots on the plot
clonename <- NA

###########################
######When no covariates!!
#############################
#print("print $E")
#print(vobj$E)

#print("printing covariates")
#print(cov(vobj$E))
##### THIS IS STEP 1 - with no covariates - basic model. 
##NO COVARIATES. 
SampleByVariable_0=t(cov(vobj$E)) #this is when you have no covariates what you use for calculating the PCs; starting at the next iteration when you add covariates to the model, you will use the residuals for this calculation and for the PCA calculation 
print("sample By Variable")
head(SampleByVariable_0)

#clonename<-rownames(SampleByVariable) #if you uncomment this, you will get the sample names next to the dots on the plot
clonename <- NA

pca_0<- prcomp(SampleByVariable_0, scale=T)
summ_0=summary(pca_0)

#print("summary PCA")
#head(summ_0)

#print("pca 1:5")
#as.data.frame(pca_0$x[,1:5])
print("printing lbpcov")
lbpcov[1:15,]

#The covariates-PCs correlation table
print("printing rescor")
resCor_0=canCorAllAgainstAll_Original(lbpcov,as.data.frame(pca_0$x[,1:5]),minimum_intersect=100)
#head(resCor_0)
print("ordered rescor")
ordered_resCor_0=do.call(order,as.data.frame(-resCor_0))
print("rewriting resCor")
resCor_0=resCor_0[ordered_resCor_0,]

print("resCor is")
#head(resCor_0)

print("data frame")
#as.data.frame(resCor_0, keep.rownames=TRUE)

print("getting rescor and corr ordered by pc1")
all_0 <- as.data.frame(resCor_0, keep.rownames=TRUE)
print("getting C2_0")
C2_0<-as.data.frame(C[rownames(C), "max_aki_outcome"])
colnames(C2_0)[1] <- "aki_corr"
C2_0$covs <- rownames(C2_0)
all_0 <- merge(C2_0, all_0, by.x="covs", by.y=0)
all_0 <- all_0[order(-all_0$PC1),]
head(all_0,50)
write.csv(all_0, file = paste(results_dir,"NOCOVARIATES_CONTRIBUTION2PCA.csv", sep=""))


############################
#STEP2: Once you add covariates!!
##############################
print("prcomp")
covariance=cov(res) #here when you have covariates added to the model, you use the residuals not the vobject to calculate the PCs
SampleByVariable=t(covariance)
pca <- prcomp(SampleByVariable, scale=T)
summ=summary(pca)
print(summ)

print("running resCor")
resCor=canCorAllAgainstAll_Original(lbpcov,as.data.frame(pca$x[,1:5]),minimum_intersect=100)
ordered_resCor=do.call(order,as.data.frame(-resCor))
resCor=resCor[ordered_resCor,]


all <- as.data.frame(resCor, keep.rownames=TRUE)
C2<-as.data.frame(C[rownames(C), "max_aki_outcome"])
colnames(C2)[1] <- "aki_corr"
C2$covs <- rownames(C2)
all <- merge(C2, all, by.x="covs", by.y=0)
all <- all[order(-all$PC1),]
head(all,70)

write.csv(all, file = paste(results_dir,"WITHCOVARIATES_CONTRIBUTION2PCA.csv",sep=""))

all.all <- all
print(all.all)

pdf(paste(plotpath, "lbp_allBatches_QC_pca-plot.pdf",sep=""))
count=1
  total=ncol(DATA)
  for(col in colnames(DATA)){
    print(paste0("on column",count,"/",total,col,"\n"))
    cat("on column",count,"/",total,col,"\n")
  # #=======pca-1 vs pca-2=======
    print("PCA 1 vs PCA 2")
    if(is.numeric(DATA[[col]])==F){
     	print("if DATA[[col]] is not numeric..")
  rn=rownames(pca$x) %in% DATA[!is.na(DATA[[col]]),]$RNA_Corrected_Sample
  a <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,1], y= pca$x[rn,2], color=factor(DATA[[col]]), label=clonename))+
    geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC1-PC2") + xlab(paste("PC1: ",round(summ$importance[2,1]*100,digits=2),"%",sep="")) +
    ylab(paste("PC2: ",round(summ$importance[2,2]*100,digits=2),"%",sep="")) + ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
  if(length(unique(DATA[[col]]))<7){
    print("if data[[col]] has unique length of 7 or less..")
    a = a + scale_colour_discrete(guide ="legend",name=substring(col, 1, 8)) +
    theme(legend.key.size=unit(0.5,"cm")) + guides(color = guide_legend(override.aes = list(size=2))) +
    theme(legend.key.width=unit(0.2,"cm")) +
  labs(fill=c(substr(col,1,12)))
  print("done with checking if col has unique length of less than 7")
  }else{
    print("data[[col]] is numeric")
    a = a + scale_color_discrete(guide =FALSE)+ labs(color=c(substring(col,1,8)))
  }
  build <- ggplot_build(a)$data
  points <- build[[1]]
  ell <- build[[3]]

  print("check which points are inside ellipse")
  # Find which points are inside the ellipse, and add this to the data
  dat <- data.frame(points[1:2], 
                    in.ell = as.logical(point.in.polygon(points$x, points$y, ell$x, ell$y)))
  outliers_3SD_PC1_PC2=points$label[which(dat$in.ell==F)]
  if(is.empty.vector(outliers_3SD_PC1_PC2)==F){
    
  }
  print("done with pca 1 vs pca2")
  # show(a)
  # #=======pca-2 vs pca-3=======
  print("PCA2 vs PCA3")
  b <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,2], y= pca$x[rn,3], color=factor(DATA[[col]]), label=clonename))+
    geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+labs(title="PC2-PC3")+ xlab(paste("PC2: ",round(summ$importance[2,2]*100,digits=2),"%",sep="")) +
    ylab(paste("PC3: ",round(summ$importance[2,3]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey") 
  if(length(unique(DATA[[col]]))<7){
    b = b + scale_colour_discrete(guide ="legend",name=substring(col, 1, 8)) +
    theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(legend.key.width=unit(0.2,"cm")) +
  labs(fill=c(substr(col,1,12)))
  }else{
    b = b + scale_color_discrete(guide =FALSE)+ labs(col=c(substring(col,1,8)))
  }

  print("pca3 vs pca 4")
  # #=======pca-3 vs pca-4=======
  c <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,3], y= pca$x[rn,4], color=factor(DATA[[col]]), label=clonename))+
    geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC3-PC4")+ xlab(paste("PC3: ",round(summ$importance[2,3]*100,digits=2),"%",sep="")) +
    ylab(paste("PC4: ",round(summ$importance[2,4]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey") 
  if(length(unique(DATA[[col]]))<7){
    c = c + scale_colour_discrete(guide ="legend",name=substring(col, 1, 8)) +
    theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(legend.key.width=unit(0.2,"cm")) +
  labs(fill=c(substr(col,1,12)))
  }else{
    c = c + scale_color_discrete(guide =FALSE)+ labs(col=c(substring(col,1,8)))
  }
  # #=======pca-4 vs pca-5=======
  print("pca4 vs PCA5")
  d <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,4], y= pca$x[rn,5], color=factor(DATA[[col]]), label=clonename))+
    geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC4-PC5")+ xlab(paste("PC4: ",round(summ$importance[2,4]*100,digits=2),"%",sep="")) +
    ylab(paste("PC5: ",round(summ$importance[2,5]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey") 
  if(length(unique(DATA[[col]]))<7){
    d = d + scale_colour_discrete(guide ="legend",name=substring(col, 1, 8)) +
    theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(legend.key.width=unit(0.2,"cm")) +
  labs(fill=c(substr(col,1,12)))
  }else{
    d = d + scale_color_discrete(guide =FALSE)+ labs(col=c(substring(col,1,8)))
  }
}else if(is.numeric(DATA[[col]])){
  rn=rownames(pca$x) %in% DATA[!is.na(DATA[[col]]),]$RNA_Corrected_Sample
  a <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,1], y= pca$x[rn,2], color=DATA[[col]][!is.na(DATA[[col]])]))+geom_point(size=3) +
    labs(title="PC1-PC2") + xlab(paste("PC1: ",round(summ$importance[2,1]*100,digits=2),"%",sep="")) +
    ylab(paste("PC2: ",round(summ$importance[2,2]*100,digits=2),"%",sep="")) + ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,1],y=pca$x[rn,2]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
    midpoint=mean(c(min(DATA[[col]][!is.na(DATA[[col]])]),max(DATA[[col]][!is.na(DATA[[col]])])))
    max=max(DATA[[col]][!is.na(DATA[[col]])])
    min=min(DATA[[col]][!is.na(DATA[[col]])])
  a = a + scale_color_gradientn(name =eval(col),colors=c("cyan","orangered1"),labels=c(max,min),breaks=c(max,min)) + labs(col=substring(col,1,8))
  #theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
  #theme(legend.key.width=unit(0.2,"cm"))
  # #=======pca-2 vs pca-3=======
  b <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,2], y= pca$x[rn,3], color=DATA[[col]][!is.na(DATA[[col]])], label=clonename))+geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC2-PC3")+ xlab(paste("PC2: ",round(summ$importance[2,2]*100,digits=2),"%",sep="")) +
    ylab(paste("PC3: ",round(summ$importance[2,3]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,2],y=pca$x[rn,3]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
    midpoint=mean(c(min(DATA[[col]][!is.na(DATA[[col]])]),max(DATA[[col]][!is.na(DATA[[col]])])))
    max=max(DATA[[col]][!is.na(DATA[[col]])])
    min=min(DATA[[col]][!is.na(DATA[[col]])])
  b = b + scale_color_gradientn(name =eval(col),colors=c("cyan","orangered1"),labels=c(max,min),breaks=c(max,min)) + labs(col=substring(col,1,8))
  #theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
  #theme(legend.key.width=unit(0.2,"cm"))
  # #=======pca-3 vs pca-4=======
  c <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,3], y= pca$x[rn,4], color=DATA[[col]][!is.na(DATA[[col]])], label=clonename))+geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC3-PC4")+ xlab(paste("PC3: ",round(summ$importance[2,3]*100,digits=2),"%",sep="")) +
    ylab(paste("PC4: ",round(summ$importance[2,4]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,3],y=pca$x[rn,4]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
    midpoint=mean(c(min(DATA[[col]][!is.na(DATA[[col]])]),max(DATA[[col]][!is.na(DATA[[col]])])))
    max=max(DATA[[col]][!is.na(DATA[[col]])])
    min=min(DATA[[col]][!is.na(DATA[[col]])])
  c = c + scale_color_gradientn(name =eval(col),colors=c("cyan","orangered1"),labels=c(max,min),breaks=c(max,min)) + labs(col=substring(col,1,8))
  #theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
  #theme(legend.key.width=unit(0.2,"cm"))
  # #=======pca-4 vs pca-5=======
  d <- ggplot(data.frame(pca$x[rn,]), aes(x= pca$x[rn,4], y= pca$x[rn,5], color=DATA[[col]][!is.na(DATA[[col]])], label=clonename))+geom_point(size=3) +geom_text(aes(label=clonename),hjust=0, vjust=0,size=1.2)+
    labs(title="PC4-PC5")+ xlab(paste("PC4: ",round(summ$importance[2,4]*100,digits=2),"%",sep="")) +ylab(paste("PC5: ",round(summ$importance[2,5]*100,digits=2),"%",sep=""))+ ggtitle(paste("Colored By:",col)) + 
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level3,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level2,linetype = "dotdash",colour="darkgrey") +
    stat_ellipse(aes(x = pca$x[rn,4],y=pca$x[rn,5]),inherit.aes=F,type=type,level=level1,linetype = "dotdash",colour="darkgrey")  
    midpoint=mean(c(min(DATA[[col]][!is.na(DATA[[col]])]),max(DATA[[col]][!is.na(DATA[[col]])])))
    max=max(DATA[[col]][!is.na(DATA[[col]])])
    min=min(DATA[[col]][!is.na(DATA[[col]])])
  d = d + scale_color_gradientn(name =eval(col),colors=c("cyan","orangered1"),labels=c(max,min),breaks=c(max,min)) + labs(col=substring(col,1,8))
  #theme(legend.key.size=unit(0.5,"cm")) + guides(colour = guide_legend(override.aes = list(size=2))) +
  #theme(legend.key.width=unit(0.2,"cm"))
}
# show(a)
# show(b)
# show(c)
# show(d)
multiplot_same_legend(a,b,c,d,cols=2) 
count=count+1
}

dev.off()

#   outliers_3SD_PC1_PC2
# }

outliers_3SD_PC1_PC2
