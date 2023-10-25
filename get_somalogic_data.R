.libPaths( c( .libPaths(), "/hpc/users/jayarp02/R/lib/4.0.3/") )
.libPaths()

library("dplyr")
library("stringr")
library("lubridate")

setwd(getwd())

### Load somalogic datas
#source("/sc/arion/projects/EHR_ML/pushkala_work/AKI_COVID_WGCNA/tag_aki_covid.R")
#write.csv(file = "somalogic_input_df.csv", somalogic)
somalogic = readRDS("somalogic_inputdf.rds")

#print(data.frame(colnames(somalogic)))
#lymphocyte_table = select(somalogic, starts_with("LAB"))
#print(colnames(lymphocyte_table))
#rnaseq_table = select(somalogic, contains("Viral"))
#print(colnames(rnaseq_table))

#print(rnaseq_table$RNA_Corrected_Sample[!is.na(rnaseq_table$RNA_Corrected_Sample)][1:500])

### Select only longest encounter for which there is a somalogic measurement
###MAIN DATAFRAME 
test <- somalogic  %>% group_by(Subject_ID) %>%
  filter(Encounter_Date %in% Encounter_Date[!is.na(RNA_Corrected_Sample)])

dim(test)
#print("checking to see if pediatric samples were removed")
#test2003 = data.frame(test[test$Subject_ID=="PICR2003", ])
#test2003$Encounter_Date[1:50]
#print("somalogic encounter date for picr2003")
#som2003 = data.frame(somalogic[somalogic$Subject_ID=="PICR2003", ])
#unique(som2003$Encounter_Date)[1:100]
#quit(status=1)
### Get covid ever positive during encounter
print("covid ever positive during encounter")

##USE THIS FOR COVID
###MAIN DATAFRAME
test <- test%>%group_by(Subject_ID, Encounter_Date) %>%
  mutate(
	COVID19_Ever_Positive_during_encounter =any(COVID19_Order_Result %in% c("DETECTED","PRESUMPTIVE POSITIVE"), na.rm = T)) 



#test <- test[test$COVID19_Ever_Positive_during_encounter==T,]
print("num unique subject IDs")
length(unique(test$Subject_ID))


#print("number encounters unique to each subject. the number totals upto the total number of unique subjects. ")
#num_encounters2 <- test %>% group_by(Subject_ID) %>%
#  summarise(n = length(unique(Encounter_Date)))
#num_encounters2 <- as.data.frame(num_encounters2)
#table(num_encounters2$n)
#num_encounters2[1:20,]

#table(num_encounters$n, num_encounters$COVID19_Ever_Positive_during_encounter)


### OLD - Select test longest encounter for each person
### NOW - If there are two encounters, then select the earliest one
print("get longest encounter date - to maximize number of time points")
###MAIN DATAFRAME
test <- test  %>% group_by(Subject_ID, Encounter_Date) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  group_by(Subject_ID) %>%
  filter(Encounter_Date == min(Encounter_Date))

test <- test %>% group_by(Subject_ID) %>%
  filter(!is.na(aki)) %>%
  bind_rows(
	test %>% group_by(Subject_ID) %>%
	filter(is.na(aki)) %>%
	mutate(aki_stage = ifelse(is.na(aki), 0, 1))
  )

#data.frame(test[test$Subject_ID=="PICR3009",])
#x <- test[test$Subject_ID %in% c("PICR7025", "PICR3005", "PICR2009", "PICR3004"),]
#x_df = data.frame(x$Subject_ID, x$aki, x$aki_stage, x$Days_Since_Admission, x$Event_Date, x$Specimen_RNAseq_Performed, x$RNA_Corrected_Sample, x$Encounter_Date, x$baseline, x$LAB_SERUM_CREATININE_Max, x$LAB_SERUM_CREATININE_Min, x$LAB_SERUM_CREATININE_Median, x$COMORBID_CHRONIC_KIDNEY_DISEASE, x$ENCOUNTER_COMORBID_ACUTE_KIDNEY_INJURY, x$RRT_Flag)
#x_df

print("****PUSH analysis****")
print("Num covid positive per subject per ncounter")

num_covid_pos_per_subject = data.frame(test$Subject_ID, test$Event_Date, test$Encounter_Date, test$aki, test$aki_stage, test$COVID19_Ever_Positive_during_encounter, test$RNA_Corrected_Sample, test$DECEASED_INDICATOR)
#head(num_covid_pos_per_subject)
print("test num covid pos per subject")
num_covid_pos_per_subject = num_covid_pos_per_subject %>%
	group_by(test.Subject_ID) %>%
        #filter(test.Subject_ID %in% c("PICR5041", "PICR5057")) %>%
	# filter(str_detect(test.Subject_ID, 'PICR3')) %>%
	mutate(n_encounter = length(unique(test.Encounter_Date)))%>%
	ungroup() %>%
	group_by(test.Subject_ID, test.Encounter_Date) %>%
        mutate(num_positive_inEnc = sum(test.COVID19_Ever_Positive_during_encounter == T, na.rm = T), 
		num_rnaseq_inEnc = sum(!is.na(test.RNA_Corrected_Sample), na.rm = T),
		num_aki_inEnc = sum(test.aki == T, na.rm = T),
		num_compositeaki_inEnc = sum((max(test.aki_stage) %in% c(2,3)) == T, na.rm = T),
		len_encounter = n()
		) 

print("order by desc n_encounter")
num_covid_pos_per_subject = num_covid_pos_per_subject[order(num_covid_pos_per_subject$n_encounter, decreasing = T),]
print("select columns")
num_covid_pos_per_subject = num_covid_pos_per_subject %>% select(c('test.Subject_ID', 'test.Event_Date', 'test.Encounter_Date', 'test.aki', 'n_encounter', 'len_encounter', 'num_positive_inEnc', 'num_rnaseq_inEnc', 'num_aki_inEnc', 'num_compositeaki_inEnc', 'test.DECEASED_INDICATOR'))
print("dataframe to write into csv")
data_for_case_ctrl = unique(data.frame(num_covid_pos_per_subject))
write.csv(data_for_case_ctrl, "data/preprocessed_data/data_for_case_ctrl.csv")

print("getting number of encounters per subject")

### count number of encounter per person
print("num encounters per person")
num_encounters <- test %>% group_by(Subject_ID) %>%
  summarise(n = length(unique(Encounter_Date)))

table(num_encounters$n) ### There is only one encounter per person. PJ upate - did not choose longest encounter. 



print("number died")
died <- test[test$DECEASED_INDICATOR==T,]
unique(died$Subject_ID)


### get dates when somalogic was taken
print("RNaseq samples - days since admission time points")
rnaseq_timepoints <- test[!is.na(test$Time_Point) &!is.na(test$RNA_Corrected_Sample),]
table(rnaseq_timepoints$Days_Since_Admission)

### check number of timepoints per person
print("num timepoints per person")
rnaseq_num_timepoints <- rnaseq_timepoints %>% group_by(Subject_ID)%>%
  summarise(count=n() )
rnaseq_num_timepoints <- as.data.frame(rnaseq_num_timepoints)
rnaseq_num_timepoints <- rnaseq_num_timepoints[order(rnaseq_num_timepoints$count, decreasing = T),]
table(rnaseq_num_timepoints$count)

print("severity score")
###MAIN DATAFRAME
test$my_severity_score <- test$Severity
test$my_severity_score[test$my_severity_score=="Moderate COVID-19"] <- 1
test$my_severity_score[test$my_severity_score=="Severe COVID-19"] <- 2
test$my_severity_score[test$my_severity_score=="Severe COVID-19 with EOD"] <- 3

test$my_severity_score <- as.numeric(test$my_severity_score)
table(test$my_severity_score)

##get table:
#serumcreat_table = select(test, contains("VENTI"))
#colnames(serumcreat_table)

### test ventilation type
###MAIN DATAFRAME
table(test$VENTILATION_TYPE)
test$ventilation <- test$VENTILATION_TYPE
test$ventilation[test$ventilation=="INVASIVE_VENTILATION"] <-3
test$ventilation[test$ventilation=="NON_INVASIVE_VENTILATION"] <-2
test$ventilation[test$ventilation=="SUPPLEMENTAL_O2"] <-1
test$ventilation[test$ventilation=="NONE_ROOM_AIR"] <-0
table(test$ventilation)


##get table:
#serumcreat_table = select(test, contains("CREAT"))
#colnames(serumcreat_table)

##print lab serum creatinine max scopre before modifier
#table(test$LAB_SERUM_CREATININE_Max)


### modify sofa score to remove creatinine contribution
###MAIN DATAFRAME
test$modifier <- case_when(
  test$LAB_SERUM_CREATININE_Max<1.2~ 0, 
  between(test$LAB_SERUM_CREATININE_Max, 1.2, 1.9)~ 1,
  between(test$LAB_SERUM_CREATININE_Max, 2.0, 3.4)~ 2, 
  between(test$LAB_SERUM_CREATININE_Max, 3.5, 4.9)~ 3, 
  test$LAB_SERUM_CREATININE_Max>4.9~4
)


print("number of subjects before modifier is modifier. and modifier table")
print(length(unique(test$Subject_ID)))
table(test$modifier)

test$modifier[is.na(test$modifier)] <- 0
table(is.na(test$modifier), is.na(test$SOFA_Score))

###MAIN DATAFRAME
test$modified_sofa <- test$SOFA_Score- test$modifier

print("summary of samples with rnaseq performed per individual")
specimen_rnaseq <- test%>% group_by(Subject_ID) %>%  
	summarize(non_na_count = sum(!is.na(RNA_Corrected_Sample), na.rm = T))

specimen_rnaseq = data.frame(specimen_rnaseq)
specimen_rnaseq


#rnaseq <- test%>% group_by(Subject_ID) %>%  summarize(non_na_count = sum(!is.na(RNA_Corrected_Sample)))
#rnaseq[30:100,]

print("comorbid CKD")
table(test$COMORBID_CHRONIC_KIDNEY_DISEASE)

print("comorbid aki")

#comorbid_table = select(test, contains("KIDNEY"))
#colnames(comorbid_table)
#table(test$ENCOUNTER_COMORBID_ACUTE_KIDNEY_INJURY)


#print("test event_date aki date")
#aki_rnaseq = test %>% group_by(Subject_ID)%>%
#            filter(!is.na(RNA_Corrected_Sample)& grepl("PICR2009", Subject_ID))%>%
#            select(c(Subject_ID, Event_Date, Time_Point, Specimen_RNAseq_Performed, aki, aki_stage, SOFA_Score, modified_sofa, my_severity_score, RNA_Corrected_Sample, Days_Since_Admission, Days_Until_Discharge))
#print(data.frame(aki_rnaseq))


##get dialysis table
##get table:
dialysis_table = select(test, contains("RRT"))
colnames(dialysis_table)

###Main TABLE!!
print("get aki encounter")
print("test!!")
###MAIN DATAFRAME
aki_encounter_date <- test %>% group_by(Subject_ID)%>%
  #filter(!str_detect(RNA_Corrected_Sample, "Pediatric")) %>%
  summarise(Subject_ID = min(Subject_ID),
            baseline_cr = min(baseline),
            Encounter_date= min(Encounter_Date),
	    max_event_date_rnaseqsamples = max(Event_Date[!is.na(RNA_Corrected_Sample)], na.rm = T),
            min_event_date_rnaseqsamples = min(Event_Date[!is.na(RNA_Corrected_Sample)], na.rm = T),
            rnaseq_first_date = min_event_date_rnaseqsamples,
	    rnaseq_last_date = max_event_date_rnaseqsamples,
            somalogic_date = min(Event_Date[Specimen_SomaLogics_Performed==T], na.rm = T),
            chosen_somalogic_date= max(Event_Date[!is.na(RNA_Corrected_Sample)], na.rm = T),
	    max_aki_stage = max(aki_stage[Event_Date<=rnaseq_last_date], na.rm = T),,
            aki_date=min(Event_Date[aki==T &!is.na(aki) ], na.rm = T),
            aki_last_date = max(Event_Date[aki==T &!is.na(aki)], na.rm = T),
            ####maybe try and find how many rnaseq performed during aki/overlap with aki?
            aki_hospital_day = min(Days_Since_Admission[aki==T &!is.na(aki) &aki_stage%in% c(2,3)], na.rm = T),
            hospital_stay = max(Days_Until_Discharge),
            ### earliest somalogic sample point
            rnaseq_sample = first(na.omit(RNA_Corrected_Sample[between(Event_Date, rnaseq_first_date, rnaseq_first_date) & !grepl("Pediatric", RNA_Corrected_Sample)])),
	    Uncorrected_Blood_Sample = min(Uncorrected_Blood_Sample[Event_Date==min_event_date_rnaseqsamples], na.rm=T),
            ###latest time point
            #somalogic_hospital_day=min(Days_Since_Admission[Event_Date==chosen_somalogic_date], na.rm=T),
            somalogic_hospital_day= max(Days_Since_Admission[Event_Date==max_event_date_rnaseqsamples], na.rm=T),
            rnaseq_hospital_day_min = min(Days_Since_Admission[Event_Date==min_event_date_rnaseqsamples], na.rm=T),
            rnaseq_hospital_day_max = max(Days_Since_Admission[Event_Date==max_event_date_rnaseqsamples], na.rm=T),
            severity = max(my_severity_score[Encounter_Date==min(Encounter_Date, na.rm = T) &between(Event_Date, max_event_date_rnaseqsamples-3, max_event_date_rnaseqsamples)] , na.rm=T),
            ventilation = max(ventilation[Encounter_Date==min(Encounter_Date, na.rm = T) &Event_Date <= max_event_date_rnaseqsamples]  , na.rm=T),
            sofa = max(SOFA_Score[between(Event_Date, max(Event_Date[!is.na(RNA_Corrected_Sample)], na.rm = T)-1, max_event_date_rnaseqsamples)]  , na.rm=T),
            #earliest_sofa_min = SOFA_Score[Encounter_Date==min(Encounter_Date, na.rm = T) &Event_Date== min(Event_Date, na.rm = T)],
            modified_sofa = max(modified_sofa[between(Event_Date, max_event_date_rnaseqsamples-1, max_event_date_rnaseqsamples)]  , na.rm=T),
            #earliest_sofa_hospital_day_min = Days_Since_Admission[Encounter_Date==min(Encounter_Date, na.rm = T) &Event_Date== min(Event_Date, na.rm = T)],
            #earliest_modified_sofa_min = modified_sofa[Encounter_Date==min(Encounter_Date, na.rm = T) &Event_Date== min(Event_Date, na.rm = T)],
            
            icu = max(ICU_Status[Event_Date==max_event_date_rnaseqsamples], na.rm=T),
            min_platelet = max(LAB_PLATELET_Min[Event_Date==max_event_date_rnaseqsamples], na.rm=T),
            max_bun = max(LAB_BUN_Max[Event_Date==max_event_date_rnaseqsamples], na.rm=T),
            min_lymphocyte = max(LAB_LYMPHOCYTE_NUMBER_Min[Event_Date==max_event_date_rnaseqsamples], na.rm=T),
            max_il1beta = max(LAB_INTERLEUKIN_1_BETA_Max[Event_Date==max_event_date_rnaseqsamples], na.rm=T),
            max_cr = max(LAB_SERUM_CREATININE_Max[Event_Date==max_event_date_rnaseqsamples], na.rm=T),

            ##flags
            aki = ifelse(any(aki, na.rm = T), 1,0),
            num_rnaseq = sum(!is.na(RNA_Corrected_Sample), na.rm = T),
            death = ifelse(any(DECEASED_INDICATOR, na.rm = T), 1,0),
            ckd = ifelse(any(COMORBID_CHRONIC_KIDNEY_DISEASE==T), 1, 0),
            comorbid_aki = ifelse(any(ENCOUNTER_COMORBID_ACUTE_KIDNEY_INJURY == T), 1, 0),
            composite_outcome = ifelse(max_aki_stage %in% c(2,3) , 1, 0),

            ###PJ for now put all dialysis to 1
            dialysis = ifelse(any(RRT_Flag==1, na.rm = T), 1,0),
            ###this one is causing problems!
	    covid_most_recent = max(COVID19_Most_Recent_Order_Result[Event_Date==max_event_date_rnaseqsamples], na.rm=T),
            covid_ever_positive_during_encounter = ifelse(any(COVID19_Order_Result[between(Event_Date, min_event_date_rnaseqsamples, max_event_date_rnaseqsamples)] %in% c("DETECTED","PRESUMPTIVE POSITIVE"), na.rm = T), 1, 0),
            #admission_viral_load = if(all(is.na(Viral_Load_N2_per_mL))) 0 else first(na.omit(Viral_Load_N2_per_mL)),
            #admission_viral_load_log = log10(admission_viral_load+10),
            #viral_load_date = min(Event_Date[Viral_Load_N2_per_mL==admission_viral_load & !is.na(Viral_Load_N2_per_mL)], na.rm = T),
            max_viral_load = ifelse (all(is.na(Viral_Load_N2_per_mL)), 0, max(na.omit(Viral_Load_N2_per_mL))),
            max_viral_load_log = log10(max_viral_load+10),
            max_viral_load_date = min(Event_Date[Viral_Load_N2_per_mL==max_viral_load & !is.na(Viral_Load_N2_per_mL)], na.rm = T)

	)

print("done test!")

aki_encounter_date <- as.data.frame(aki_encounter_date)
print("raw aki encounter_date")
#aki_encounter_date[aki_encounter_date$Subject_ID=="PICR3005",]
#aki_encounter_date
dim(aki_encounter_date)

print("aki encounter date earliest. remove NA")
head(aki_encounter_date[is.na(aki_encounter_date$Uncorrected_Blood_Sample),])

print("aki encounter date again")
aki_encounter_date <- aki_encounter_date[!is.na(aki_encounter_date$Uncorrected_Blood_Sample),]
dim(aki_encounter_date)

print("table aki encounter date for aki stats")
print("aki")
table(aki_encounter_date$aki)
print("aki_date")
table(aki_encounter_date$aki_date)
print("chosen_somalogic_date")
table(aki_encounter_date$chosen_somalogic_date)
aki_encounter_date$aki_date_minus_somalogic_date <- aki_encounter_date$aki_date-aki_encounter_date$chosen_somalogic_date
aki_encounter_date$aki <- ifelse(aki_encounter_date$aki_date_minus_somalogic_date==Inf, 0, 1)

print("min aki date minus max date that rnaseq was done")
table(aki_encounter_date$aki_date_minus_somalogic_date)

print("table aki again")
table(aki_encounter_date$aki)

print("severity")
table(aki_encounter_date$severity)

print("cs stats for aki and aki encounter date where difference between encounter date and somalogics/rnaseq date is less than 0")
table(aki_encounter_date$aki, aki_encounter_date$aki_date_minus_somalogic_date<0)

print("stats for difference between encounter date and somalogics/rnaseq date")
table(aki_encounter_date$aki_date_minus_somalogic_date)

print("stats for somalogic hospital day - min of days since admission for when rnaseq was perfomed")
table(aki_encounter_date$somalogic_hospital_day)

print("stats for aki hospital day - min number of days since admission for which aki = T")
table(aki_encounter_date$aki_hospital_day)

print("stats for aki encounter date where aki flag shows something")
table(aki_encounter_date$aki)

print("num unique rows for subject ID - i.e. num patients in this table")
length(unique(aki_encounter_date$Subject_ID))

print("aki encounter date where aki date minus somalogic date is Inf - did not have AKI or it should have somthing done within 3 days")
#aki_encounter_date <- aki_encounter_date[aki_encounter_date$aki_date_minus_somalogic_date=="Inf" | aki_encounter_date$aki_date_minus_somalogic_date<3 | aki_encounter_date$aki_date_minus_somalogic_date>-3, ]
aki_encounter_date <- aki_encounter_date[aki_encounter_date$aki_date_minus_somalogic_date==Inf |(aki_encounter_date$aki_date_minus_somalogic_date> -3 &aki_encounter_date$aki_date_minus_somalogic_date< 3),]
print("last dim aki encounter date")
dim(aki_encounter_date)
print("severity table")
table(aki_encounter_date$severity)
print("aki table")
table(aki_encounter_date$aki)
print("sofa table")
table(aki_encounter_date$sofa)
print("cs stats covid positive during encounter and aki table")
table(aki_encounter_date$covid_ever_positive_during_encounter,aki_encounter_date$aki )
print("dialysis table")
table(aki_encounter_date$dialysis)
print("cs stats aki and max aki stage table")
table(aki_encounter_date$aki, aki_encounter_date$max_aki_stage)

print("find index of duplicated records..?")
aki_encounter_date[duplicated(aki_encounter_date$Subject_ID),]
#quit(status=1)
#print("looking at sample PICR7118")
#aki_encounter_date[aki_encounter_date$Subject_ID=="PICR7118",]


print("remove duplicated elements")
aki_encounter_date <-aki_encounter_date[!duplicated(aki_encounter_date$Subject_ID),]
dim(aki_encounter_date)
print("-------get eskd-------")

eskd_dx_codes <- c("ICD10_N18.6")
dialysis_codes <- c("ICD10_Z99.2", "ICD10_T85.691S", "ICD10_T85.691A","ICD10_T85.71XD")
dx_subset <- test[match(aki_encounter_date$Subject_ID, test$Subject_ID), c("Subject_ID", "ICD10_N18.6")]
dx_subset <- as.data.frame(dx_subset)
eskd_dx_results <- data.frame(Subject_ID= dx_subset$Subject_ID, eskd_code=sapply(dx_subset[,-1],function(x) any(!is.na(x) &x!=0)))

dialysis_subset <- test[match(aki_encounter_date$Subject_ID, test$Subject_ID), c("Subject_ID", "ICD10_Z99.2", "ICD10_T85.691S", "ICD10_T85.691A","ICD10_T85.71XD")]
dialysis_subset <- as.data.frame(dialysis_subset)
dialysis_results <- data.frame(Subject_ID= dialysis_subset$Subject_ID, dialysis_code=apply(dialysis_subset[,-1], MARGIN = 1, function(x) any(!is.na(x) &x!=0)))
combined <- merge(eskd_dx_results, dialysis_results, by="Subject_ID")
table(combined$eskd_code, combined$dialysis_code)
combined$eskd <- combined$eskd_code==T &combined$dialysis_code==T
#combined[combined$Subject_ID=="PICR5022", ]
print("dim before removing")
dim(aki_encounter_date)
#aki_encounter_date[aki_encounter_date$Subject_ID%in%c("PICR5022", "PICR7073", "PICR7157"), ]
### remove eskd patients
#load( file="data/preprocessed_data/eskd.RData")
aki_encounter_date$eskd <- combined$eskd[match(aki_encounter_date$Subject_ID, combined$Subject_ID)]
aki_encounter_date <- aki_encounter_date[aki_encounter_date$eskd==F, ]
print("sim after removing")
dim(aki_encounter_date)
save(aki_encounter_date, file="data/preprocessed_data/somalogic_rnaseq_after_aki_encounter_clinical_data.RData")
# save(aki_encounter_date, file="data/preprocessed_data/somalogic_before_severe_aki_encounter_clinical_data.RData")

print("those who died****")
died <- aki_encounter_date[aki_encounter_date$death==T,]
table(died$aki, died$covid_ever_positive_during_encounter)
notdied <- aki_encounter_date[aki_encounter_date$death==F,]
table(notdied$aki, notdied$covid_ever_positive_during_encounter)
print("those with rnaseq data")
rna <- aki_encounter_date[aki_encounter_date$num_rnaseq>0,]
table(rna$composite_outcome, rna$covid_ever_positive_during_encounter)

quit(status=1)

### Filter by cutoff date
print("filter by cutoff date Encounter date")
summary(test$Encounter_Date)
table(aki_encounter_date$aki_date_minus_somalogic_date )

cutoff_date <- median(aki_encounter_date$Encounter_date)
cutoff_date <- quantile(aki_encounter_date$Encounter_date, probs = 0.5, type = 1, na.rm=T)

# cutoff_date <- quantile(test$Encounter_Date ,0.7)
table(aki_encounter_date$Encounter_date ==cutoff_date)
table(aki_encounter_date$Encounter_date > cutoff_date)
table(aki_encounter_date$Encounter_date <= cutoff_date, aki_encounter_date$aki)

print("those who died****")
died <- aki_encounter_date[aki_encounter_date$death==T,]
table(died$aki, died$covid_ever_positive_during_encounter)
notdied <- aki_encounter_date[aki_encounter_date$death==F,]
table(notdied$aki, notdied$covid_ever_positive_during_encounter)
print("those with rnaseq data")
rna <- aki_encounter_date[aki_encounter_date$num_rnaseq>0,]
table(rna$composite_outcome, rna$covid_ever_positive_during_encounter)



quit(status=1)


### Split data into discovery/validation
encounter_dates <- unique(aki_encounter_date$Encounter_date)
encounter_dates <- encounter_dates[order(encounter_dates)]
validation <- seq(1, length(encounter_dates), 3)
validation_dates <- encounter_dates[validation]

discovery_data <- aki_encounter_date[!aki_encounter_date$Encounter_date %in% validation_dates,]
table(discovery_data$aki)
validation_data <- aki_encounter_date[aki_encounter_date$Encounter_date %in% validation_dates,]
table(validation_data$aki)
save(discovery_data, file="data/preprocessed_data/aki_discovery.RData")

save(validation_data, file="data/preprocessed_data/aki_validation.RData")


discovery_data <- aki_encounter_date[aki_encounter_date$Encounter_date <=cutoff_date,]
table(discovery_data$aki)
validation_data <- aki_encounter_date[aki_encounter_date$Encounter_date >cutoff_date,]
table(validation_data$aki)

save(discovery_data, file="data/preprocessed_data/aki_discovery.RData")

save(validation_data, file="data/preprocessed_data/aki_validation.RData")


discovery_data <- aki_encounter_date[aki_encounter_date$Encounter_date %in% discovery_dates,]
table(discovery_data$aki)
validation_data <- aki_encounter_date[!aki_encounter_date$Encounter_date %in% discovery_dates,]
table(validation_data$aki)

save(discovery_data, file="data/preprocessed_data/aki_discovery_composite.RData")

save(validation_data, file="data/preprocessed_data/aki_validation_composite.RData")


aki_encounter_date <- aki_encounter_date[aki_encounter_date$Encounter_date >cutoff_date,]
table(aki_encounter_date$death, aki_encounter_date$max_aki_stage)
table(aki_encounter_date$composite_outcome )

save(aki_encounter_date, file="data/preprocessed_data/aki_discovery_composite.RData")

save(aki_encounter_date, file="data/preprocessed_data/aki_validation_composite.RData")
length(unique(aki_encounter_date$Subject_ID))

#table(aki_encounter_date$c)
table(aki_encounter_date$ckd)
table(aki_encounter_date$severity)
table(aki_encounter_date$somalogic_hospital_day)
### only take patients with somalogic within 48 hours of hospital admissionn
# aki_encounter_date <- aki_encounter_date[aki_encounter_date$somalogic_hospital_day<=2,]
# table(aki_encounter_date$aki)



table(aki_encounter_date$aki)
table(aki_encounter_date$somalogic_hospital_day)

print("those who died****")
died <- aki_encounter_date[aki_encounter_date$death==T,]
died[1:150]



subset <- somalogic[somalogic$Subject_ID==died$Subject_ID[1],]
data.frame(subset$Encounter_Date, subset$Event_Date, subset$DECEASED_INDICATOR)


