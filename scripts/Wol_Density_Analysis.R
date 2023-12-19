# MNE ANALYSIS INTRO ------------------------------------------------------
# Analysing Wolbachia Density data ## 
# Created by: Angelique Asselin 
# Created: 10/08/2023
# Edited: 14/08/2023 
# Edited: 19/12/2023 


# EDIT LOG  ---------------------------------------------------------------

# 19/12/23 
# Updating the MNE and TCID50 with the clean stuff.

# PREPARATION -------------------------------------------------------------

# Writes today' date for the file names 
today <- Sys.Date()
date <- format(today, format = "%y%m%d")

# Installs packages needed (can be skipped if you have already downloaded them)
# install.packages("readxl")
# install.packages("tidyverse")
# install.packages("strex")
# install.packages("moments")


# Borrows packages from the library of those installed
library(readxl)
library(tidyverse)
library(strex)
library(ggplot2)
library(ggtext)
library(ggpubr)
library(moments)
library(MASS)
library(rstatix)
library(lmtest)
library(car)

# What are your genes
house_gene <- "rpl" # house keeping gene
target_gene <- "WOL" # target gene

# Imports log of all the samples to be tested for Wolbachia density 
# Only necessary if you have a spreadsheet with all your samples and their details. 
# The list of sample names 

# IMPORTING MNE AND TCID data  --------------------------------------------

# Importing data 
Log <- (read_csv("/Users/angeliqueasselin/Documents/GitHub/modelling/data/model/231213_Data_MNE_TCID_model.csv")) %>% 
  subset(Initial_V == 50) %>% # Selecting only the dataset for which density was measured
  subset(Fly_Line == "Wol") # Selecting only the samples that have Wolbachia

# IMPORTING RAW DATA ------------------------------------------------------

# Compiles a list of all the raw data files (exported from Qiagen Machine)
## CHANGE THIS TO BE THE DIRECTORY OF YOUR FOLDER ##
file.list_MNE <- list.files( path = "data/raw", pattern='*.csv', full.names = TRUE) # lists all the .csv files in the folder

# Creates a data frame for each file
df.list_MNE <- lapply(file.list_MNE, function (x){
  df <- read.csv(x, as.is = TRUE ,stringsAsFactors = T, sep = ',')
  df$Take.Off<- as.double(df$Take.Off)
  df$Amplification<- as.double(df$Amplification)
  df$X2nd.Deriv..Max<- as.double(df$X2nd.Deriv..Max)
  df$Rep..Takeoff<- as.double(df$Rep..Takeoff)
  df$Rep..Takeoff..Std..Dev..<- as.double(df$Rep..Takeoff..Std..Dev..)
  df$File.Name <- as.character(x)
  return(df)
}) # turns each file in the list into a dataframe

# Compiles each file data frame into one compiled data frame
df <- bind_rows(df.list_MNE)

# Checks the dataframe

tail(df)
head(df)

# TIDYING RAW DATA --------------------------------------------------------

# Gets rid of excess columns if present
df <- df[,1:7]

# Splits the names of the samples rpl_405 to rpl in one column and the sample number 405 in another 
df_new <- df %>%
  mutate(Number = str_extract(Name, "\\d+"),
         Letter = str_extract(Name, "[A-Za-z]+")) %>%
  mutate(Number = ifelse(is.na(Number), str_extract(Name, "(?<=_)[0-9]+"), Number),
         Letter = ifelse(is.na(Letter), str_extract(Name, "[A-Za-z]+"), Letter)) %>%
  select(-Name)

# Renaming all the Columns 
df_tidy_1 <- df_new %>% 
  rename(Gene = Letter, Sample = Number, File = File.Name ,  Take_off = Take.Off, Second_Der_Max = X2nd.Deriv..Max, Rep_Take_off = Rep..Takeoff, Std_dev = Rep..Takeoff..Std..Dev.. ) 

# Changing any formatting to be consistent rpl -> rpl 
df_tidy_2 <- df_tidy_1 %>% 
  mutate(Sample = as.numeric(Sample)) %>%  # Makes sure the samples are considered numeric 
  mutate(Sample_number_correct = ifelse(!is.na(Sample), TRUE, FALSE)) %>% # Checks that all the sample numbers are correct
  filter(Sample_number_correct == TRUE) %>% # Removes and NA samples (those NTC and other control samples)
  unite(Sample_file, Sample, File, sep = "_", remove = FALSE) %>% # Combining the sample and file name to make sure samples are paired with their run (rpl and WOL)
  select(-File) %>% # Removes just the file name column 
  mutate(Gene_formated = ifelse(Gene %in% c(house_gene, target_gene), TRUE, FALSE))%>% # codes when Gene is positive or neg controls
  dplyr::filter(Gene_formated == TRUE) # Filter out any NTC or NRT controls

# Creating a data frame of all the house keeping gene data
house_gene_data<- df_tidy_2 %>% 
  filter(Gene == house_gene) %>% 
  group_by(Sample_file) %>%
  mutate(Tech_Rep = row_number())

# Creating a data frame of all the target gene data
target_gene_data <- df_tidy_2 %>% 
  filter(Gene == target_gene) %>% 
  group_by(Sample_file) %>%
  mutate(Tech_Rep = row_number()) 


# Joining the samples so the gene take offs are across one row not two separate rows

joined_df <- full_join(house_gene_data, target_gene_data, by = c("Sample_file", "Tech_Rep"), keep = TRUE) 

# Re-formats to get rid of repeated information but also label according to the gene
house_gene_name <- paste0("_", house_gene)
target_gene_name <- paste0("_", target_gene)


joined_df_1 <-  joined_df %>%
  rename_with(~paste0(gsub("\\.x$", house_gene_name, .)),
              ends_with(".x")) %>%
  rename_with(~paste0(gsub("\\.y$", target_gene_name, .)),
              ends_with(".y")) %>% 
  mutate(Sample_file = ifelse(Sample_file_rpl == Sample_file_WOL, # creates Sample_file
                              Sample_file_rpl, 
                              NA))  %>% 
  mutate(Filter_Reason = ifelse(is.na(Sample_file),  "Unpaired", NA)) # adds in a section to mention if the sample is unpaired


# Making sure there are no NA in the data and that each technical replicate contains the correct sd and rep take off information 

joined_df_2 <- joined_df_1 %>%
  group_by(Sample_file) %>%
  mutate(Std_dev_rpl = ifelse(is.na(Std_dev_rpl), Std_dev_rpl[!is.na(Std_dev_rpl)], Std_dev_rpl)) %>%
  mutate(Std_dev_WOL = ifelse(is.na(Std_dev_WOL), Std_dev_WOL[!is.na(Std_dev_WOL)], Std_dev_WOL)) %>%
  mutate(Rep_Take_off_rpl = ifelse(is.na(Rep_Take_off_rpl), Rep_Take_off_rpl[!is.na(Rep_Take_off_rpl)], Rep_Take_off_rpl)) %>%
  mutate(Rep_Take_off_WOL = ifelse(is.na(Rep_Take_off_WOL), Rep_Take_off_WOL[!is.na(Rep_Take_off_WOL)], Rep_Take_off_WOL)) %>%
  ungroup() 

# QUALITY CONTROL RAW DATA 1.0 ---------------------------------------------------------

# Ensuring that the data meets a few quality control thresholds 
# Take off ct value for each gene is above minimu < 11.2 and below maximum >= 30 
# Std deviation between technical reps is < 0.5
# Amplification efficiency of each primer is > 1.4 

QC_data_1 <- joined_df_2 %>% 
  mutate(Filter_Reason = ifelse(!is.na(Take_off_rpl) & Take_off_rpl > 30, paste(Filter_Reason, "Take off rpl too high repeat with more concentrated sample", sep = " and "), Filter_Reason)) %>% 
  mutate(Filter_Reason = ifelse(!is.na(Take_off_rpl) & Take_off_rpl < 11.2 , paste(Filter_Reason, "Take off rpl too low repeat with more diluted sample", sep = " and "), Filter_Reason)) %>% 
  mutate(Filter_Reason = ifelse(!is.na(Take_off_WOL) & Take_off_WOL < 11.2 , paste(Filter_Reason, "Take off WOL too low repeat with more diluted sample", sep = " and "), Filter_Reason)) %>% 
  mutate(Filter_Reason = ifelse(!is.na(Take_off_WOL) & Take_off_WOL > 30 , paste(Filter_Reason, "Take off WOL too high repeat with more concentrated sample", sep = " and "), Filter_Reason)) %>% 
  mutate(Filter_Reason = ifelse(!is.na(Amplification_rpl) & Amplification_rpl < 1.4 , paste(Filter_Reason, "Amplification Efficiency (rpl)", sep = " and "), Filter_Reason)) %>% 
  mutate(Filter_Reason = ifelse(!is.na(Amplification_WOL) & Amplification_WOL < 1.4 , paste(Filter_Reason, "Amplification Efficiency (WOL)", sep = " and "), Filter_Reason)) %>% 
  mutate(Filter_Reason = ifelse(!is.na(Std_dev_rpl) & Std_dev_rpl >= 0.5 , paste(Filter_Reason, "Standard Deviation (rpl)", sep = " and "), Filter_Reason)) %>% 
  mutate(Filter_Reason = ifelse(!is.na(Std_dev_WOL) & Std_dev_WOL >= 0.5 , paste(Filter_Reason, "Standard Deviation (WOL)", sep = " and "), Filter_Reason)) %>% 
  mutate(Quality_Control = ifelse(is.na(Filter_Reason), "PASS", "FAIL")) %>%
  group_by(Sample_file) %>% 
  mutate(Filter_Reason = ifelse(Quality_Control[1] == Quality_Control[2], Filter_Reason, paste(Filter_Reason, "Only one replicate passed QC", sep = " and "))) %>% 
  ungroup() %>% 
  mutate(Quality_Control = ifelse(is.na(Filter_Reason), "PASS", "FAIL")) %>% # updating the pass and fail 
  mutate(Sample_rpl = as.numeric(Sample_rpl)) %>% 
  mutate(Sample_WOL = as.numeric(Sample_WOL)) %>% 
  mutate(Sample = ifelse(Sample_rpl == Sample_WOL, Sample_rpl, NA)) %>% 
  select(-Sample_rpl, -Sample_WOL) 


# FORMAT RAW DATA FOR ANALYSIS --------------------------------------------

# Creates data frame of only data that have passed quality control QC

data  <- QC_data_1 %>% 
  filter(is.na(Filter_Reason))

# Formats the QC data frame for analysis

data_analyse <- data %>% 
  select(Take_off_rpl, Take_off_WOL, Sample_file) %>%
  rename(Sample_Name = Sample_file, Reference_Takeoff = Take_off_rpl, Target_Takeoff = Take_off_WOL) %>%  
  mutate(Number = row_number()) %>% 
  select(Number, Sample_Name, Target_Takeoff, Reference_Takeoff) # puts in the correct order


# ANALYSIS ----------------------------------------------------------------

# Any problems let me know and I can sort it! - Tyson
# This analysis utilises Equation (2) and (3) from DOI 10.1093/bioinformatics/btg157
# Before using: This is for Simplex PCR samples (reference and target primers run in different tubes)
# For Multiplex PCR samples, use other Script
# NOTE THIS WILL STILL ANALYSE DATA IF YOU HAVE A DIFFERENT NUMBER OF REPS 

Data <- data_analyse

#Sets % SEM cutoff for samples (default 20%)
SEMCutoff = 20

# Do you wish to automatically censor data over the above SEM% threshold? 0 = no, 1 = yes (default yes)
UseCutoff = 1

# Change std Amp. efficiencies if necessary
TarAmp = 1.873 # WOL
RefAmp = 1.64 # rpl32

# If you need to delete the dataframe \/ (do at start of each analysis)
rm(Normalised)
rm(list)

# Run me (housekeeping stuff to sort the data/prepare for analysis)
Normalised <- data.frame("Sample","MNE","SE of MNE", "SE of MNE as % (x<=20%)")
i=0
val = 1:nrow(Data)
val <- val[seq(1,(length(val)), 2)]
w = unlist(c(Data[4]))
x = unlist(c(Data[3]))
y = unlist(c(Data[2]))
s = matrix(x,nrow=1,ncol=nrow(Data))
v = matrix(w,nrow=1,ncol=nrow(Data))
t = matrix(y,nrow=nrow(Data),ncol=1)
list <- data.frame("Censored Samples","MNE","SE of MNE", "SE of MNE as % (x>25%)")

# For loop processes all the data and adds it to the dataframe
for (i in val) {
  
  # Calculates MNE of Duplicates
  MNE2 = (RefAmp^((v[1,i]+v[1,i+1])/2))/(TarAmp^((s[1,i]+s[1,i+1])/2))
  
  # Calculates SE of MNE
  SEofCTTAR = sd(c(s[1,i],s[1,i+1]))/sqrt(2)
  SEofCTREF = sd(c(v[1,i],v[1,i+1]))/sqrt(2)
  
  SEofMNE = MNE2*(((log(TarAmp)*SEofCTTAR)^2)+((log(RefAmp)*SEofCTREF)^2))^.5
  
  # Calculates SE of MNE as %
  SEofMNEasP = (SEofMNE/MNE2)*100
  
  # Adds each set of data to the data frame (you may wish to change the digits= option to be smaller if you wish)
  if (UseCutoff == 1) {
    if (SEofMNEasP > SEMCutoff) {
      cat("Sample:",t[i,1],"was censored due to exceeding SEM cutoff of",SEMCutoff, "with a value of",(round(SEofMNEasP,digits=2)),"\n")
      list[nrow(list) + 1,] = c(t[i,1], round(MNE2,digits=10), round(SEofMNE,digits=5), round(SEofMNEasP,digits=2)) }
    else {
      Normalised[nrow(Normalised) + 1,] = c(t[i,1], round(MNE2,digits=10), round(SEofMNE,digits=5), round(SEofMNEasP,digits=2))}
  }
  else {
    Normalised[nrow(Normalised) + 1,] = c(t[i,1], round(MNE2,digits=10), round(SEofMNE,digits=5), round(SEofMNEasP,digits=2))}
}

# Formats the analysed dataset 

MNE <- Normalised %>% 
  rename(Sample_file = X.Sample., WOL_MNE = X.MNE., WOL_SE_MNE = X.SE.of.MNE., WOL_SE_Percent = X.SE.of.MNE.as....x..20...) %>% 
  slice(-1) %>% 
  mutate(Sample = str_first_number(Sample_file)) %>% # Formats to include sample number
  select(Sample, everything())

# Visualize the data frame

MNE

# Formats the censored data set
SEM_censored <- list %>% 
  rename(Sample_file = X.Censored.Samples., MNE = X.MNE., SE_MNE = X.SE.of.MNE., SE_Percent = X.SE.of.MNE.as....x.25...) %>% 
  mutate(Sample_file = as.character(Sample_file)) %>% 
  slice(-1)

# Visualize the data set
SEM_censored

# QUALITY CONTROL ANALYSED DATA 2.0 -----------------------------------------------------

# Creates a data frames adding filter reason for SEM

QC_data_2 <- QC_data_1 %>% # Data set before all these had been filtered out 
  mutate(Filter_Reason = ifelse(Sample_file %in% SEM_censored$Sample_file, paste(Filter_Reason, "SEM", sep = " and "), Filter_Reason)) %>% # adds filter reason for SEM 
  mutate(Quality_Control = ifelse(is.na(Filter_Reason), "PASS", "FAIL"))

# REPEATS -----------------------------------------------------------------

# Looks at all the samples that need to be done and then determines which ones don't
# have Wolbachia density data. These ones need to be completed/repeated. 
# If they need to be repeated this will determine the appropriate approach
# This is based on why the samples failed QC

# Creates a vector of all the samples for the experiment 

sample <- Log %>% # from the data log
  select(Sample) %>% 
  mutate(Progress = "TO BE DONE") %>% 
  mutate(Sample = as.numeric(Sample)) %>% 
  as.data.frame() %>% 
  arrange(Sample)

# Create a vector of unique sample names that have been analysed

COMPLETE_samples <- as.data.frame(unique(MNE$Sample)) %>% 
  rename(COMPLETE = "unique(MNE$Sample)") %>% 
  separate(COMPLETE, into = c("Sample", "File"), sep = "_") %>%  #This will prompt a warning but it is fine it is just getting rid of the _raw and .csv following the file 
  select(-File) %>% 
  mutate(Progress = "COMPLETE") %>% 
  mutate(Sample = as.numeric(Sample)) %>% 
  arrange(Sample) # arrange in ascending order based on sample number

# Joining the two data frames together and any that are missing should have an NA in the Sample 

joined_repeats <- left_join(x = sample, y = COMPLETE_samples, by = "Sample", keep = TRUE) %>% # combine the two by matching up the sample if the sample doesn't match it will have NA but keep the data
  mutate(Progress = if_else(!is.na(Progress.y), "COMPLETE", "TO BE DONE")) %>% # where y = MNE data samples that have been COMPLETED
  select(-Progress.x, -Progress.y) %>% 
  mutate(Sample =  Sample.x) %>% 
  select(-Sample.x, -Sample.y) %>% 
  rename(Sample = "Sample") %>% 
  mutate(Sample = as.numeric(Sample)) %>% 
  select(Sample, Progress) %>% 
  arrange(Sample)

# Selecting only those that still need to be done 

repeats <- joined_repeats %>% 
  filter(Progress == "TO BE DONE")

# Using the repeats data frame to filter out the reason for the repeats

repeats_df <- left_join(x = repeats, y = QC_data_2, by = "Sample") %>% # left join because you only want those that need to be repeated
  mutate(Filter_Reason = ifelse(is.na(Filter_Reason), yes = paste(Filter_Reason, "No raw data", sep = " and "), no = Filter_Reason))

# Using the filter reason to determine the appropriate repeat approach

repeats_full_1 <- repeats_df %>% 
  mutate(Approach = ifelse(str_detect(Filter_Reason, "dilute"), yes =  "qPCR Dilute", no = NA)) %>% 
  mutate(Approach = ifelse(is.na(Approach) & str_detect(Filter_Reason, "No raw data"),yes =  "qPCR Original", no = Approach)) %>% 
  mutate(Approach = ifelse(is.na(Approach), "Check", Approach)) %>% 
  arrange(Sample, Approach) %>% 
  separate(Sample_file, c("Duplicate_sample", "File"), sep = "_") %>% 
  select(-Duplicate_sample) %>% 
  arrange(Approach)


# Sub-setting the repeats from the log to get all the associated data with them 
repeats_full_2 <- left_join(x = repeats_full_1, y = Log, by = join_by(Sample)) %>% 
  select(Sample, Approach, Filter_Reason, everything()) %>% 
  arrange(Sample)

# Creates a sub set of the data frame that is just those that need to be checked

repeats_check <- repeats_full_1 %>% 
  filter(Approach == "Check")

# Wolbachia Status 
wol_stat <- repeats_full_2 %>% 
  mutate(Wol_status = ifelse(Approach == "qPCR Original", yes = "M", no= NA)) %>% 
  mutate(Wol_status = ifelse(is.na(Wol_status) & Approach == "Check", yes = "ND", Wol_status)) %>% 
  select(Sample, Approach, Wol_status) %>% 
  unique()

wol_nd <- wol_stat %>% 
  subset(Wol_status == "ND")

wol_M <- wol_stat %>% 
  subset(Wol_status == "M")


# FULL DATASET (with new data) -----------------------------------------------------------
data_final_graph <- full_join(Log, MNE, by = join_by(Sample)) %>%
  mutate(Wol_status = ifelse(Sample %in% wol_M$Sample, "M", NA)) %>%
  mutate(Wol_status = ifelse (Sample %in% wol_nd$Sample, "ND", Wol_status)) %>%
  mutate(
    Time_group = case_when(
      Time %in% c(0, 6, 10, 24, 34, 48) ~ 1,
      Time %in% c(58, 72, 82, 96, 106, 120) ~ 2,
      Time %in% c(130, 144, 154, 168, 178, 192) ~ 3,
      Time %in% c(202, 216, 226, 240, 250, 264) ~ 4,
      Time %in% c(288, 312, 336, 360, 384, 408) ~ 5,
      TRUE ~ NA_integer_
    )
  ) %>%
  arrange(Time) %>%
  mutate(Time_group = as.character(Time_group)) %>%
  mutate(WOL_MNE = as.numeric(WOL_MNE)) %>% 
  mutate(log_WOL_MNE = log10(WOL_MNE))%>%
  mutate(log_IU_fly = as.numeric(log_IU_fly)) %>% 
  mutate(log_MNE = as.numeric(log_MNE)) %>% 
  select(Sample,
         Fly_Line,
         Titre,
         Initial_V,
         Time,
         Time_group,
         log_IU_fly,
         log_MNE,
         log_WOL_MNE, 
         everything()) 

data_final_stat <- data_final_graph %>% 
  mutate(Time = as.numeric(Time))

# Calculating the sample number (not including the NA samples)
sample_IU <- data_final_stat %>% 
  drop_na(log_IU_fly, WOL_MNE)

sample_number_IU <- length(sample_IU$log_IU_fly)

sample_RNA <- data_final_stat %>% 
  drop_na(log_MNE, WOL_MNE)

sample_number_RNA <- length(sample_RNA$log_MNE)



# PRESENTING DATA  --------------------------------------------------------


# Graphing WOl and the virus particles (IU)
graph_IU <- ggplot(data_final_graph) +
  geom_point(aes(x = WOL_MNE, y = log_IU_fly, col = as.factor(Time_group)), size = 5) +  # Adjust size here
  geom_smooth(aes(x = WOL_MNE, y = log_IU_fly), method = "glm") +
  scale_color_manual(
    values = c("lightblue", "blue", "pink", "orange", "red"),
    labels = c("24-48 Hours", "58-120 Hours", "130-192 Hours", "202-264 Hours", "288-408 Hours"), 
    name = "Time"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 35),
    axis.title.x = element_markdown(),
    axis.title.y = element_text(size = 35, face = "bold"),
    axis.text = element_text(size = 30),
    legend.text = element_text(size = 20),  # Adjust legend text size
    legend.title = element_text(size = 25)  # Adjust legend title size
  ) +
  ylab(expression(paste(log[10], "DCV " , "Infectious " , "Virus ", "Particles "))) +
  xlab(expression("Wolbachia Density")) +
  guides(fill = guide_legend(title = "Time", label.theme = element_text(size = 20)))  # Adjust legend title size



# Graphing WOL and the viral RNA (MNE)
graph_RNA <- ggplot(data_final_graph) +
  geom_point(aes(x = WOL_MNE, y = log_MNE, col = Time_group), size = 5)+
  geom_smooth(aes(x = WOL_MNE, y = log_MNE), method = "glm") +
  scale_color_manual(
    values = c("lightblue", "blue", "pink", "orange", "red"),  # Specify the colors
    labels = c("24-48 Hours", "58-120 Hours", "130-192 Hours", "202-264 Hours", "288-408 Hours"), 
    name = "Time"# Specify the labels
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 35),
    axis.title.x = element_markdown(),
    axis.title.y = element_text(size = 35, face = "bold"),
    axis.text = element_text(size = 30),
    legend.text = element_text(size = 20),  # Adjust legend text size
    legend.title = element_text(size = 25)  # Adjust legend title size
  ) +
  ylab(expression(paste(log["10"], "DCV " , "Mean " , "Normalised ", "Expression "))) +
  xlab(expression("Wolbachia Density")) +
  guides(fill = "none")

# Graphing Wolbachia vs Time 
graph_Wol_Time <- ggplot(data_final_graph) +
  geom_point(aes(x = Time, y = WOL_MNE), size = 5)+
  geom_smooth(aes(x = Time, y = WOL_MNE), method = "glm")+
  theme_bw() +
  theme(
    axis.title = element_text(size = 35, face = "bold"),
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
    axis.text = element_text(size = 30),
    legend.text = element_text(size = 20),  # Adjust legend text size
    legend.title = element_text(size = 25)  # Adjust legend title size
  ) +
  ylab(expression("Wolbachia Density")) +
  xlab(expression("Time (hours)")) +
  guides(fill = "none") 

# Runs all the plots
graph_IU
graph_RNA
graph_Wol_Time

# ADDRESSING NON-DETECT  --------------------------------------------------
# Ιnitial approach. For log_MNE and log_IU_fly that are NA add in the minimum 
# recorded value  

# Replacing NA for non-detects with the minimum value recorded (just to get a vibe)
data_NON_DETECT <-  data_final_graph 
data_NON_DETECT$log_MNE[is.na(data_NON_DETECT$log_MNE)] <- min(data_NON_DETECT$log_MNE, na.rm = T)
data_NON_DETECT$log_IU_fly[is.na(data_NON_DETECT$log_IU_fly)] <- min(data_NON_DETECT$log_IU_fly, na.rm = T)
data_NON_DETECT <-  data_NON_DETECT %>% 
  subset( Sample != 159) # Remove sample that is missing not ND 

# Graphing WOL and the virus particles (IU) NON-DETECTS
graph_IU_NON_DETECT <- ggplot(data_NON_DETECT) +
  geom_point(aes(x = WOL_MNE, y = log_IU_fly, col = Time_group), size = 5)+
  geom_smooth(aes(x = WOL_MNE, y = log_IU_fly), method = "glm") +
  scale_color_manual(
    values = c("lightblue", "blue", "pink", "orange", "red"),  # Specify the colors
    labels = c("0-48 Hours", "58-120 Hours", "130-192 Hours", "202-264 Hours", "288-408 Hours"), 
    name = "Time"# Specify the labels
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 35, face = "bold"),
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
    axis.text = element_text(size = 30),
    legend.text = element_text(size = 20),  # Adjust legend text size
    legend.title = element_text(size = 25)
  ) +
  ylab(expression("log(Viral Particles)")) +
  xlab(expression("Wolbachia Density")) +
  guides(fill = "none") 

# Graphing WOL and the viral RNA (MNE)
graph_RNA_NON_DETECT <- ggplot(data_NON_DETECT) +
  geom_point(aes(x = WOL_MNE, y = log_MNE, col = Time_group), size = 5)+
  geom_smooth(aes(x = WOL_MNE, y = log_MNE), method = "glm") +
  scale_color_manual(
    values = c("lightblue", "blue", "pink", "orange", "red"),  # Specify the colors
    labels = c("0-48 Hours", "58-120 Hours", "130-192 Hours", "202-264 Hours", "288-408 Hours"), 
    name = "Time"# Specify the labels
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 35, face = "bold"),
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
    axis.text = element_text(size = 30),
    legend.text = element_text(size = 20),  # Adjust legend text size
    legend.title = element_text(size = 25)
  ) +
  ylab(expression("log(Viral RNA)")) +
  xlab(expression("Wolbachia Density")) +
  guides(fill = "none") 

# Run graphs 
graph_IU_NON_DETECT
graph_RNA_NON_DETECT

# STATISTICAL ANALYSIS  ---------------------------------------------------

# CORRELATION TEST --------------------------------------------------------
# Correlation between WOL and RNA 
cor_RNA <-
  cor.test(
    data_final_stat$log_MNE, 
    data_final_stat$WOL_MNE, 
    method = "spearman")

# Correlation between WOL and IU 
cor_IU <-
  cor.test(
    data_final_stat$log_IU_fly,
    data_final_stat$WOL_MNE,
    method = "spearman",
    exact = FALSE
  )

# Correlation between WOL and Time
cor_time <-
  cor.test(
    data_final_stat$Time,
    data_final_stat$WOL_MNE,
    method = "spearman",
    exact = FALSE
  )

cor_pval <- rbind(cor_IU$p.value, cor_RNA$p.value, cor_time$p.value)
cor_rho <- rbind(cor_IU$estimate, cor_RNA$estimate, cor_time$estimate)
cor_names <- rbind("Viral Particles", "Viral RNA", "Time")

cor_all <- as.data.frame(cbind(cor_rho, cor_pval, cor_names)) %>% 
  rename(pvalue = V2) %>% 
  rename(Correlation = V3)


# Definitely a correlation 29/09/23 
# Once you account for time with the GLM then Wolbachia does not have a sig impact. 
# Don't need to remove NA from data 

# GLM Wolbachia vs Time  --------------------------------------------------

# WOL vs time 
model_1_Time <- glm(formula =  WOL_MNE ~  Time,
                    data = data_final_stat)

# log10 WOL vs time 
model_2_Time <- glm(formula =  log_WOL_MNE ~  Time,
                    data = data_final_stat)

# Comparing models 
AIC(model_1_Time, model_2_Time)

# output 
summary(model_2_Time)

# GLM RNA  ----------------------------------------------------------------
# Analysis done considering time as a continuous variable 

## Un-transformed ##

# without time 

model_1_RNA <- glm(formula =  log_MNE ~ WOL_MNE,
                   data = data_final_stat)

# Without interaction 
model_2_RNA <- glm(formula =  log_MNE ~ WOL_MNE + Time,
                   data = data_final_stat)

# With interaction
model_3_RNA <- glm(formula =  log_MNE ~ WOL_MNE * Time,
                   data = data_final_stat)

# summary of models 
summary(model_1_RNA)
summary(model_2_RNA)
summary(model_3_RNA)

# Best model AIC 
AIC(model_1_RNA, model_2_RNA, model_3_RNA) # model_2_RNA

## Log transformed ## 
# Without time
model_1_RNA_log <- glm(formula =  log_MNE ~ log_WOL_MNE ,
                       data = data_final_stat)

# Without interaction 
model_2_RNA_log <- glm(formula =  log_MNE ~ log_WOL_MNE + Time,
                       data = data_final_stat)

# With interaction
model_3_RNA_log <- glm(formula =  log_MNE ~ log_WOL_MNE * Time,
                       data = data_final_stat)

# Summary 
summary(model_1_RNA_log)
summary(model_2_RNA_log)
summary(model_3_RNA_log)

# Best model AIC 
AIC(model_1_RNA_log, model_2_RNA_log, model_3_RNA_log) # model_2_RNA_log

# GLM VIRUS PARTICLE ------------------------------------------------------

## Un-transformed ##

# without time 

model_1_IU <- glm(formula =  log_IU_fly ~ WOL_MNE,
                  data = data_final_stat)

# Without interaction 
model_2_IU <- glm(formula =  log_IU_fly ~ WOL_MNE + Time,
                  data = data_final_stat)

# With interaction
model_3_IU <- glm(formula =  log_IU_fly ~ WOL_MNE * Time,
                  data = data_final_stat)

# summary of models 
summary(model_1_IU)
summary(model_2_IU)
summary(model_3_IU)

# Best model AIC 
AIC(model_1_IU, model_2_IU, model_3_IU) # model_2_IU

## Log transformed ## 
# Without time
model_1_IU_log <- glm(formula =  log_IU_fly ~ log_WOL_MNE ,
                      data = data_final_stat)

# Without interaction 
model_2_IU_log <- glm(formula =  log_IU_fly ~ log_WOL_MNE + Time,
                      data = data_final_stat)

# With interaction
model_3_IU_log <- glm(formula =  log_IU_fly ~ log_WOL_MNE * Time,
                      data = data_final_stat)

# Summary 
summary(model_1_IU_log)
summary(model_2_IU_log)
summary(model_3_IU_log)

# Best model AIC 
AIC(model_1_IU_log, model_2_IU_log, model_3_IU_log) # model_2_IU_log

# GLM ASSUMPTIONS ---------------------------------------------------------
# Plotting the assumptions 
plot(model_2_RNA)
plot(model_2_IU)
plot(model_2_RNA_log)
plot(model_2_IU_log)

# Stick with un-transformed and I don't think there is another distribution that really explains 
# apart from normal so keep it? Transforming it really skews the scale location 
# point 77 seems like an outlier 

# 1. Normal distributed residuals # FAIL (RNA)
# Create a QQ plot of residuals
ggqqplot(residuals(model_2_RNA)) # Within CI 
ggqqplot(residuals(model_2_IU)) # Within CI 
ggqqplot(residuals(model_2_RNA_log)) # Within CI 
ggqqplot(residuals(model_2_IU_log)) # Within CI 


# Create a plot of the data to see if it is normal
ggdensity(data_final_stat, x = "WOL_MNE", fill = "grey") +
  stat_overlay_normal_density(color = "red", linetype = "dashed")
skewness(data_final_stat$WOL_MNE, na.rm = TRUE) # moderately skewed

ggdensity(data_final_stat, x = "log_WOL_MNE", fill = "grey")+
  stat_overlay_normal_density(color = "red", linetype = "dashed") 
skewness(data_final_stat$log_WOL_MNE, na.rm = TRUE) # moderately skewed

# Conduct shapiro test on residuals
shapiro_test(residuals(model_2_RNA)) # FAIL
shapiro_test(residuals(model_2_IU)) # PASS 
shapiro_test(residuals(model_2_RNA_log)) # FAIL 
shapiro_test(residuals(model_2_IU_log)) # PASS 

## THEREFORE: Regardless of transformation the distribution of the residuals is 
# still not normal keep it un-transformed to make conclusions easier## 
# Failed the formal shapiro_test but visually it looks pretty good 

# 2. Homoscedesiticity of residuals - PASS
bptest(model_2_RNA, studentize = FALSE) # PASS
bptest(model_2_IU, studentize = FALSE) # PASS 


gqtest(model_2_RNA, 
       point = 0.5, 
       fraction = 0, 
       alternative = c("greater", "two.sided", "less"),
       order.by = NULL, 
       data = data_final_stat)

gqtest(model_2_IU, 
       point = 0.5, 
       fraction = 0, 
       alternative = c("greater", "two.sided", "less"),
       order.by = NULL, 
       data = data_final_stat)

# 3. Outliers - FAIL
outlierTest(model_2_RNA) # 77
outlierTest(model_2_IU)  # 77 

# See what sample the outlier is
outlier77 <- data_final_stat[77,] %>% 
  select(Sample, Time, Fly_Line, Titre, Replicate, log_IU_fly, log_MNE, log_WOL_MNE)
sample <- data_final_stat[54,]


# Filter out the outlier 
data_out <- data_final_stat[-77,]
data_out_rep <- rbind(data_out,sample) 

# Filter out a random sample to see if impact is from removing a sample or the outlier specifically
data_out_2 <- data_final_stat[-54,]


# fitting the outlier

model_2_RNA_out <- glm(formula =  log_MNE ~ WOL_MNE + Time,
                       data = data_out)
model_2_RNA_out_2<- glm(formula =  log_MNE ~ WOL_MNE + Time,
                        data = data_out_2)
model_2_RNA_out_3<- glm(formula =  log_MNE ~ WOL_MNE + Time,
                        data = data_out_rep)
summary(model_2_RNA_out_3)
summary(model_2_RNA_out_2)
summary(model_2_RNA_out)
summary(model_2_RNA)

# re-testing assumption 
# 1. Normal distributed residuals 
shapiro_test(residuals(model_2_RNA_out))
shapiro_test(residuals(model_2_RNA_out_2))

# Create a QQ plot of residuals
ggqqplot(residuals(model_2_RNA_out)) # Within CI 

# Create a plot of the data to see if it is normal
ggdensity(data_out, x = "WOL_MNE", fill = "grey") +
  stat_overlay_normal_density(color = "red", linetype = "dashed")
skewness(data_out$WOL_MNE, na.rm = TRUE) # moderately skewed

# 2. Homoscedesiticity of residuals 
bptest(model_2_RNA_out, studentize = FALSE)
gqtest(model_2_RNA_out, point = 0.5, fraction = 0, alternative = c("greater", "two.sided", "less"),order.by = NULL, data = data_grow_M)

# 3. Outliers 
outlierTest(model_2_RNA_out)

## CONCLUSION ## 
# Leave the outlier in and conduct the analysis
# 4.  Residuals independent - PASS 
bgtest(model_2_RNA, order = 1, order.by = NULL, type = c("Chisq", "F"), data = data_final_stat) # PASS pvalue nonsig
bgtest(model_2_IU, order = 1, order.by = NULL, type = c("Chisq", "F"), data = data_final_stat) # PASS pvalue nonsig


# STATS FOR NON-DETECTS ---------------------------------------------------
# Replaced the non-detectable MNE and TCID50 data with the minimum amount observed 
data_NON_DETECT 

# CORRELATION TEST --------------------------------------------------------
# Correlation between WOL and RNA - NO
cor_RNA_ND <-
  cor.test(
    data_NON_DETECT$log_MNE, 
    data_NON_DETECT$WOL_MNE, 
    method = "spearman", 
    exact = FALSE)

# Correlation between WOL and IU - NO
cor_IU_ND <-
  cor.test(
    data_NON_DETECT$log_IU_fly,
    data_NON_DETECT$WOL_MNE,
    method = "spearman",
    exact = FALSE
  )

# Correlation between WOL and Time - NO 
cor_time_ND <-
  cor.test(
    data_NON_DETECT$Time,
    data_NON_DETECT$WOL_MNE,
    method = "spearman",
    exact = FALSE
  )

cor_pval_nd <- rbind(cor_IU_ND$p.value, cor_RNA_ND$p.value, cor_time_ND$p.value)
cor_rho_nd <- rbind(cor_IU_ND$estimate, cor_RNA_ND$estimate, cor_time_ND$estimate)
cor_names_nd <- rbind("Viral Particles", "Viral RNA", "Time")

cor_all_nd <- as.data.frame(cbind(cor_rho_nd, cor_pval_nd, cor_names_nd)) %>% 
  rename(pvalue = V2) %>% 
  rename(Correlation = V3)

# No statistically sig correlations

# GLM Wolbachia vs Time  --------------------------------------------------

# WOL vs time - No relationship
model_1_Time_nd <- glm(formula =  WOL_MNE ~  Time,
                       data = data_NON_DETECT)

# log10 WOL vs time 
model_2_Time_nd <- glm(formula =  log_WOL_MNE ~  Time,
                       data = data_NON_DETECT)

# Comparing models 
AIC(model_1_Time_nd, model_2_Time_nd)

# output 
summary(model_2_Time_nd) # No impact 

# GLM RNA  ----------------------------------------------------------------
# Analysis done considering time as a continuous variable 

## Un-transformed ##

# without time 

model_1_RNA_nd <- glm(formula =  log_MNE ~ WOL_MNE,
                      data = data_NON_DETECT)

# Without interaction 
model_2_RNA_nd <- glm(formula =  log_MNE ~ WOL_MNE + Time,
                      data = data_NON_DETECT)

# With interaction
model_3_RNA_nd <- glm(formula =  log_MNE ~ WOL_MNE * Time,
                      data = data_NON_DETECT)

# summary of models 
summary(model_1_RNA_nd)
summary(model_2_RNA_nd)
summary(model_3_RNA_nd)

# Best model AIC 
AIC(model_1_RNA_nd, model_2_RNA_nd, model_3_RNA_nd) # model_2_RNA_nd

## Log transformed ## 
# Without time
model_1_RNA_log_nd <- glm(formula =  log_MNE ~ log_WOL_MNE ,
                          data = data_NON_DETECT)

# Without interaction 
model_2_RNA_log_nd <- glm(formula =  log_MNE ~ log_WOL_MNE + Time,
                          data = data_NON_DETECT)

# With interaction
model_3_RNA_log_nd <- glm(formula =  log_MNE ~ log_WOL_MNE * Time,
                          data = data_NON_DETECT)

# Summary 
summary(model_1_RNA_log_nd)
summary(model_2_RNA_log_nd)
summary(model_3_RNA_log_nd)

# Best model AIC 
AIC(model_1_RNA_log_nd, model_2_RNA_log_nd, model_3_RNA_log_nd) # model_2_RNA_log_nd

# GLM VIRUS PARTICLE ------------------------------------------------------

## Un-transformed ##

# without time 

model_1_IU_nd <- glm(formula =  log_IU_fly ~ WOL_MNE,
                     data = data_NON_DETECT)

# Without interaction 
model_2_IU_nd <- glm(formula =  log_IU_fly ~ WOL_MNE + Time,
                     data = data_NON_DETECT)

# With interaction
model_3_IU_nd <- glm(formula =  log_IU_fly ~ WOL_MNE * Time,
                     data = data_NON_DETECT)

# summary of models 
summary(model_1_IU_nd)
summary(model_2_IU_nd)
summary(model_3_IU_nd)

# Best model AIC 
AIC(model_1_IU_nd, model_2_IU_nd, model_3_IU_nd) # model_2_IU_nd

## Log transformed ## 
# Without time
model_1_IU_log_nd <- glm(formula =  log_IU_fly ~ log_WOL_MNE ,
                         data = data_NON_DETECT)

# Without interaction 
model_2_IU_log_nd <- glm(formula =  log_IU_fly ~ log_WOL_MNE + Time,
                         data = data_NON_DETECT)

# With interaction
model_3_IU_log_nd <- glm(formula =  log_IU_fly ~ log_WOL_MNE * Time,
                         data = data_NON_DETECT)

# Summary 
summary(model_1_IU_log_nd)
summary(model_2_IU_log_nd)
summary(model_3_IU_log_nd)

# Best model AIC 
AIC(model_1_IU_log_nd, model_2_IU_log_nd, model_3_IU_log_nd) # model_2_IU_log_nd

# GLM ASSUMPTIONS ---------------------------------------------------------
# Plotting the assumptions 
plot(model_2_RNA_nd)
plot(model_2_IU_nd)
plot(model_2_RNA_log_nd)
plot(model_2_IU_log_nd)

# Stick with un-transformed and I don't think there is another distribution that really explains 
# apart from normal so keep it? Transforming it really skews the scale location 


# 1. Normal distributed residuals # FAIL (RNA)
# Create a QQ plot of residuals
ggqqplot(residuals(model_2_RNA_nd)) # Some outside 
ggqqplot(residuals(model_2_IU_nd)) # Quite a lot outside
ggqqplot(residuals(model_2_RNA_log_nd)) # end samples outside
ggqqplot(residuals(model_2_IU_log_nd)) # Quite a lot outside


# Create a plot of the data to see if it is normal (not really important not part of the assumption)
ggdensity(data_NON_DETECT, x = "WOL_MNE", fill = "grey") +
  stat_overlay_normal_density(color = "red", linetype = "dashed")
skewness(data_NON_DETECT$WOL_MNE, na.rm = TRUE) # moderately skewed

ggdensity(data_NON_DETECT, x = "log_WOL_MNE", fill = "grey")+
  stat_overlay_normal_density(color = "red", linetype = "dashed") 
skewness(data_NON_DETECT$log_WOL_MNE, na.rm = TRUE) # moderately skewed

# Conduct shapiro test on residuals
shapiro_test(residuals(model_2_RNA_nd)) # FAIL
shapiro_test(residuals(model_2_IU_nd)) # FAIL 
shapiro_test(residuals(model_2_RNA_log_nd)) # FAIL 
shapiro_test(residuals(model_2_IU_log_nd)) # FAIL

## THEREFORE: Regardless of transformation the distribution of the residuals is 
# still not normal keep it un-transformed to make conclusions easier## 
# Failed the formal shapiro_test but visually it looks pretty good 

# 2. Homoscedesiticity of residuals - PASS
bptest(model_2_RNA_nd, studentize = FALSE) # PASS
bptest(model_2_IU_nd, studentize = FALSE) # PASS 


gqtest(model_2_RNA_nd, # PASS
       point = 0.5, 
       fraction = 0, 
       alternative = c("greater", "two.sided", "less"),
       order.by = NULL, 
       data = data_NON_DETECT)

gqtest(model_2_IU_nd, # PASS (if you round pvalue)
       point = 0.5, 
       fraction = 0, 
       alternative = c("greater", "two.sided", "less"),
       order.by = NULL, 
       data = data_NON_DETECT)

# 3. Outliers - FAIL (Same approach as above don't remove as there is no valid justification)
outlierTest(model_2_RNA_nd) # 77
outlierTest(model_2_IU_nd)  # 77 


# 4.  Residuals independent - PASS 
bgtest(model_2_RNA_nd, order = 1, order.by = NULL, type = c("Chisq", "F"), data = data_NON_DETECT) # PASS pvalue nonsig
bgtest(model_2_IU_nd, order = 1, order.by = NULL, type = c("Chisq", "F"), data = data_NON_DETECT) # PASS pvalue nonsig

summary(model_2_RNA)
summary(model_2_RNA_nd)

summary(model_2_IU)
summary(model_2_IU_nd)


# EXPORTING DATA ---------------------------------------

# Saving Analysed Data ----------------------------------------------------

# ANALYSED saves the MNE of analysed data

name <- paste(date, "WOL", "Data", "Analysed", sep = "_")
csv_name <- paste(name, "csv", sep = ".")
write.csv(data_final_graph, file = paste("data/analysed/",csv_name, sep =""), row.names = FALSE)

# REPEATS saves the samples that after QC need to be repeated

name <- paste(date, "WOL", "repeats", "final", sep = "_")
csv_name <- paste(name, "csv", sep = ".")
csv_name_directory <- paste("data/repeats/",csv_name, sep = "")
write.csv(repeats_full_2, file = csv_name_directory, row.names = FALSE)

# Saving plots --------------------------------------------------------------

# Saves Infectious Units vs Wolbachia
name <- paste(date, "results", "plot","glm","IU", sep = "_")
pdf_name <- paste(name, "pdf", sep = ".")
ggsave(filename = pdf_name,
       path = "results",
       plot = graph_IU, 
       width = 10, 
       height = 10)

# Saves RNA vs Wolbachia 
name <- paste(date, "results","plot","glm", "RNA", sep = "_")
pdf_name <- paste(name, "pdf", sep = ".")
ggsave(filename = pdf_name, 
       path = "results", 
       plot = graph_RNA,
       width = 10, 
       height = 10)

# Saves Time vs Wolbachia 
name <- paste(date, "results", "plot", "glm", "Time", sep = "_")
pdf_name <- paste(name, "pdf", sep = ".")
ggsave(filename = pdf_name, 
       path = "results", 
       plot = graph_Wol_Time ,
       width = 10, 
       height = 10)

## NON-DETECT plots
# Saves Infectious Units vs Wolbachia
name <- paste(date, "results", "plot", "glm", "IU",  "NON", "DETECT", sep = "_")
pdf_name <- paste(name, "pdf", sep = ".")
ggsave(filename = pdf_name, path = "results", plot = graph_IU_NON_DETECT,
       width = 10, 
       height = 10)

# Saves RNA vs Wolbachia 
name <- paste(date, "results","plot", "glm", "RNA",  "NON", "DETECT", sep = "_")
pdf_name <- paste(name, "pdf", sep = ".")
ggsave(filename = pdf_name, path = "results", plot = graph_RNA_NON_DETECT,
       width = 10, 
       height = 10)

# Saving Correlation  -----------------------------------------------------

sink("results/231219_results_output_spearman_correlation_RNA.txt") # change file name
print(cor_RNA) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**

sink("results/231219_results_output_spearman_correlation_IU.txt") # change file name
print(cor_IU) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**

sink("results/231219_results_output_spearman_correlation_Time.txt") # change file name
print(cor_time) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**

write.csv(cor_all,
          file = "/Volumes/SENV/Research/Johnson/Lab documents/11 Research Students 2022-2023/04 Results /MA 05 Wolbachia Density /R analysis/results/231219_results_output_spearman_correlation_all.csv",
          row.names = FALSE)

# Saving GLM  -------------------------------------------------------------

### TIME ### 
#GLM TIME vs WOL 
sink("results/231006_TIME_WOL_GLM_output.txt") # change file name
print(summary(model_1_Time)) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**

### RNA #### 

# GLM Wol vs TIME 
sink("results/230919_results_output_GLM_Time.txt") # change file name
print(summary(model_1_Time)) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**


# GLM RNA output without time 
sink("results/230919_results_output_GLM_RNA_no_time.txt") # change file name
print(summary(model_1_RNA)) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**

# GLM RNA with time
sink("results/231005_results_output_GLM_RNA_time_THIS.txt") # change file name
print(summary(model_2_RNA)) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**

# GLM RNA with time and interaction
sink("results/230919_results_output_GLM_RNA_time_interaction.txt") # change file name
print(summary(model_3_RNA)) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**

#### IU ####
### IU #### 
# GLM IU output without time 
sink("results/230919_results_output_GLM_IU_no_time.txt") # change file name
print(summary(model_1_IU)) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**

# GLM IU with time
sink("results/231005_results_output_GLM_IU_time_THIS.txt") # change file name
print(summary(model_2_IU)) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**

# GLM IU with time and interaction
sink("results/230919_results_output_GLM_IU_time_interaction.txt") # change file name
print(summary(model_3_IU)) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**


# NON-DETECTS GLM OUTPUT --------------------------------------------------

### NON-Detects ### 
### TIME ### 
#GLM TIME vs WOL 
sink("results/231006_results_output_GLM_Time_ND.txt") # change file name
print(summary(model_1_Time_nd)) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**

### RNA #### 
# GLM RNA output without time 
sink("results/230919_results_output_GLM_RNA_no_Time_ND.txt") # change file name
print(summary(model_1_RNA_nd)) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**

# GLM RNA with time
sink("results/231005_results_output_GLM_RNA_Time_ND.txt") # change file name
print(summary(model_2_RNA_nd)) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**

# GLM RNA with time and interaction
sink("results/230919_results_output_GLM_RNA_Time_interaction-ND.txt") # change file name
print(summary(model_3_RNA_nd)) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**

### IU #### 
# GLM IU output without time 
sink("results/230919_results_output_GLM_IU_no_Time_ND.txt") # change file name
print(summary(model_1_IU_nd)) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**

# GLM IU with time
sink("results/231005_results_output_GLM_IU_Time_ND.txt") # change file name
print(summary(model_2_IU_nd)) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**

# GLM IU with time and interaction
sink("results/230919_results_output_GLM_IU_Time_interaction_ND.txt") # change file name
print(summary(model_3_IU_nd)) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**

# old ---------------------------------------------------------------------
# FULL DATASET ------------------------------------------------------------

# To create a data frame of the finalized samples
# These have passed QC and been analysed
Check <- left_join(MNE, QC_data_2, by = join_by(Sample, Sample_file))
# Left join because I only want the QC data for those that have been analysed
data_MNE_final <- left_join(MNE, QC_data_2, by = join_by(Sample, Sample_file)) %>% 
  filter(Quality_Control == "PASS") # Gets rid of the SEM ones 

# Adding all the relevant information to those samples 
data_MNE_ALL <- left_join(Log, data_MNE_final, by = join_by(Sample)) %>% 
  select(-Sample_file_rpl, 
         -Sample_file_WOL, 
         -Gene_formated_rpl, 
         -Gene_formated_WOL, 
         -Sample_number_correct_rpl,
         -Sample_number_correct_WOL, 
         -Gene_rpl, 
         -Gene_WOL) %>% 
  mutate(IU_fly = as.numeric(IU_fly)) %>% 
  mutate(MNE = as.numeric(MNE)) %>% 
  mutate(log_IU_fly = log10(IU_fly)) %>% 
  mutate(log_MNE = log10(MNE)) %>% 
  mutate(WOL_MNE = as.numeric(WOL_MNE)) %>% 
  mutate(log_WOL_MNE = log10(WOL_MNE)) %>% 
  mutate(Time_group = case_when(
    Time %in% c(0, 6, 10, 24, 34, 48) ~ 1,
    Time %in% c(58, 72, 82, 96, 106, 120) ~ 2,
    Time %in% c(130, 144, 154, 168, 178, 192) ~ 3,
    Time %in% c(202, 216, 226, 240, 250, 264) ~ 4,
    Time %in% c(288, 312, 336, 360, 384, 408) ~ 5,
    TRUE ~ NA_integer_
  )) %>%
  arrange(Time) %>% 
  mutate(Time_group = as.character(Time_group))

# Getting rid of the duplicate rows that contain the ct for the Wol

data_MNE_final_log <- data_MNE_ALL %>% 
  select(Sample, Fly_Line, Titre, Time, Time_group, log_IU_fly, log_MNE, WOL_MNE,log_WOL_MNE) %>% 
  mutate(log_IU_fly = as.numeric(log_IU_fly)) %>% 
  mutate(WOL_MNE = as.numeric(WOL_MNE)) %>% 
  mutate(log_MNE = as.numeric(log_MNE)) %>% 
  distinct()

# Coding Time as a continous variable 
data_MNE_final_log_FOR_TIME <- data_MNE_final_log %>% 
  mutate(Time = as.numeric(Time))


# Defining 30 different unique colours 
library(RColorBrewer)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_sample <- col_vector %>% 
  sample(size = 5) %>% 
  as.data.frame() %>% 
  rename(colour = ".")

# Defining the time points 
times_sampled <- as.data.frame(unique(data_MNE_final_log$Time_group)) 

time.colors <- times_sampled %>% 
  bind_cols(col_sample) %>% 
  rename(Time = `unique(data_MNE_final_log$Time_group)`, Colour = colour) %>% 
  mutate(Time = as.numeric(Time)) %>% 
  arrange(Time)

time_colour  <- setNames(as.character(time.colors$Colour), time.colors$Time)

