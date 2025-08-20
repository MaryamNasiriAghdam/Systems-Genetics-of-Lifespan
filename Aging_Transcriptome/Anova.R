### Code for Pre-SAS Table Preparation, ANOVA in SAS, and FDR Correction (Benjamini–Hochberg)

###Pre SAS table preparation

Data_Anova <- data.table(fread("TMM_Transposed.csv", header = TRUE))


Data_Anova[substr(Data_Anova$Samples,1,2)=="D5","Age"]<-"Day5"
Data_Anova[substr(Data_Anova$Samples,1,2)=="W1","Age"]<-"Day7"
Data_Anova[substr(Data_Anova$Samples,1,2)=="W2","Age"]<-"Day14"
Data_Anova[substr(Data_Anova$Samples,1,2)=="W3","Age"]<-"Day21"
Data_Anova[substr(Data_Anova$Samples,1,2)=="W4","Age"]<-"Day28"
Data_Anova[substr(Data_Anova$Samples,1,2)=="W5","Age"]<-"Day35"
Data_Anova[substr(Data_Anova$Samples,1,2)=="W6","Age"]<-"Day42"
Data_Anova[substr(Data_Anova$Samples,1,2)=="W7","Age"]<-"Day49"
Data_Anova[substr(Data_Anova$Samples,1,2)=="W8","Age"]<-"Day56"
Data_Anova[substr(Data_Anova$Samples,1,2)=="W9","Age"]<-"Day63"
Data_Anova[substr(Data_Anova$Samples,1,3)=="W10","Age"]<-"Day70"

# relocate the Age column
Data_Anova <- Data_Anova %>% relocate(Age, .after = Samples)

# create Sex column
library("stringr")   
Data_Anova[str_sub(Data_Anova$Samples,-3,-3)=="F","Sex"]<-"Female"
Data_Anova[str_sub(Data_Anova$Samples,-3,-3)=="M","Sex"]<-"Male"

# relocate the Sex column
Data_Anova <- Data_Anova %>% relocate(Sex, .after = Age)

# create Replicate column
Data_Anova[str_sub(Data_Anova$Samples,-1,-1)=="1","Replicate"]<-"1"
Data_Anova[str_sub(Data_Anova$Samples,-1,-1)=="3","Replicate"]<-"2"

# relocate the Replicate column
Data_Anova <- Data_Anova %>% relocate(Replicate, .after = Age)

# create Tissue column
Data_Anova[str_sub(Data_Anova$Samples,-5,-5)=="B","Tissue"]<-"Body"
Data_Anova[str_sub(Data_Anova$Samples,-5,-5)=="H","Tissue"]<-"Head"
Data_Anova[str_sub(Data_Anova$Samples,-5,-5)=="O","Tissue"]<-"Reproductive"
Data_Anova[str_sub(Data_Anova$Samples,-5,-5)=="T","Tissue"]<-"Reproductive"


# relocate the Tissue column
Data_Anova <- Data_Anova %>% relocate(Tissue, .after = Replicate)



# Long table

Data_Anova_long <- Data_Anova %>% 
  pivot_longer(cols = starts_with("FBgn"),
               names_to = "FBgn", values_to = "Expression")

Data_Anova_long_df <- as.data.frame(Data_Anova_long)

write.table(Data_Anova_long_df, file = "TMM_PRE_SAS.txt", sep = "\t", row.names = TRUE)
write.csv(Data_Anova_long_df, file = "TMM_PRE_SAS.csv", row.names = TRUE)


### Code performed in SAS
/*Omnibus ANOVA by proc glm*/
  
  LIBNAME Anova '/home/maryamn/Anova';

PROC IMPORT DATAFILE="/home/maryamn/TMM_PRE_SAS.csv" DBMS=CSV
OUT=Anova.Omni_glm REPLACE;
GETNAMES=YES;
RUN;

DATA Anova.Omni_glm;
SET Anova.Omni_glm (DROP=VAR1);
RUN;

PROC CONTENTS DATA=Anova.Omni_glm; 
RUN;

proc sort data=Anova.Omni_glm out=Anova.glm_sorted;
by FBgn;
RUN;

/* format output tables*/
  PROC template;
edit Stat.GLM.ProbF; format=E12.; end;
run;
PROC template;
edit Stat.GLM.LSMSlice; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Contrasts; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Tests; format=E12.; end;
run;
proc template;
edit Stat.Mixed.Tests3; parent = Stat.Mixed.FTests; end;
run;


/*no html output or note writing to decrease run time*/
  ODS RESULTS OFF;
options nonotes;

/*Naming output datasets*/
  ods noproctitle;
ods graphics / imagemap=off;
ods output FitStatistics=FS_omni;
ods output LSMeans=LSMeans_omni;
ods output ClassLevels= ClassLevels_omni;
ods output NObs= NObs_omni;
ods output OverallANOVA= OverallANOVA_omni;
ods output ModelANOVA= ModelANOVA_omni;
ods output Contrasts=Contrasts_omni;


/*Defining and running the model for each gene (geneid)*/
  proc glm data=Anova.glm_sorted plots=none;
class Age Sex Tissue;
model Expression = Age|Sex|Tissue;
by FBgn;
lsmeans Age Age*Sex*Tissue / ; 
run;


### Correct for FDR (BH)
library(data.table)
library(tidyr)
library(dplyr)

library(bit64)

Data_Anova <- data.table(fread("Omni_Type3.csv", header=TRUE)) #hypothesis type 3

# Pivot Wider
Data_Mixed_Wider <- Data_Anova %>% 
  pivot_wider(names_from = Source, values_from = c(Source, DF, SS, MS, FValue, ProbF))

Data_Mixed_Wider_df <- as.data.frame(Data_Mixed_Wider)

write.table(Data_Mixed_Wider_df, file = "Omni_GLM_Maryam/Glm_type3_Maryam_Wider.txt", sep = "\t", row.names = TRUE)

Anova_pvalue <- read.table("Glm_type3_Maryam_Wider.txt", header = TRUE, sep = "\t")

library(rstatix)

Anova_Padj1 <- adjust_pvalue(Anova_pvalue, p.col = "ProbF_Age", output.col = "padj_Age", method = "BH")

Anova_Padj2 <- adjust_pvalue(Anova_Padj1, p.col = "ProbF_Sex", output.col = "padj_Sex", method = "BH")

Anova_Padj3 <- adjust_pvalue(Anova_Padj2, p.col = "ProbF_Tissue", output.col = "padj_Tissue", method = "BH")

Anova_Padj4 <- adjust_pvalue(Anova_Padj3, p.col = "ProbF_Age.Sex", output.col = "padj_Age*Sex", method = "BH")

Anova_Padj5 <- adjust_pvalue(Anova_Padj4, p.col = "ProbF_Age.Tissue", output.col = "padj_Age*Tissue", method = "BH")

Anova_Padj6 <- adjust_pvalue(Anova_Padj5, p.col = "ProbF_Sex.Tissue", output.col = "padj_Sex*Tissue", method = "BH")

Anova_Padj7 <- adjust_pvalue(Anova_Padj6, p.col = "ProbF_Age.Sex.Tissue", output.col = "padj_Age*Sex*Tissue", method = "BH")



write.table(Anova_Padj7, file = "GLM_Type3_Padj.txt", sep = "\t", row.names = TRUE)



