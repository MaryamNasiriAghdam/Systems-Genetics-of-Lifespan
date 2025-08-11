# Code for pre-processing Count Table and RNA-seq Normalization

raw_counts <- read.delim("combined_counts.txt", header = TRUE, sep = "\t", row.names = "Geneid")

## Calculate sums
raw_counts$X5D_O_F_3 <-apply(raw_counts[,11:12],1,sum)
raw_counts$W10_B_F_3 <-apply(raw_counts[,16:18],1,sum) # the one with three columns
raw_counts$W4_B_M_1 <-apply(raw_counts[,67:68],1,sum)
raw_counts$W4_H_M_1 <-apply(raw_counts[,72:73],1,sum)
raw_counts$W5_H_M_1 <-apply(raw_counts[,85:86],1,sum)
raw_counts$W6_B_M_1 <-apply(raw_counts[,94:95],1,sum)
raw_counts$W7_B_M_1 <-apply(raw_counts[,107:108],1,sum)
raw_counts$W7_H_F_1 <-apply(raw_counts[,110:111],1,sum)
raw_counts$W7_H_M_1 <-apply(raw_counts[,113:114],1,sum)
raw_counts$W8_H_F_1 <-apply(raw_counts[,124:125],1,sum)
raw_counts$W8_H_M_1 <-apply(raw_counts[,127:128],1,sum)
raw_counts$W9_B_M_1 <-apply(raw_counts[,136:137],1,sum)
raw_counts$W9_H_F_3 <-apply(raw_counts[,140:141],1,sum)

## Remove that columns after calculating sums
raw_counts <- raw_counts[,!names(raw_counts) %in% c("X5D_O_F_3_AGTCAA_L004_R1_C9D5KANXX.bam", 
                                                    "X5D_O_F_3_AGTCAA_L004_R1_C9EMEANXX.bam", 
                                                    "W10_B_F_3_TAGCTT_L005_R1_C9EMEANXX.bam", 
                                                    "W10_B_F_3_TAGCTT_L005_R1_C9MM3ANXX.bam", 
                                                    "W10_B_F_3_TAGCTT_L007_R1_C9D5KANXX.bam", 
                                                    "W4_B_M_1_ACTTGA_L005_R1_C9MM3ANXX.bam", 
                                                    "W4_B_M_1_ACTTGA_L007_R1_C9EMEANXX.bam", 
                                                    "W4_H_M_1_CCGTCC_L005_R1_C9MM3ANXX.bam", 
                                                    "W4_H_M_1_CCGTCC_L006_R1_C9EMEANXX.bam", 
                                                    "W5_H_M_1_GTGAAA_L005_R1_C9MM3ANXX.bam", 
                                                    "W5_H_M_1_GTGAAA_L006_R1_C9EMEANXX.bam", 
                                                    "W6_B_M_1_TGACCA_L005_R1_C9MM3ANXX.bam", 
                                                    "W6_B_M_1_TGACCA_L007_R1_C9EMEANXX.bam", 
                                                    "W7_B_M_1_CGATGT_L005_R1_C9MM3ANXX.bam", 
                                                    "W7_B_M_1_CGATGT_L007_R1_C9EMEANXX.bam", 
                                                    "W7_H_F_1_TCGGCA_L005_R1_C9MM3ANXX.bam", 
                                                    "W7_H_F_1_TCGGCA_L006_R1_C9EMEANXX.bam",
                                                    "W7_H_M_1_GAGTGG_L005_R1_C9MM3ANXX.bam", 
                                                    "W7_H_M_1_GAGTGG_L006_R1_C9EMEANXX.bam", 
                                                    "W8_H_F_1_CGGAAT_L005_R1_C9MM3ANXX.bam", 
                                                    "W8_H_F_1_CGGAAT_L007_R1_C9EMEANXX.bam", 
                                                    "W8_H_M_1_GGTAGC_L005_R1_C9MM3ANXX.bam", 
                                                    "W8_H_M_1_GGTAGC_L008_R1_C9EMEANXX.bam", 
                                                    "W9_B_M_1_CATGGC_L005_R1_C9MM3ANXX.bam", 
                                                    "W9_B_M_1_CATGGC_L008_R1_C9EMEANXX.bam",
                                                    "W9_H_F_3_CAGATC_L005_R1_C9MM3ANXX.bam",
                                                    "W9_H_F_3_CAGATC_L007_R1_C9D5KANXX.bam")]
## Clean column names
# Remove everything after -1 from Column name
names(raw_counts) = gsub("(_1).*", "\\1", x = names(raw_counts))

# Remove everything after -3 from Column name
names(raw_counts) = gsub("(_3).*", "\\1", x = names(raw_counts))

# Relocate Columns

library('dplyr')
raw_counts <- raw_counts %>% relocate(X5D_O_F_3, .after = X5D_O_F_1)
raw_counts <- raw_counts %>% relocate(W10_B_F_3, .after = W10_B_F_1)
raw_counts <- raw_counts %>% relocate(W4_B_M_1, .after = W4_B_F_3)
raw_counts <- raw_counts %>% relocate(W4_H_M_1, .after = W4_H_F_3)
raw_counts <- raw_counts %>% relocate(W5_H_M_1, .after = W5_H_F_3)
raw_counts <- raw_counts %>% relocate(W6_B_M_1, .after = W6_B_F_3)
raw_counts <- raw_counts %>% relocate(W7_B_M_1, .after = W7_B_F_3)
raw_counts <- raw_counts %>% relocate(W7_H_F_1, .after = W7_B_M_3)
raw_counts <- raw_counts %>% relocate(W7_H_M_1, .after = W7_H_F_3)
raw_counts <- raw_counts %>% relocate(W8_H_F_1, .after = W8_B_M_3)
raw_counts <- raw_counts %>% relocate(W8_H_M_1, .after = W8_H_F_3)
raw_counts <- raw_counts %>% relocate(W9_B_M_1, .after = W9_B_F_3)
raw_counts <- raw_counts %>% relocate(W9_H_F_3, .after = W9_H_F_1)

raw_counts <- raw_counts %>% relocate(W10_B_F_1, .after = last_col())
raw_counts <- raw_counts %>% relocate(W10_B_F_3, .after = last_col())
raw_counts <- raw_counts %>% relocate(W10_B_M_1, .after = last_col())
raw_counts <- raw_counts %>% relocate(W10_B_M_3, .after = last_col())
raw_counts <- raw_counts %>% relocate(W10_H_F_1, .after = last_col())
raw_counts <- raw_counts %>% relocate(W10_H_F_3, .after = last_col())
raw_counts <- raw_counts %>% relocate(W10_H_M_1, .after = last_col())
raw_counts <- raw_counts %>% relocate(W10_H_M_3, .after = last_col())
raw_counts <- raw_counts %>% relocate(W10_O_F_1, .after = last_col())
raw_counts <- raw_counts %>% relocate(W10_O_F_3, .after = last_col())
raw_counts <- raw_counts %>% relocate(W10_T_M_1, .after = last_col())
raw_counts <- raw_counts %>% relocate(W10_T_M_3, .after = last_col())

# Rename X5D to D5
raw_counts <- raw_counts %>% rename(D5_B_F_1 = X5D_B_F_1)
raw_counts <- raw_counts %>% rename(D5_B_F_3 = X5D_B_F_3)
raw_counts <- raw_counts %>% rename(D5_B_M_1 = X5D_B_M_1)
raw_counts <- raw_counts %>% rename(D5_B_M_3 = X5D_B_M_3)
raw_counts <- raw_counts %>% rename(D5_H_F_1 = X5D_H_F_1)
raw_counts <- raw_counts %>% rename(D5_H_F_3 = X5D_H_F_3)
raw_counts <- raw_counts %>% rename(D5_H_M_1 = X5D_H_M_1)
raw_counts <- raw_counts %>% rename(D5_H_M_3 = X5D_H_M_3)
raw_counts <- raw_counts %>% rename(D5_O_F_1 = X5D_O_F_1)
raw_counts <- raw_counts %>% rename(D5_O_F_3 = X5D_O_F_3)
raw_counts <- raw_counts %>% rename(D5_T_M_1 = X5D_T_M_1)
raw_counts <- raw_counts %>% rename(D5_T_M_3 = X5D_T_M_3)


## write clean raw count table
write.table(raw_counts, file = "New_Combined_Cleaned_Raw_Counts.txt", sep = "\t", row.names = TRUE)


###RNA_seq Normalization

library(dplyr)
library(matrixStats)
library(data.table)
library(tidyr)
library(edgeR)
library(limma)


combined_counts<-read.table(file="New_Combined_Cleaned_Raw_Counts.txt", header = TRUE)

combined_counts <- as.data.frame(combined_counts)

Data<-combined_counts

#Remove the length column from Data and call it Data_nolength
Data_nolength<-Data[,-1]



#Data filtered by Median count>2
Datamatrix<-data.matrix(Data_nolength)
Medians<-rowMedians(Datamatrix)
Medians<-data.frame(Medians)
row.names(Medians)<-row.names(Data_nolength)
DatacbindMedian<-cbind(Data_nolength,Medians)
filterMedian<-filter(DatacbindMedian,Medians>=2)
filterMedian<-subset(filterMedian,select=-c(Medians))

# Print the number of rows
print(nrow(filterMedian))
head(filterMedian)

write.table(filterMedian, file = "filterMedian.txt", sep = "\t", row.names = TRUE)


#Filter out genes with proportion of non zero samples greater than 0.25

# Calculate the proportion of non-zero values for each gene
non_zero_counts <- rowSums(filterMedian != 0)
non_zero_proportion <- non_zero_counts / ncol(filterMedian)

# Filter out genes where the proportion of non-zero values is greater than 0.25
filtered_genes <- filterMedian[non_zero_proportion <= 0.25, ]

# Print the number of rows in the filtered data frame
print(nrow(filtered_genes))



##Double check
# Check the distribution of non-zero proportions
hist(non_zero_proportion, breaks = 20, main = "Distribution of Non-zero Proportions", xlab = "Proportion")


## TMM normalization of log transformation

d <- DGEList(counts=filterMedian)
dim(d)
apply(d$counts, 2, sum)

d$samples$lib.size <- colSums(d$counts)
d$samples
d <- calcNormFactors(d, method = "TMM")
TMM <- cpm(d, log=T)
write.table(TMM, file = "New_TMM_Normalized_Logtransformed.txt", sep = "\t", row.names = TRUE)
TMM_1 <- cpm(d)
write.table(TMM_1, file = "New_TMM_Normalized.txt", sep = "\t", row.names = TRUE)
write.csv(TMM_1, file = "New_TMM_Normalized.csv", row.names = TRUE)

## Expression box plots (Repeat for all 6 groups)


All_TMM <- read.table("New_TMM_Normalized.txt", sep = "\t", header = TRUE, row.names = 1)
print(names(All_TMM))

# Select columns with names containing "H_M"
selected_columns <- All_TMM[, grepl("H_M", names(All_TMM))]
names <- row.names(selected_columns)

replicate_groups <- gsub("_\\d$", "", names(selected_columns))
replicate_groups


selected_columns_log <- log (selected_columns+1)


boxplot(selected_columns_log, 
        main = "Gene Expression Boxplot (Male Head)", 
        ylab = "Expression Value",
        las = 2,
        cex.axis = 0.8,
        col = replicate_colors[as.numeric(factor(replicate_groups))])
