############################################################################
# This script runs the development and evaluation of Methylation Risk Score.
# Using effect sizes generated from MOA and MOMENT.
# This script applies to all pairs of test and training sets, including European ancestry only subset for sensitivity analyses.

# Author: Li Ying Thong
# Jan 2025
# Software: R
############################################################################

# Load libraries
library(broom)
library(pROC)
library(dplyr)

# Read probe effect sizes and p-values generated using MOA or MOMENT from training samples
MWAS <- read.table("/file_path/MOA_Endometrium_training_CIR.moa", header = TRUE)
# Read SVA adjusted DNAm beta values of test samples
load("/file_path/myprofile_data_test_ORMadj.Robj")
# Read PRS
PRS <- read.table("/file_path/PRS_Endometrium.sscore")
# Read samplesheet
samplesheet <- read.csv("/file_path/Samplesheet_for_analysis.csv", header = TRUE)

# Divide probe effect sizes based on p-value thresholds 
## write table of DMPs with p<0.5
index <- which(MWAS$p < 0.5)
MWAS_removed_0.5 <- MWAS[index,]
head(MWAS_removed_0.5)

## write table of DMPs with p<0.2
index <- which(MWAS$p < 0.2)
MWAS_removed_0.2 <- MWAS[index,]
head(MWAS_removed_0.2)

## write table of DMPs with p<0.1
index <- which(MWAS$p < 0.1)
MWAS_removed_0.1 <- MWAS[index,]
head(MWAS_removed_0.1)

## write table of DMPs with p<1e-2
index <- which(MWAS$p < 1e-2)
MWAS_removed_e2 <- MWAS[index,]
head(MWAS_removed_e2)

## write table of DMPs with p<1e-3
index <- which(MWAS$p < 1e-3)
MWAS_removed_e3 <- MWAS[index,]
head(MWAS_removed_e3)

## write table of DMPs with p<1e-4
index <- which(MWAS$p < 1e-4)
MWAS_removed_e4 <- MWAS[index,]
head(MWAS_removed_e4)

## write table of DMPs with p<1e-5
index <- which(MWAS$p < 1e-5)
MWAS_removed_e5 <- MWAS[index,]
head(MWAS_removed_e5)

# Get Endometriosis case-control status
phenotype_cov_data <- samplesheet
unique(phenotype_cov_data$Endometriosis..Yes.No.)
index <- which(phenotype_cov_data$Endometriosis..Yes.No. == "Yes")
phenotype_cov_data$Endometriosis..Yes.No.[index] <- "1"
index <- which(phenotype_cov_data$Endometriosis..Yes.No. == "No")
phenotype_cov_data$Endometriosis..Yes.No.[index] <- "0"
phenotype_cov_data$Endometriosis..Yes.No. <- as.factor(phenotype_cov_data$Endometriosis..Yes.No.)

# Format PRS data frame for AUC calculation
PRS <- PRS[, c(1, 6)]
phenotype_cov_data <- merge(phenotype_cov_data, PRS, by.x = "Epic_Complete.Bar.code", by.y = "V1",)
colnames(phenotype_cov_data)[which(colnames(phenotype_cov_data) == "V6")] <- "PRS"

# Check beta values of test samples
index <- which(samplesheet$Institute.for.Analysis == "CIR") # Select test samples according to institute 
training <- samplesheet[-index,] 
test <- samplesheet[index,] 
myprofile_data <- myprofile_data_test
myprofile_data <- as.data.frame(myprofile_data, fix.empty.names = TRUE)
myprofile_data_test_IID <- myprofile_data[myprofile_data$IID %in% test$Epic_Complete.Bar.code,]

# Calculate MRS using p-value threshold < 1e-5
weights_e5 <- MWAS_removed_e5[, c(2, 6)]
index <- which(colnames(myprofile_data_test_IID) %in% weights_e5$Probe)
myprofile_data_test_IID_e5 <- myprofile_data_test_IID[, index]
colnames_seq <- colnames(myprofile_data_test_IID_e5)
weights_e5 <- weights_e5[match(colnames_seq, weights_e5$Probe),]
myprofile_data_test_IID_e5 <- as.matrix(myprofile_data_test_IID_e5)
myprofile_data_test_IID_e5_num <- matrix(as.numeric(myprofile_data_test_IID_e5), ncol = ncol(myprofile_data_test_IID_e5))
weights_e5_matrix <- weights_e5[, 2]
weights_e5_matrix <- as.matrix(weights_e5_matrix)
weights_e5_matrix_num <- matrix(as.numeric(weights_e5_matrix), ncol = ncol(weights_e5_matrix))
MRS_e5_test <- myprofile_data_test_IID_e5_num %*% weights_e5_matrix_num
MRS_e5_test_data <- as.data.frame(MRS_e5_test, fix.empty.names = TRUE)
colnames(MRS_e5_test_data) <- "MRS_e5_test"
MRS_e5_test_data$IID <- myprofile_data_test_IID$IID
MRS_e5_test_data <- MRS_e5_test_data[, c(2, 1)]

# Evaluate the performance of MRS using AUC
AUC_test_data <- merge(MRS_e5_test_data, phenotype_cov_data, by.x = "IID", by.y = "Epic_Complete.Bar.code")
colnames(AUC_test_data)[2] <- "MRS"
## logistic regression
model_test <- glm(Endometriosis..Yes.No.~MRS, family="binomial", data=AUC_test_data)
summary_model <- summary(model_test)
summary_model$coefficients[2, 1] # MRS estimate
summary_model$coefficients[2, 2] # MRS estimate SE 
summary_model$coefficients[2, 4] # MRS pvalue 
summary_model$coefficients[1, 1] # intercept estimate 
summary_model$coefficients[1, 2] # intercept estimate SE 
summary_model$coefficients[1, 4] # intercept pvalue 
## AUC 
predicted_test <- predict(model_test, AUC_test_data, type="response")
roc_object_test <- roc(AUC_test_data$Endometriosis..Yes.No., predicted_test)
AUC_test <- auc(roc_object_test) # AUC
a <- ci.auc(roc_object_test)
AUC_summary <- summary(a)
AUC_summary[6] # upper 95% CI limit of the AUC
AUC_summary[1] # lower 95% CI limit of the AUC

# Evaluate the performance of MRS and PRS combined using AUC
AUC_test_data <- merge(MRS_e5_test_data, phenotype_cov_data, by.x = "IID", by.y = "Epic_Complete.Bar.code")
colnames(AUC_test_data)[2] <- "MRS"
## logistic regression
model_test <- glm(Endometriosis..Yes.No.~MRS+PRS, family="binomial", data=AUC_test_data)
summary_model <- summary(model_test)
summary_model$coefficients[2, 1] # MRS estimate 
summary_model$coefficients[2, 2] # MRS estimate SE 
summary_model$coefficients[3, 1] # PRS estimate 
summary_model$coefficients[3, 2] # PRS estimate SE 
summary_model$coefficients[2, 4] # MRS pvalue 
summary_model$coefficients[3, 4] # PRS pvalue
## AUC
predicted_test <- predict(model_test, AUC_test_data, type="response")
roc_object_test <- roc(AUC_test_data$Endometriosis..Yes.No., predicted_test)
AUC_test <- auc(roc_object_test) # AUC
a <- ci.auc(roc_object_test)
AUC_summary <- summary(a)
AUC_summary[6] # upper 95% CI limit of the AUC
AUC_summary[1] # lower 95% CI limit of the AUC

# Repeat "Calculate MRS using p-value threshold < 1e-5", "Evaluate the performance of MRS using AUC" and "Evaluate the performance of MRS using AUC" on p-value threshold < 1e-4, 1e-3, 1e-2, 0.1, 0.2, 0.5

