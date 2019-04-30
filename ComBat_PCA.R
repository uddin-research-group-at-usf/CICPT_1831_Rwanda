####################################################################################
# CICPT_1831 (Rawanda) Background Corrected & QC'd: ComBat
####################################################################################


library(tibble)
library(data.table)
library(sva)
library(impute)


## Scientific to decimal
options(scipen = 999)


## convert columns to factor
make_factor <- function(fact_cols, data){
  data[fact_cols] <- lapply(data[,fact_cols], factor)
  return(data)
}


########################################################################
# Running PCA Analysis

run_pca <- function(combat_data, file_name, pca_num  ){
  
  #transpose the data
  betas <- t(combat_data)
  
  # get PCs
  PCobj <- prcomp(betas, retx = T, center = T, scale. = T)
  PCs <- PCobj$x
  
  propvar <- summary(PCobj)$importance["Proportion of Variance", 1:pca_num]
  cummvar <- summary(PCobj)$importance["Cumulative Proportion", 1:pca_num]
  
  save(PCs, file = paste0("CICPT_1831_bkgd_beta_Imputed_", file_name,
                          "_PCs_Removed Cross Hybrid_CpGOnly.Rdata"))
  
  
  # Generate Variance Plots
  
  # Plot of the proportion of variability explained by the top  PCs
  
  pdf(paste0("CICPT_1831_Imputed_bkgd_beta_", file_name,
             "_propvar_Removed Cross Hybrid_CpGOnly.pdf"))
  par(mar = c(5,5,4,2))
  barplot(propvar*100, xlab = paste("Top", 5, "PCs", sep = " "), ylab = "Variation Explained (%)",
          cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5, ylim=c(0,50))
  dev.off()
  
  
  # Plot of the cummulative proportion of variability explained by the top R PCs
  pdf(paste0("CICPT_1831_Imputed_bkgd_beta_", file_name,
             "_cummvar_Removed Cross Hybrid_CpGOnly.pdf"))
  par(mar = c(5,5,4,2))
  barplot(cummvar*100, xlab = paste("Top", pca_num, "PCs", sep = " "), 
          ylab = "Cummulative Variation Explained (%)",
          cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5, ylim=c(0,75))
  abline(a = 100*cummvar[pca_num], b = 0, lwd = 3, col = "red", lty = "dashed")
  dev.off()
  
  
  # Stop if any missing value
  stopifnot(!sum(is.na(match(rownames(samps_with_sex), rownames(PCs)))))
  
  PCs <- PCs[rownames(samps_with_sex), ]
  
  stopifnot(all(rownames(samps_with_sex) == rownames(PCs)))
  
  return(PCs)
}


## Fit an analysis for Variance model
check_variance <- function(pc.list, data, column_name){
  pvals <- NULL
  for(i in seq_along(pc.list)) {
    fit <- aov(data[, pc.list[i] ] ~ data[, column_name])
    p.val <- summary(fit)[[1]]["Pr(>F)"][[1]][1]
    pvals <- append(pvals, p.val)
  }
  return(pvals)
}



# Calculat Pvalues
calc_pvals <- function(pc.list, sampsPC){
  
  chip.pvals <- check_variance(pc.list, sampsPC, "Sentrix_ID")
  plate.pvals <- check_variance(pc.list, sampsPC, "Sample_Plate" )
  pos.pvals <- check_variance(pc.list, sampsPC, "Sentrix_Position")
  sex.pvals <- check_variance(pc.list, sampsPC, "Gender")
  pvals <- rbind(chip.pvals, plate.pvals, pos.pvals, sex.pvals)
  
  return(pvals)
}



## plot heatmap 
plot_heatmap <- function(pvals, file_name){
  rownames(pvals) <- c("Chip", "Plate", "Position", "Sex") # "Rundate"
  colnames(pvals) <- pc.list
  
  pvals.m <- melt(pvals)
  pvals.m$Var1 <- as.character(pvals.m$Var1)
  pvals.m$Var2 <- as.character(pvals.m$Var2)
  
  ## create a new column and assing some dummy value
  pvals.m$level <- "p > 0.05"
  
  ## assign actual values
  pvals.m[pvals.m[, "value"] < 0.05 & pvals.m[, "value"] > 0.01, "level"] <-"p < 0.05"
  pvals.m[pvals.m[, "value"] < 0.01 & pvals.m[, "value"] > 0.001, "level"] <-"p < 0.01"
  pvals.m[pvals.m[, "value"] < 0.001 & pvals.m[, "value"] > 0.00001, "level"] <-"p < 0.001"
  pvals.m[pvals.m[, "value"] < 0.00001, "level"] <- "p < 0.00001"
  pvals.m$level <- factor(pvals.m$level, levels = c("p > 0.05", "p < 0.05", "p < 0.01", "p < 0.001", "p < 0.00001"))
  
  
  ## assign color
  myColors <- c("white", "pink", "orange", "red", "darkred")
  names(myColors) <- levels(pvals.m$level)
  colScale <- scale_fill_manual(name = "level", values = myColors)
  
  #Make PDF with the heatmap (variables and 6 PCs)
  plot_name = paste0("CICPT_1831_bkgd_beta_", file_name,  
                     "_PCA_HEATMAP_Removed Cross Hybrid_CpGOnly.pdf")
  
  ggplot(pvals.m, aes(x = Var2, y = Var1)) + geom_tile(aes(fill = factor(level))) +
    colScale +
    xlab("Principal Components") + ylab("Technical Artifacts") +
    ggtitle(paste0("Principal Components, CICPT_1831_", file_name, "_Beta Matrix")) +
    theme(panel.border = element_blank()) + theme_bw()
  ggsave(plot_name)
  #dev.off()
  
  
  #Save the pvalues of each variable and 6 PCs to a spreadsheet
  write.csv(pvals, paste0("CICPT_1831_bkgd_beta_", file_name,  
                          "_PCA_HEATMAP_pvals_Removed Cross Hybrid_CpGOnly.csv"))
  
  
  #Create a PNG file with the heatmap
  # png("CICPT_1831_bkgd_beta_RAQ_PCA_HEATMAP_Removed Cross Hybrid_CpGOnly.png",
  #     width = 800, height = 480, bg = "transparent")
  # ggplot(pvals.m, aes(x = Var2, y = Var1)) + geom_tile(aes(fill = factor(level))) +
  #   colScale +
  #   xlab("Principal Components") + ylab("Technical Artifacts") +
  #   theme(panel.border = element_blank()) + theme_bw() +
  #   theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24),
  #         legend.text = element_text(size = 20), legend.title = element_text(size = 24))
  # dev.off()
}



# Heat Map
# 1. Sentrix_ID = chip
# 2. Sample_Plate = plate
# 3. Sentrix_position = where on the chip a sample is across all chips
# 4. Sex
# 5. Rundate (Error on Machine resulted in several samples being run on different days)

# ------------------------------------------------------------------ #
post_combat_BE <- function(combat_data, data_wd_sex, col_name, pcs){
  
  ## Post combat
  PCs <- run_pca(combat_data, col_name, pcs)
  
  ## merge data and generates PCs
  sampsPC <- merge(data_wd_sex, PCs, by = "row.names")
  
  rownames(sampsPC) <- sampsPC[, "Row.names"]
  
  all(rownames(sampsPC) == sampsPC[, "Row.names"])
  
  
  # Need all the categorical variables to be factors
  
  # After Combat
  pvals <- calc_pvals(pc.list, sampsPC)
  
  #post combat
  plot_heatmap(pvals, col_name)
}

# ------------------------------------------------------------------- #


####################################################################################
# Step 1: Load and Prep Data
####################################################################################

## Load sample sheet
samps <- read.csv("UddinMethylationSampleSheet2019.csv", header = TRUE, skip = 7)


## Combine Sentrix_ID and Position to make rownames
rownames(samps) <- paste0(samps$Sentrix_ID, "_", samps$Sentrix_Position)


#Load in the BMIQ normalized data
bmiq <- local(get(load("CICPT_1831_beta_BMIQ_FULL_bySample_RemovedCrossHybrid.Rdata")))

## Check for missing values
sum(is.na(bmiq)) #150623

# Creates flag for missing observations
missing <- which(is.na(bmiq)) 
length(missing) #150623


# Note: if there are any beta values that are 0 or 1, they will have to be modified
# before log transforming in order to prevent "-Inf"

min(bmiq, na.rm = TRUE) # no zero 
max(bmiq, na.rm = TRUE) # there is at least 1 beta value of 1
bmiq[which(bmiq == 1)] <- 0.9999999 #Changes the value of 1 to 0.9999999

bmiq <- log2(bmiq / (1-bmiq)) # log transforming
min(bmiq, na.rm = TRUE) #-73.68111
max(bmiq, na.rm = TRUE) #53

sum(bmiq == "-Inf", na.rm = TRUE) # should be zero
range(bmiq, na.rm = TRUE)

all(rownames(samps) == colnames(bmiq)) # Should be True

#IF FALSE
sampe <- samps[colnames(bmiq), ]

#check, view the files
all(rownames(sampe) == colnames(bmiq)) # should be TRUE now
samps <- sampe
all(rownames(samps) == colnames(bmiq))#Should be True
 
sum(is.na(match(rownames(samps), colnames(bmiq)))) #0
samps <- samps[colnames(bmiq),]
all(rownames(samps)==colnames(bmiq)) #TRUE

head(samps)[,1:4]
head(bmiq)[,1:4]



####################################################################################
# Step 2: Impute Data
####################################################################################

beta.imputed <- impute.knn(as.matrix(bmiq))
beta.imputed <- beta.imputed$data
save(beta.imputed, file = "CICPT_1831_betaImputed_LogTransformed_preComBat_Removed_Cross_Hybrid_CpGOnly.RData")
head(samps)[, 1:4]
head(beta.imputed)[, 1:4]

all(rownames(samps) == colnames(beta.imputed)) # should be TRUE
str(samps)

# cols to be in factor
fact_cols <- c("Sample_Name", "Sample_Well", "Sample_Plate", "Sample_Group",
               "Sentrix_ID", "Sentrix_Position")

samps <- make_factor(fact_cols, samps)


head(samps)[, 1:4]
head(beta.imputed)[, 1:4]


## load data with gender information
gender_data <- local(get(load("CICPT_1831_data_with_gender.Rdata")))
gender_data <- gender_data[, c("Gender", "sentid_sentpos")]

## merge by rownames
gender_data <- merge(samps, gender_data, by.x = "row.names", by.y = "sentid_sentpos")

## rownames to column
samps_with_sex <- column_to_rownames(gender_data, var = "Row.names")


#add exposure status
non_expo <- which(grepl("CT", samps_with_sex$Sample_Name))
expo <- which(!grepl("CT", samps_with_sex$Sample_Name))


## creat dummy and and yes/no exposure status
samps_with_sex$Exposure <- "dummy"
samps_with_sex$Exposure[non_expo] <- "No"
samps_with_sex$Exposure[expo] <- "Yes"



## PCs to plot
pc.list <- c("PC1", "PC2", "PC3", "PC4", "PC5")


## Before running Combat
samps_with_sex <- make_factor(c(fact_cols, "Gender"), samps_with_sex)
preCombat_PCs <- run_pca(beta.imputed, "PreCombat_Chip", 10)
preCombat_sampsPC <- merge(samps_with_sex, preCombat_PCs, by = "row.names")
rownames(preCombat_sampsPC) <- preCombat_sampsPC[, "Row.names"]
all(rownames(preCombat_sampsPC) == preCombat_sampsPC[, "Row.names"]) #should be true


# calc pvalues
preCombat_pvals <- calc_pvals(pc.list, preCombat_sampsPC)

# generate heatmap
plot_heatmap(preCombat_pvals, "PreCombat")




####################################################################################
# Step 3: Run ComBat
####################################################################################

mod <- model.matrix(~as.factor(Exposure), data = samps_with_sex)
dim(mod) #  1 intercept, 1 Sex column

all(colnames(beta.imputed) == rownames(mod))
all(colnames(beta.imputed) == rownames(samps_with_sex))


# Run ComBat for Chip Effects
combat_beta <- ComBat(dat = beta.imputed, batch = samps_with_sex$Sentrix_ID, mod = mod)
save(combat_beta, file = "CICPT_1831_bkgd_beta_ComBAT_chipAdj_Exposure_Removed Cross Hybrid_CpGOnly.Rdata")


post_combat_BE(combat_data = combat_beta, data_wd_sex = samps_with_sex,
            col_name = "PostCombat_Chip", pcs = 6 )




# Run ComBat for Plate Effects
combat_beta <- ComBat(dat = combat_beta, batch = samps_with_sex$Sample_Plate, mod = mod)
save(combat_beta, file = "CICPT_1831_112_bkgd_beta_ComBat_step2_adjPos_Exposure_Removed Cross Hybrid_CpGOnly.Rdata")

min(combat_beta)
max(combat_beta)

post_combat_BE(combat_data = combat_beta, data_wd_sex = samps_with_sex,
               col_name = "PostCombat_Plate", pcs = 6 )




# Run ComBat for Position Effects
combat_beta <- ComBat(dat = combat_beta, batch = samps_with_sex$Sentrix_Position, mod = mod)
save(combat_beta, file = "CICPT_1831_112_bkgd_beta_ComBat_step2_adjPos_Exposure_Removed Cross Hybrid_CpGOnly.Rdata")

min(combat_beta)
max(combat_beta)

post_combat_BE(combat_data = combat_beta, data_wd_sex = samps_with_sex,
               col_name = "Chip_pos", pcs = 6 )




# Run ComBat for Sex Effects
## Not valid if adjusted for sex
samps_with_sex$Gender <- as.factor(samps_with_sex$Gender)
combat_beta <- ComBat(dat = combat_beta, batch = samps_with_sex$Gender, mod = mod)
save(combat_beta, file = "CICPT_1831_112_bkgd_beta_ComBat_step2_adjSex_Exposure_Removed Cross Hybrid_CpGOnly.Rdata")

min(combat_beta)
max(combat_beta)

post_combat_BE(combat_data = combat_beta, data_wd_sex = samps_with_sex,
               col_name = "Sex", pcs = 6 )




# from lumi package methylation_preprocess.R: beta <- 2^m/(2^m + 1)
reversbeta <- 2 ^ combat_beta / (2 ^ combat_beta + 1)
min(reversbeta) 
max(reversbeta) #1


# checking that that the rownames and column names are the same so that the missing index will be the same
all(rownames(reversbeta) == rownames(bmiq))
all(colnames(reversbeta) == colnames(bmiq))

reversbeta[missing] <- NA # putting missing values back in
sum(is.na(reversbeta[missing])) #150623
reversbetarnd <- round(reversbeta, digits = 4)

save(reversbeta, missing, file="CICPT_1831_bkgd_beta_postComBat_Removed Cross Hybrid_CpGOnly_notround.Rdata")
save(reversbetarnd, missing, file="CICPT_1831_bkgd_beta_postComBat_Removed Cross Hybrid_CpGOnly.Rdata")


#END

