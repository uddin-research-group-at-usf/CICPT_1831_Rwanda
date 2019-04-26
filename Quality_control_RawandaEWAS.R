
library(data.table)
library(RPMM)


# Read methylation data from Genome Studio

methy_files <- c("CICPT-1831SampleMethy.txt") 
methy_data <- fread(methy_files, header = TRUE, sep = "\t", check.names = FALSE)
methy_data <- data.frame(methy_data, check.names = FALSE)


# All integer class to numeric
classes <- sapply(methy_data, class)
classes[which(classes == "integer")] <- "numeric"	


# Assign rownames
rownames(methy_data) <- methy_data[, "TargetID"]
head(methy_data)


# Percentage of missing values
mean(is.na(methy_data))


# Total number of missing values
sum(is.na(methy_data))



# Extract information from the main data file
extract_info <- function(ext_col, data, out_file){
  
  ext_data <- methy_data[, grep(ext_col, colnames(methy_data))]
  ext_data <- as.matrix(ext_data)
  
  if(ext_col == ".Detection")
    ext_col <- paste0(ext_col, ".Pval")
  
  colnames(ext_data) <- gsub(ext_col, "", colnames(ext_data))
  save(ext_data, file = out_file)
}



## Extract beta values
## avg data is saved to match with other data
## to make sure the col names and row names are same

# Extract beta values
extract_info(".AVG_Beta", methy_data, "CICPT_1831_beta_RAW.Rdata")


## Extract signal A and save
extract_info(".Signal_A", Methy_data, "CICPT_1831_signalA_RAW.Rdata")


## Extract signal B
extract_info(".Signal_B", Methy_data, "CICPT_1831_signalB_RAW.Rdata")


## Extract detection pvalues
extract_info(".Detection", Methy_data, "CICPT_1831_detectP_RAW.Rdata")



# load .Rdata file
# no need to load the data if you have run the above code in the same session
# if running this section of the code without running above code
# should have generated the data and saved

beta_vals = local(get(load('CICPT_1831_beta_RAW.Rdata')))

signal_A = local(get(load('CICPT_1831_signalA_RAW.Rdata')))

signal_B = local(get(load('CICPT_1831_signalB_RAW.Rdata')))

dect_pval = local(get(load('CICPT_1831_detectP_RAW.Rdata')))



# Checking the order of the matrices. These should all be true
if(!all(colnames(beta_vals) == colnames(signal_A)) |
   !all(rownames(beta_vals) == rownames(signal_A)) |
   !all(colnames(beta_vals) == colnames(signal_B)) |
   !all(rownames(beta_vals) == rownames(signal_B)) |
   !all(colnames(beta_vals) == colnames(dect_pval))|
   !all(rownames(beta_vals) == rownames(dect_pval)))
  stop("colnames names or row names do not match")




# ----------------------------------------------------------------------------- #

library(CpGassoc)


#remove probes with low signal intensity and too many missing values
# We set detection p-value threshold to 0.001
# We remove any CpGs with more than 10% missing data across the samples
beta.qc <- cpg.qc(beta.orig = beta_vals, siga = signal_A, 
                  sigb = signal_B, pval = dect_pval, p.cutoff = 0.001, cpg.miss = 0.1)

densityPlot(beta.qc, main = "Pre-Normalization CICPT-1831-Methylation EPIC-ARRAY 850K NS= 80")


# "Removed 0 samples with low signal"
# "Removed 6092 CpG sites with missing data for > 0.1 of samples"
save(beta.qc, file="CICPT-1831_beta_QC.Rdata")



# ----------------------------------------------------------------------------- #
library(wateRmelon)
library(FDb.InfiniumMethylation.hg19)


## Read annotation file from illumina
annot_file <- read.csv("MethylationEPIC_v-1-0_B4.csv", header = TRUE, skip = 7, row.names = 1)


## cross hybridization CpG and non-CPG 
cross_hyb_probes = local(get(load('cross_hyb_data-CpG--nonCpG-targeting.Rdata')))
colnames(cross_hyb_probes)[1 ] <- "TargetID"


## Probes to remove
rm_probes <- as.character(cross_hyb_probes[, "TargetID"])

sum(is.na(match(rm_probes, rownames(beta.qc)))) # 419 #Some may have already been removed in previous steps
rm_probes <- rm_probes[!is.na(match(rm_probes, rownames(beta.qc)))]
sum(is.na(match(rm_probes, rownames(beta.qc)))) # should be 0

beta.qc <- beta.qc[-match(rm_probes, rownames(beta.qc)),]




# ----------------------------------------------------------------------------- #
# Step 3: Creating design_v. for BMIQ normalization. From BMIQ documentation:
# "design_v: corresponding vector specifying probe design type (1=type1,2=type2). 
# This must be of the same length as beta.v and in the same order."

dim(annot_file) 
annot_file <- annot_file[rownames(beta.qc), ] #Keep only cg sites in beta.qc


beta.qc[which(beta.qc == 0)] <- 0.000001
beta.qc[which(beta.qc == 1)] <- 0.999999


#The matrix created below will contain all particpants and have identical cells across the matrix,
#indicating either "I" or "II" for type of probe. We will create this matrix then convert "I" to "1" and 
# "II" to "2". However, it is better when running BMIQ to have a vector with this information so we 
#will isolate 1 column (doesn't matter which) and verify it is a numeric vector and proceed.
#doing so will save us time. I tried a different way an got a lot of errors with NAs, so this code was 
#updated to speed up the processing time.

design_v <- matrix(nrow = nrow(beta.qc), ncol = ncol(beta.qc)) #create a matrix with set row and column names
rownames(design_v) <- rownames(beta.qc) # Make sure the row names match
colnames(design_v) <- colnames(beta.qc) # Make sure the column names match


## Probnames matching type I
type1_probes <- rownames(annot_file[annot_file[, "Infinium_Design_Type"] == "I", ])


## Flag with 1
#populate matrix with 1 into cells where the CpG sites are detected by Type I probes
design_v[type1_probes, ] <- 1 

min(design_v[type1_probes, ]) #1 #checking to make sure only 1 was inserted
max(design_v[type1_probes, ]) #1 #checking to make sure only 1 was inserted
sum(is.na(design_v[type1_probes, ])) #0 #check to see if there are any missing data points



type2_probes<-rownames(annot_file[annot_file[, "Infinium_Design_Type"] == "II", ])

#populate matrix with 2 into cells where the CpG sites are detected by Type II probes
design_v[type2_probes,] <- 2 

min(design_v[type2_probes, ]) # 2 # checking to make sure only 2 was inserted
max(design_v[type2_probes, ]) # 2 # checking to make sure only 2 was inserted
sum(is.na(design_v[type2_probes, ])) # 0 # check to see if there are any missing data points

#Convert the data matrix to a vector with the probe type in only 1 column. 
# pulling out only the first column (which column is called doesn't matter)
design_v2 <- design_v[, 1] 


class(design_v2) # check the class to make sure it is numeric
design_v <- design_v2 # convert to the variable we need design_v



# ----------------------------------------------------------------------------- #
# Step 4: Running BMIQ


# Checking order of beta and design matrices before running BMIQ
all(rownames(design_v) == rownames(beta.qc))
all(colnames(design_v) == colnames(beta.qc))
dim(design_v) 
dim(beta.qc) 

# if any na
sum(is.na(beta.qc)) 
length(design_v) 


sum(is.na(beta.qc)) 


#loop function to run BMIQ on each sample. sampleID is defined so the output should generate 
#actual participantID (sentrixID or BarcodeID) in the title of each plot

for(sample in 1:ncol(beta.qc)){
  tryCatch({
    beta.v <- beta.qc[, sample]
    samp_ID <- colnames(beta.qc)[sample]
    bmiqTemp <- BMIQ(beta.v = beta.v, design.v = design_v, sampleID = samp_ID) 
    beta.qc[, sample] <- bmiqTemp$nbeta
    message("Processing sample ", sample, ": ", colnames(beta.qc)[sample],  "... done\n" )
  },error = function(err){
    message("Normalization failed for sample ", sample, ": ", colnames(beta.qc)[sample] )
    message(err, "\n")
  }
  )
  }


bmiq_data <- beta.qc

# save
save(bmiq_data, design_v, file="CICPT_1831_beta_BMIQ_FULL_bySample_RemovedCrossHybrid.Rdata")


# plot
densityPlot(bmiq_data, main = "Normalized CICPT-1831-Methylation EPIC-ARRAY 850K NS= 80")





