
library(tibble)


# Scripts to perfomr discordant sex check


# Scientific to decimal
# Pay attention to this, as some R version append 0's to everthing
options(scipen = 999)


# read gender file
# read the sheet having phenotype information
gender <- read.csv(file.choose(), header = TRUE) # Sample_with_Age_Sex.csv


## Excel reads F as False
## change FALSE to F
gender[gender == FALSE] <- "F"
head(gender)


## extract sample and gender information
df1 <- gender[, c(2, 3, 5)]
df2 <- gender[, c(7, 8, 10)]

head(df1)
head(df2)

colnames(df1) <- c("Sample_Name", "Gender", "Age")
colnames(df2) <- c("Sample_Name", "Gender", "Age")

combined <- rbind(df1, df2)
head(combined)


# As original sample sheet and the output sample sheet has different notations of 
# samples names, to make them uniform, remove some special characters
replace_special <- function(data, col_name){

  data[[col_name]] <- gsub("/", "", data[[col_name]], fixed = TRUE)
  data[[col_name]] <- gsub(" ", "", data[[col_name]], fixed = TRUE)
  data[[col_name]] <- gsub("(", "", data[[col_name]], fixed = TRUE)
  data[[col_name]] <- gsub(")", "", data[[col_name]], fixed = TRUE)
  data[[col_name]] <- gsub("-", "", data[[col_name]], fixed = TRUE)
  return(data)
}


original_sheet <- replace_special(combined[, 1:2], "Sample_Name")
head(original_sheet)


# read sample sheet to match with
# orininal sample sheet containg the sentrix ids and samples
# find which samples are matching and pull out
sample_sheet <- read.csv(file.choose()) # UddinMethylationSampleSheet2019.csv-2.csv
head(sample_sheet)
form_sample_sheet <- replace_special(sample_sheet, "Sample_Name")

head(form_sample_sheet)


# After making the two sample sheets unifrom, we need to match to find if all samples in 
# original sample sheet are matching with the data sample sheet (on which chips are run)
## find match, and take the intersection (the samples which are actually run on chips)

matching_samp <- original_sheet[original_sheet$Sample_Name %in% form_sample_sheet$Sample_Name, ]


# if any samples are not matching
not_matching_samp <- original_sheet[!original_sheet$Sample_Name %in% form_sample_sheet$Sample_Name, ]


## finding unique samples, not to allow any duplicates
unique(matching_samp$Sample_Name)


## duplicates
if(any(duplicated(matching_samp[, 1]))){
  duplicate_samps <- matching_samp[duplicated(matching_samp[, 1]),]
  warning("some duplicate samples detected, ", duplicate_samps )
}
  
  
## removing the duplicated one
matching_samp_wo_dupli <- matching_samp[!duplicated(matching_samp[, 1]), ]



## some samples are form DNHS data
## load phenotype data of DNHS
DNHS_samples <- read.csv(file.choose())

# extract DNHS samples from the smaple sheet
DNHS_8_samp <- form_sample_sheet[73:80, ]

# extract info of these 8 samples from DNHS original sample sheet
DNHS_samps_gender <- DNHS_samples[match(DNHS_8_samp$Sample_Name, DNHS_samples$BloodID), c(3,5)]


DNHS_samps_gender$Gender <- as.character(DNHS_samps_gender$Gender)
DNHS_samps_gender$Gender[DNHS_samps_gender$Gender == "Female"] <- "F"
DNHS_samps_gender$Gender[DNHS_samps_gender$Gender == "Male"] <- "M"


colnames(DNHS_samps_gender) <- c("Sample_Name", "Gender")
DNHS_samps_gender

## combine all rawanda and DNHS samples
all_samples <- rbind(matching_samp_wo_dupli, DNHS_samps_gender)



## contorls samples are from males, so add controls to the data
control_samples <- form_sample_sheet[grepl("PooledMaleDNA", 
                                           form_sample_sheet$Sample_Name, ignore.case = TRUE), ]
control_samples$Gender <- c("M", "M")
control_samples <- control_samples[, c("Sample_Name", "Gender")]



## combine all samples and control samples
all_eighty_samples <- rbind(all_samples, control_samples)


# merge all_eighty_samples and the output data sample sheet to get 
# Sentrix_ID to match it with output of minfi for discordant check
combined_data <- merge(all_eighty_samples, form_sample_sheet, by = "Sample_Name")

if(nrow(all_eighty_samples) != nrow(combined_data)){
  warning("Some duplicates are introducted while merging data")
}

# find and remove the duplicates
combined_data <- combined_data[1:80, ]

# new column to make it compatable with sample names to find match sex match
combined_data$sentid_sentpos <- paste0(combined_data$Sentrix_ID, "_", combined_data$Sentrix_Position)
save(combined_data , file = "CICPT_1831_data_with_gender.Rdata")




# Now when data is prepared for discordant check
# load data and run the following code
combined_data <- local(get(load("CICPT_1831_data_with_gender.Rdata")))



# ---------------------------------------------------------------------------- #

library(minfi)

# read idat files, put in a directory and read recursivelyset
# set directory to the location of idat files 
RGset <- read.metharray.exp("Data", recursive = T)


# if any pehnotype data, otherwise it will not display anything
pd <- pData(RGset)
head(pd)



# The RGChannelSet stores also a manifest object that contains the 
# probe design information of the array:
manifest <- getManifest(RGset)
manifest

head(getProbeInfo(manifest))



#A MethylSet objects contains only the methylated and unmethylated signals. You create this by

MSet <- preprocessRaw(RGset)
MSet



# The accessors getMeth and getUnmeth can be used to get the
# methylated and unmethylated intensities matrices:
head(getMeth(MSet)[,1:3])
head(getUnmeth(MSet)[,1:3])


# The functions getBeta, getM and getCN return respectively the 
# Beta value matrix, M value matrix and the Copy Number matrix.

RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
RSet


beta <- getBeta(RSet)



GRset <- mapToGenome(RGset)
GRset


# To access the full annotation, one can use the command getAnnotation:
annotation <- getAnnotation(GRset)
names(annotation)



#The functions getQC and plotQC are designed to extract and plot the quality 
# control information from the MethylSet:

qc <- getQC(MSet)
head(qc)
plotQC(qc)




## Sex prediction
predictedSex <- getSex(GRset, )$predictedSex
head(predictedSex)

# Total males and females
sum(predictedSex == "F")
sum(predictedSex == "M")

sampleNames <- sampleNames(GRset)
probeNames <- featureNames(GRset)
pheno <- pData(GRset)

pred_sex_data <- data.frame("Sample_Name" = sampleNames, "gender" = predictedSex)
save(pred_sex_data, file =  "Sex_prediction_minfi.Rdata")


## split to group the chips and find ratio of gender
library(dplyr)
library(tidyr)
library(plyr)

sex_data <- pred_sex_data %>% separate(Sample_Name, c("Sample", "chip"), "_")

sex_data_mf <- ddply(sex_data, .(Sample), summarise, nfml = length(gender[gender == "F"]), 
                     nmale = length(gender[gender == "M"]))



# plot male and female distribution on the chips
boxplot(sex_data_mf$nfml, sex_data_mf$nmale,
        horizontal=TRUE,
        names=c("Females","Males"),
        col=c("turquoise","tomato"),
        xlab=" Number of Males and Females",
        main="Distribution of males and females on chips")



## matching the data for discordant sex
sex_pheno <- merge(combined_data, pred_sex_data, 
                   by.x = "sentid_sentpos", by.y = "Sample_Name", all.x = T)
View(sex_pheno[, c("sentid_sentpos", "Sample_Name", "Gender", "gender")])


## mismatcches
mis_match <- sex_pheno[sex_pheno$Gender != sex_pheno$gender, ]
mis_match_subset <- mis_match[, c("sentid_sentpos", "Sample_Name", "Gender", "gender")]
write.csv(mis_match, "Sex_mismatch_samples.csv")

# rename columns 
colnames(mis_match_subset)[3] <- "Gender_input_file"
colnames(mis_match_subset)[4] <- "Gender_output_file"

# write subset of informtion coltaining 
# sentid_sentpos", "Sample_Name", "Gender", "gender"
write.csv(mis_match_subset, "Sex_mismatch_samples_subset.csv")



# number of males and females in input file
table(sex_pheno$Gender)


# number of males and females in output file
table(sex_pheno$gender)



# To choose the cutoff to separate the two gender clusters, one can plot med(Y) against med(Y)
# with the function plotSex:
#plotSex(getSex(GRset, cutoff = -2))

estSex <- getSex(GRset, cutoff = -2)

plot(
  x = estSex$xMed,
  y = estSex$yMed,
  type = "n",
  xlab = "X chr, median total intensity (log2)",
  ylab = "Y chr, median total intensity (log2)")
text(
  x = estSex$xMed,
  y = estSex$yMed,
  col = ifelse(estSex$predictedSex == "M", "deepskyblue", "deeppink3"))
legend(
  "bottomleft",
  c("M", "F"),
  col = c("deepskyblue", "deeppink3"),
  pch = 16)



# ------------------------------------------------------------------------------------#
## find how many samples of discordant sex are in the flagged file based on lab metrics

# read flagged file
flagged_file <- read.csv(file.choose()) # Microarray_addition_3_14_2019_ZG_MU (Wani, Agaz).csv
flagged_samps <- replace_special(flagged_file, "Sample.Name")

# Intersection of samples that are present in discordant sex and flagged file
flagged_DSex_inters <- flagged_samps[!is.na(match(mis_match$Sample_Name, samp_names$Sample.Name)), ]
write.csv(flagged_DSex_inters, "Samples found in discordant sex and flagged fil.csv")


# to bind with discordant sex
flagged_two_col <- flagged_samps[, c("No.", "Sample.Name")]
colnames(flagged_two_col) <- c("S.No", "Sample_Name")


# Union of samples in discordant sex and flagged file
discordant_sex_mm <- rownames_to_column(mis_match_subset[, c(1, 2)]) 
colnames(discordant_sex_mm)[1] <- "S.No"
discordant_sex_mm <- discordant_sex_mm[, c("S.No", "Sample_Name")]

flag_DSex_together <- unique(rbind(flagged_two_col, discordant_sex_mm))
write.csv(flag_DSex_together, "Samples (union) in discordant sex and flagged file.csv" )


# End
