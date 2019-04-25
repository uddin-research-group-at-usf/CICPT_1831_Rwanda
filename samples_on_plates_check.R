
# Intension of writing this code is to check if all the
# samples are at desired location on the chip.
# To make sure there is no mismatch between the samples 
# given for platting and output after running the data on chips


# load input sample sheet given to run on chips
# this sample sheet is prepared to run data on chips
samples <-  read.csv(file.choose(), header = TRUE)
head(sample)

## extract sample.name and plate.well information
samples_subs <- samples[, c("Sample.Name", "Plate.Well")]
head(samples_subs)



# load output sample sheet after running chips
# this sample sheet will be with the idat files 
# after running chips
sample_sheet <- read.csv(file.choose(), header = TRUE)
head(sample_sheet)


# the input and output sample sheets have sometimes different notations of 
# representing sampples, make them uniform by removing any special characters
# in sample names
replace_special <- function(data, col_name){

  data[[col_name]] <- gsub("/", "", data[[col_name]], fixed = TRUE)
  data[[col_name]] <- gsub(" ", "", data[[col_name]], fixed = TRUE)
  data[[col_name]] <- gsub("(", "", data[[col_name]], fixed = TRUE)
  data[[col_name]] <- gsub(")", "", data[[col_name]], fixed = TRUE)
  data[[col_name]] <- gsub("-", "", data[[col_name]], fixed = TRUE)
  return(data)
}


# replace special characters to make it uniform
original_sheet <- replace_special(samples_subs, "Sample.Name")
original_sheet$Plate.Well <- as.character(original_sheet$Plate.Well)


# if any samples are empty, remove those
empty_ids <- which(original_sheet$Sample.Name == "")
original_sheet <- original_sheet[-empty_ids, ]
colnames(original_sheet) <- c("Sample_Name", "Plate_Well")


# replace '0' with '-' and then remove '-'
substr(original_sheet$Plate.Well[1:72], 2, 2) <- '-'
original_sheet$Plate.Well <- gsub("-", "", original_sheet$Plate.Well)

# replace in output sample sheet
out_samp_sheet <- replace_special(sample_sheet, "Sample_Name")
head(out_samp_sheet)


# merge by sample well / plate well
merged_dat <- unique(merge(out_samp_sheet, original_sheet, by.x = "Sample_Well", by.y = "Plate_Well"))

# if true all samples are at correct place
if(!all(merged_dat$Sample_Name.x == merged_dat$Sample_Name.y))
  stop("Samples are not at desired location on chips")

