#BiocManager::install("minfi", version = "3.8")
#BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", version = "3.8")
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)



# Set the path to idat files directory and read the files
# in recursive fashion
RGset <- read.metharray.exp("/Users/hussainwani/Rawanda_CICPT_1831/Data", 
                            recursive = TRUE, verbose = TRUE)



## Preprocess noob
RGset_noob <- preprocessNoob(RGset, offset = 15, dyeCorr = TRUE, 
                             verbose = TRUE, dyeMethod="single")
save(RGset_noob, file = "RGset_noob.Rdata")


# Extract pval, beta, methylated and unmethylated signal files
# to extract pvalues, 'RGChannelSEt'/'RGChannelSetExtended' is needed
pval <- detectionP(RGset)
beta <- getBeta(RGset_noob)
signalA <- getUnmeth(RGset_noob)
signalB <- getMeth(RGset_noob)


# plot beta output of preprocessNoob
batches <-  ceiling(length(colnames(beta)) / 8)

samp_num <- 1
for (i in 1:batches) {
  pdf(paste0("Density plot of chip", i,  "samples.pdf"))
  
  if(i == 1)
    chip_st <- 1
  
  # 8 samples per chip
  chip_stp <- i * 8
  densityBeanPlot(beta[, chip_st:chip_stp], 
                  main =  paste("Density plot of chip", i,  "samples"))
  
  # add sample index to chip samples
  # 8 samples on each chip
  # x-axis, y-axis and number to write
  for (i in 1:8) {
    text(1, i, samp_num)
    samp_num <-  samp_num + 1
  }

  dev.off()
  chip_st <- chip_stp + 1
}


# save pval, beta, signalA and SignalB
save(pval, beta, signalA, signalB, file = "PreprocessNoobOutput.RData")


# clear objects
rm(list=ls())


## Quality control using CpGassoc
library(CpGassoc)

# run this line if you want to load the saved data
# otherwise running about code should should generate the necessary data
load("PreprocessNoobOutput.RData")

# run QC to remove CpG sites
beta.qcd <- cpg.qc(beta, signalA, signalB, pval, 
                   p.cutoff = .01, cpg.miss = .1, sample.miss = .1)

write.csv(beta.qcd, file = "preprocessNoob_qcd.csv", quote = FALSE,row.names = TRUE)


# check if any sample is removed
if(length(colnames(beta) != length(colnames(beta.qcd)))){
  indx <- which(colnames(beta) %in% colnames(beta.qcd) == FALSE)
  warning("Samples at index: ", indx, " removed" )
}
  

# Remove cross reactive probes: 
# Got the cross reactive probes from 
# a paper (http://www.sciencedirect.com/science/article/pii/S221359601630071X). 
# Processing QCd data to remove only Cross Reactive Probes
# load data which is already prepared by MkCrossHybridData.R

cross_hyb_data <- local(get(load("cross_hyb_data-CpG--nonCpG-targeting.Rdata")))
beta_aftr_rm_crosshbd <- beta.qcd[!(row.names(beta.qcd) %in% cross_hyb_data), ]

write.csv(beta_aftr_rm_crosshbd, 
          file ="noob_qcd_crossReactiveProbesRemoved.csv", quote = FALSE, row.names = TRUE)
