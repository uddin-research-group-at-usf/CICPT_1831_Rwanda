
# Cross hybridized probes are in two seprate files 
# found in paper (http://www.sciencedirect.com/science/article/pii/S221359601630071X). 
# combine the data in a single file and then this file is used in another program 
# to remove these probes

files <- c("cross-hybridising CpG-targeting probes.txt", "cross-hybridising non-CpG-targeting probes.txt")

# read data
cross_hyb_data <- lapply(files, function(x) read.csv(x, header = FALSE))

# combine data
comb_cross_hyb_data <- do.call(rbind, cross_hyb_data)

#save
save(comb_cross_hyb_data, file =  "cross_hyb_data-CpG--nonCpG-targeting.Rdata")


