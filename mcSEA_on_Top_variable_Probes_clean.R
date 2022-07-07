#' Code to perform DMR analysis
#'

library(feather)
library(tibble)
library(mCSEA)
library(xlsx)
library(ggpubr)
library(clusterProfiler)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(data.table)

set.seed(42)

beta_vals <- read_feather("E:/Rawanda EWAS/Data/Rwanda_Noob_QCd_Combat_adj.feather")
beta_vals <- column_to_rownames(beta_vals, var = 'rowname')
beta_vals <- as.data.frame(beta_vals, stringsAsFactors = FALSE)
dim(beta_vals)

# we will select the top cpgs based on variance
# First we need to calculated the variance for each row (cpg)
get_variance <- function(beta){
  beta$var <- apply(beta, 1, var)
  x <- beta[order(beta$var, decreasing = TRUE), ]
  return(x)
}


# As we want top 5% probes with high variance
# Now select the top cpgs based on variance
get_top <- function(beta, per){
  top <- (nrow(beta) * per)/100
  x <- beta[1:top, ]
  return(x)
}


# Call functions, using only the top 5% probes
beta_wd_var <- get_variance(beta = beta_vals)
beta_top <- get_top(beta = beta_wd_var, per = 5)
# max(beta_top$var)
# min(beta_top$var)

# In the review, it was asked to perform analysis using
# all the probes
# beta_top = beta_vals


dim(beta_top)


#free some space
# rm(beta_vals)
# gc()


# function to subset the data
subset_data <- function(df, pat1, pat2){
  x <- df[which(grepl(pat1, df[[pat2]])), ]
  return(x)
}


# subset beta values
subset_beta <- function(df, pheno){
  x <- df[, which(colnames(df) %in% rownames(pheno))]
  stopifnot(all(colnames(x) == rownames(pheno)))
  return(x)
}


# pull out required columns
get_required_cols <- function(pheno, cols){
  x <- pheno[, cols]
  return(x)
}


# Run the test
# we will look for DMRs in promoter regions
run_mCSEA <- function(b_rank, beta, pheno){
  x <- mCSEATest(b_rank, beta, pheno,
                 regionsTypes = "promoters", platform = "EPIC",
                 nproc = 10, minCpGs = 5)
  return(x)
}



# Select those will padj < 0.05
get_significant <- function(results){
  x <- results$promoters[results$promoters$padj < 0.05, ]
  return(x)
}


# Function to plot the significant genes
# the output of mCSEA
# it will save pdfs
plot_genes <- function(path, results, genelis){
  setwd(path)
  for(i in 1:length(genelis)){
    g <- genelis[i]
    message('Saving plot for gene ', g)
    output <- tryCatch({
      mCSEAPlot(results, regionType = "promoters",
                dmrName = g, leadingEdge = TRUE,
                CGI = TRUE,
                transcriptAnnotation = "symbol", makePDF = TRUE)
    }, error = function(err){
      message("Gene not ploted ", g)
    })
  }
}


# As we get the cpgs from mCSEA into a vector
# Function to split the cpgs of one gene
clean_cpgs <- function(cpgs){
  x <- unlist(strsplit(cpgs, ","))
  x <- gsub(" ", "", x, fixed = TRUE) # remove any spaces
  return(x)
}


# Function to plot cpgs
plot_cpgs <- function(b_vals, g_name=NULL, s_name ){
  pl <- list()
  cpgs <- colnames(b_vals)[which(grepl('^cg', colnames(b_vals)))]
  for(i in 1:length(cpgs)){
    p <- ggboxplot(b_vals, x = "Exposure", y = cpgs[i],
                   color = "Exposure", palette = c("#0072B5FF", "#BC3C29FF"),# "nejm",
                   size = 1,
                   add = "dotplot"
                   )
    #  Add p-value
    p <- p + stat_compare_means(method ="t.test", label = "p.format") +
      font("title", size = 20)+ #face = "bold"
      font("xlab", size = 16)+
      font("ylab", size = 16)+
      font("xy.text", size = 13)+
      font("caption", size =15, face= "bold")+
      theme(legend.position = "none",
            axis.title.x =element_blank(),
            axis.title.y =element_blank())

    pl[[length(pl) + 1]] <- p
  }
  names(pl) <- cpgs
  return(pl)
}



# load phenotype data
pheno <- read.csv("E:/Rawanda EWAS/Data/pheno_final_59.csv")
pheno <- column_to_rownames(pheno, var = 'Row.names') # convert column to rownames
pheno$Exposure <- ifelse(pheno$Exposure == 1, "Yes", "No")
pheno$Exposure <- as.factor(pheno$Exposure)
dim(pheno)


# get children and mother separate
children <- subset_data(pheno, pat1 = 'chd', pat2 = 'Group')
mothers <- subset_data(pheno, pat1 = 'moth', pat2 = 'Group')
table(children$Exposure)
table(mothers$Exposure)

# Check how many males and females are in exposed and unexposed in children
# Mothers will be all females
ch_exp <- subset_data(df = children, pat1 = "exp", pat2 = "Group")
ch_contrl <- subset_data(df = children, pat1 = "con", pat2 = "Group")
dim(ch_exp)
dim(ch_contrl)
table(ch_exp$Gender)
table(ch_contrl$Gender)

get_sumry <- function(beta_df){
  print(dim(beta_df))
  print("Min: ")
  print(min(beta_df, na.rm = T))
  print("Max: ")
  print(max(beta_df, na.rm = T))
  col_mean <- colMeans(beta_df, na.rm = T)
  mean_ch <- sum(col_mean)/length(col_mean)
  print("mean: ")
  print(mean_ch)
  return(mean_ch)
}


# Now get the beta values for children and mothers
beta_children <- subset_beta(df = beta_top,
                             pheno = children)
dim(beta_children)
# exposed children beta values
exp_ch_beta <- beta_children[, which(colnames(beta_children) %in% rownames(ch_exp))]
con_ch_beta <- beta_children[, which(colnames(beta_children) %in% rownames(ch_contrl))]
dim(beta_children)
dim(exp_ch_beta)
get_sumry(beta_df = beta_children)
sink("E:/Rawanda EWAS/Data/Mean_methylation_levels.txt")
print("Exposed child")
get_sumry(beta_df = exp_ch_beta)

print("Control child")
get_sumry(beta_df = con_ch_beta)



# Mothers
exp_moth <- subset_data(df = mothers, pat1 = "exp", pat2 = "Group")
con_moth <- subset_data(df = mothers, pat1 = "con", pat2 = "Group")
dim(exp_moth)
dim(con_moth)

beta_mothers <- subset_beta(df = beta_top,
                            pheno = mothers)
dim(beta_mothers)
exp_mo_beta <- beta_mothers[, which(colnames(beta_mothers) %in% rownames(exp_moth))]
con_mo_beta <- beta_mothers[, which(colnames(beta_mothers) %in% rownames(con_moth))]
dim(exp_mo_beta)
dim(con_mo_beta)
dim(beta_mothers)

print("Exposed mothers")
# get_sumry(beta_df = beta_mothers)
get_sumry(beta_df = exp_mo_beta)

print("Control mothers")
get_sumry(beta_df = con_mo_beta)

sink()

# Now mean methylation of exposed vs un-exposed
exp <- subset_data(pheno, pat1 = 'exp', pat2 = 'Group')
cont <- subset_data(pheno, pat1 = 'cont', pat2 = 'Group')
table(exp$Exposure)
table(cont$Exposure)
exp_beta <- beta_top[, which(colnames(beta_top) %in% rownames(exp))]
cont_beta <- beta_top[, which(colnames(beta_top) %in% rownames(cont))]
dim(exp_beta)
dim(cont_beta)
get_sumry(beta_df = exp_beta)
get_sumry(beta_df = cont_beta)


# Check if cell types sum to 1
cell_sum <- pheno %>%
  select(c("CD8T", "CD4T", "NK", "Bcell",
                             "Mono", "Neu")) %>%
  rowSums()


# Get only the information we need
req_cols <- c("Exposure", "CD8T", "CD4T", "NK", "Bcell",
              "Mono", "Neu", "Gender") # , "PTSD_Tot", "BDI", "Neu"
pheno_children <- get_required_cols(pheno = children, cols =  req_cols)
head(pheno_children)

pheno_mothers <- get_required_cols(pheno = mothers, cols = req_cols)
head(pheno_mothers)


# ----------------------------------- Run the below code only one time
# If the results are saved, skip this part

# Check if information in beta values and phenotype is matching
all(colnames(beta_children) == rownames(pheno_children))
all(colnames(beta_mothers) == rownames(pheno_mothers))


# Rank the probes in mothers and children
myRank_ch <- rankProbes(beta_children, pheno_children,
                        refGroup = "Yes", typeInput = "beta")
myRank_mo <- rankProbes(beta_mothers, pheno_mothers,
                        refGroup = "Yes", typeInput = "beta")


# Now run the test and select the significant
# First for children
myResults_ch <- run_mCSEA(b_rank = myRank_ch, beta = beta_children,
                          pheno = pheno_children)
myResults_ch_sig <- get_significant(results = myResults_ch) # < 0.05
dim(myResults_ch_sig)
dim(beta_children)


# For mothers
myResults_mo <- run_mCSEA(b_rank = myRank_mo, beta = beta_mothers,
                          pheno = pheno_mothers)
myResults_mo_sig <- get_significant(results = myResults_mo)
dim(myResults_mo_sig)
dim(beta_mothers)

head(myResults_mo_sig)

# Common genes in children and mothers
rownames(myResults_ch_sig)[rownames(myResults_ch_sig) %in%
                             rownames(myResults_mo_sig)]




# Plot genes functions and save it
path <- "E:/Rawanda EWAS/Plots/mCSEA/Children/"
plot_genes(path = path, results = myResults_ch,
           genelis = rownames(myResults_ch_sig))

path1 <- "E:/Rawanda EWAS/Plots/mCSEA/Mothers/"
plot_genes(path = path1, results = myResults_mo,
           genelis = rownames(myResults_mo_sig))


# This was for one gene
mCSEAPlot(myResults_mo, regionType = "promoters",
          dmrName = "FGFR2", leadingEdge = TRUE,
          CGI = TRUE, transcriptAnnotation = "symbol",
          genes = TRUE, makePDF = FALSE)



# write.xlsx(myResults_ch_sig, file = "E:/Rawanda EWAS/mCSEA Output/Children.xlsx",
#            sheetName = "Children", col.names = TRUE,
#            row.names = TRUE)


write.xlsx(myResults_ch_sig, file = "E:/Rawanda EWAS/mCSEA Output/Children_using_all_probes.xlsx",
           sheetName = "Children", col.names = TRUE,
           row.names = TRUE)

# write.xlsx(myResults_mo_sig, file = "E:/Rawanda EWAS/mCSEA Output/Mothers.xlsx",
#            sheetName = "Mothers", col.names = TRUE,
#            row.names = TRUE)

write.xlsx(myResults_mo_sig, file = "E:/Rawanda EWAS/mCSEA Output/Mothers_using_all_probes.xlsx",
           sheetName = "Mothers", col.names = TRUE,
           row.names = TRUE)


# common in children and mothers combined
common_ch_mo <- merge(myResults_ch_sig, myResults_mo_sig,
                      by = 'row.names')

write.xlsx(common_ch_mo, file = "E:/Rawanda EWAS/mCSEA Output/Common Children and Mothers_allProbes.xlsx",
           sheetName = "children&Mothers", col.names = TRUE,
           row.names = FALSE)


# Now check how many are common between the analysis using top 5%
# and all the probes
ch <- read.xlsx("E:/Rawanda EWAS/mCSEA Output/Children.xlsx", sheetIndex = 1)
mo <- read.xlsx("E:/Rawanda EWAS/mCSEA Output/Mothers.xlsx", sheetIndex = 1)
co <- read.xlsx("E:/Rawanda EWAS/mCSEA Output/Common Children and Mothers.xlsx", sheetIndex = 1)


# we need to sort
ch <- ch[order(ch$pval), ]
mo <- mo[order(mo$pval), ]


# As asked in the revisions to get the whole DNAm mean difference
# of the DMRs between cases and controls
length(ch$leadingEdge)
cal_mean <- function(df, df1, df2){
  cpgs <- df$leadingEdge
  mean_vals <- list()
  t_vals <- list()
  p_vals <- list()
  for (i in 1:length(cpgs)){
    x <- strsplit(cpgs[[i]], split = ",")
    x <- trimws(x[[1]])
    print(x)
    a <- df1[which(rownames(df1) %in% x), ] # exposed
    b <- df2[which(rownames(df2) %in% x), ] # control
    message("Exposed")
    # m_e <- get_sumry(beta_df = a)
    a_col_s <- colSums(a, na.rm = T)/ nrow(a)

    message("Unexposed")
    # m_ue <- get_sumry(beta_df = b)
    b_col_s <- colSums(b, na.rm = T)/ nrow(b)

    # dif <- m_e - m_ue
    t_test <- t.test(a_col_s, b_col_s)
    # mean_vals[i] <- dif
    t_vals[i] <- t_test$statistic
    p_vals[i] <- t_test$p.value
  }
  # names(mean_vals) <- df$Gene
   # df_r <- data.frame("DMR" = df$Gene, "Mean_diff" = unlist(mean_vals ))
  df_r <- data.frame("DMR" = df$Gene, "t" = unlist(t_vals ),
                     "p.value" = unlist(p_vals))
  return(df_r)
  }



mean_diff_ch <- cal_mean(df = ch, df1 = exp_ch_beta,
         df2 = con_ch_beta)

mean_diff_mo <- cal_mean(df = mo, df1 = exp_mo_beta,
                         df2 = con_mo_beta)

mean_diff_ch
write.csv(mean_diff_ch, "E:/Rawanda EWAS/Data/Children_DMRs_meth_difference.csv")

mean_diff_mo
write.csv(mean_diff_mo, "E:/Rawanda EWAS/Data/Mothers_DMRs_meth_difference.csv")

# update the dmrs with mean difference
upd_ch_exp <- merge(ch, mean_diff_ch,
                    by.x = "Gene",
                    by.y = "DMR")
dim(upd_ch_exp)
write.xlsx(upd_ch_exp, "E:/Rawanda EWAS/Data/Children_with_meth_difference.xlsx")

upd_mo_exp <- merge(mo, mean_diff_mo,
                    by.x = "Gene",
                    by.y = "DMR")
dim(upd_mo_exp)

write.xlsx(upd_mo_exp, "E:/Rawanda EWAS/Data/Mothers_with_meth_difference.xlsx")




# also check whole methylation difference between exp vs unexp
# sum across samples
# exp_rsum <- exp_beta %>% mutate(row_sum =
#                                   rowSums(across(everything()),
#                                           na.rm = T)/ncol(.))
# cont_rsum <- cont_beta %>% mutate(row_sum =
#                                     rowSums(across(everything()),
#                                             na.rm = T)/ncol(.))
#
# exp_cont_beta <- data.frame("Cpgs" = rownames(exp_rsum),
#                             "exp_sum" = exp_rsum$row_sum,
#                             "cont_sum" = cont_rsum$row_sum)
# head(exp_cont_beta)
#
# cont_vs_exp <- t.test(exp_cont_beta$exp_sum, exp_cont_beta$cont_sum,
#                       conf.level = .95)
# cont_vs_exp

# sum across samples
exp_col_s <- colSums(exp_beta, na.rm = T)/ nrow(exp_beta)
cont_col_s <- colSums(cont_beta, na.rm = T)/ nrow(exp_beta)
t.test(exp_col_s, cont_col_s)

# how many genes are  present
table(ch$Gene %in% rownames(myResults_ch_sig))
table(mo$Gene %in% rownames(myResults_mo_sig))
table(co$Gene %in% common_ch_mo$Row.names)

# which genes are not present
ch$Gene[which(!ch$Gene %in% rownames(myResults_ch_sig))]
mo$Gene[which(!mo$Gene %in% rownames(myResults_mo_sig))]


# -------------------The plots generated for presentation -------------------
# Load the saved results

path <- "E:/Rawanda EWAS/mCSEA Output/"
myResults_ch_sig <- read.xlsx(paste0(path, "Children.xlsx"), sheetIndex = 1)
View(myResults_ch_sig)

myResults_mo_sig <- read.xlsx(paste0(path, "Mothers.xlsx"), sheetIndex = 1)
View(myResults_mo_sig)

common_ch_mo <- read.xlsx(paste0(path, "Common Children and Mothers.xlsx"), sheetIndex = 1)
View(common_ch_mo)

# get chromosome locations
annotation_file <- fread("E:/Rawanda EWAS/Data/MethylationEPIC_v-1-0_B4.csv",
                         data.table = F, skip = 7,
                         fill = T)
dim(annotation_file)

View(head(annotation_file))

# Get the information of the genes that are matching
get_chr_info <- function(annot_file, gene_names){

  gene_names_ptr <- c(paste0("^", gene_names,";"), paste0("^", gene_names,"$"),
                  paste0(";", gene_names,"$")) # create pattern to get all combinations

  # Now get all the patterns from the annotation file
  anno <- annot_file[which(grepl(paste(c(gene_names_ptr), collapse = "|"),
                                         ignore.case = T, annot_file$UCSC_RefGene_Name)), ]
  anno <- anno[which(!duplicated(anno$UCSC_RefGene_Name)), ]
  anno <- anno[, c('CHR', 'UCSC_RefGene_Name')]

  unq_genes <- lapply(strsplit(anno$UCSC_RefGene_Name, ';'), unique) # get unique

  anno$gene_anno <- lapply(unq_genes, function(x) x[x %in% gene_names][1]) # get those matching the sig list
  anno <- anno[which(!duplicated(anno$gene_anno)), ]

}

# function call
anno_ch <- get_chr_info(annot_file = annotation_file, gene_names = myResults_ch_sig$Gene)
myResults_ch_sig <- merge(myResults_ch_sig, anno_ch[, c(1,3)],
                          by.x = 'Gene', by.y = 'gene_anno')


View(myResults_ch_sig)
dim(myResults_ch_sig)
write.xlsx(myResults_ch_sig[, c(1,9,2:8)], paste0(path, "Children_with_chr.xlsx"),
           row.names = F)


anno_mo <- get_chr_info(annot_file = annotation_file, gene_names = myResults_mo_sig$Gene)
# annot_file <- annotation_file
# gene_names <- myResults_mo_sig$Gene

View(anno_mo)

myResults_mo_sig <- merge(myResults_mo_sig, anno_mo[, c(1,3)],
                          by.x = 'Gene', by.y = 'gene_anno')

View(myResults_mo_sig)
write.xlsx(myResults_mo_sig[, c(1,9,2:8)], paste0(path, "Mothers_with_chr.xlsx"),
           row.names = F)

# function to subset cpgs
subset_cpgs <- function(b_vals, cpg_lis, pheno){
  x <- b_vals[which(rownames(b_vals) %in% cpg_lis), ]
  x <- merge(t(x), pheno, by='row.names')
  return(x)
}


# Now get the genes to plot
# We will plot two genes
genes_to_plot <- c(common_ch_mo$Gene[1], common_ch_mo$Gene[5],
                   common_ch_mo$Gene[6]
                   )


# Now get the cpgs of common genes in children and mothers
# Now we will get the cpgs to plot from the significant ones

# children
cpgs_com_ch <- lapply(seq_along(genes_to_plot), function(i){
  clean_cpgs(cpgs = common_ch_mo$leadingEdge.x[
    common_ch_mo$Gene == genes_to_plot[i]])
})
names(cpgs_com_ch) <- genes_to_plot


# Mothers
cpgs_com_mo <- lapply(seq_along(genes_to_plot), function(i){
  clean_cpgs(cpgs = common_ch_mo$leadingEdge.y[
    common_ch_mo$Gene == genes_to_plot[i]])
})
names(cpgs_com_mo) <- genes_to_plot




# Now get the cpgs from the data
df_com_ch <- lapply(seq_along(cpgs_com_ch), function(i){
  subset_cpgs(b_vals = beta_children,
              cpg_lis = cpgs_com_ch[[i]],pheno = pheno_children )

})
names(df_com_ch) <- genes_to_plot


df_com_mo <- lapply(seq_along(cpgs_com_mo), function(i){
  subset_cpgs(b_vals = beta_mothers,
              cpg_lis = cpgs_com_mo[[i]],pheno = pheno_mothers )

})
names(df_com_mo) <- genes_to_plot



# Test if we have all cpgs
Map(function(x, y) x %in% colnames(y),
    cpgs_com_ch, df_com_ch)

Map(function(x, y) x %in% colnames(y),
    cpgs_com_mo, df_com_mo)



# Children

com_ch_plot <- lapply(seq_along(genes_to_plot), function(i){
  plot_cpgs(b_vals = df_com_ch[[i]],
            g_name = genes_to_plot[i],
            s_name = "Children")
})
names(com_ch_plot) <- genes_to_plot


# Mothers
com_mo_plot <- lapply(seq_along(genes_to_plot), function(i){
  plot_cpgs(b_vals = df_com_mo[[i]],
            g_name = genes_to_plot[i],
            s_name = "Mothers")
})
names(com_mo_plot) <- genes_to_plot

common_cpgs <- names(com_ch_plot$BCOR)[which(names(com_ch_plot$BCOR) %in% names(com_mo_plot$BCOR))]

# Cpgs to plot based on consistency with the estimated marginal means (janelle)
# BCOR:  cg15627188, cg02931660, cg02932805:
# PRDM8: cg18073471: , cg05059566, cg06307913:
# VWDE: cg03579179,cg06484146, cg20607287).

plot_grid(com_ch_plot$BCOR$cg15627188, com_mo_plot$BCOR$cg15627188) # for gene BCOR
plot_grid(com_ch_plot$BCOR$cg02931660, com_mo_plot$BCOR$cg02931660) # for gene BCOR
plot_grid(com_ch_plot$BCOR$cg02932805, com_mo_plot$BCOR$cg02932805) # for gene BCOR

t1 <- textGrob("Children", gp = gpar(fontsize = 16))
t2 <- textGrob("Mothers", gp = gpar(fontsize = 16))

p1 <- grid.arrange(arrangeGrob(com_ch_plot$BCOR$cg15627188, top = t1))
p2 <- grid.arrange(arrangeGrob(com_mo_plot$BCOR$cg15627188, top = t2))

bcor_p1 <- plot_grid(p1, p2, rel_heights=c(0.2, 2))
# bcor_p1 <- plot_grid(com_ch_plot$BCOR$cg15627188, com_mo_plot$BCOR$cg15627188,
#                      rel_heights=c(0.2, 2)) # for gene BCOR
y.grob <- textGrob(expression(atop("Methylation (%)", ~"cg15627188" ~ italic(" (BCOR)"))),
                   gp=gpar(fontsize=16), rot=90)

bcor_p1 <- grid.arrange(arrangeGrob(bcor_p1, left = y.grob))



plot_grid(com_ch_plot$PRDM8$cg18073471, com_mo_plot$PRDM8$cg18073471,labels = LETTERS[1:2]) # gene PRDM8
plot_grid(com_ch_plot$PRDM8$cg05059566, com_mo_plot$PRDM8$cg05059566,labels = LETTERS[1:2])
plot_grid(com_ch_plot$PRDM8$cg06307913, com_mo_plot$PRDM8$cg06307913,labels = LETTERS[1:2])


# ch_p.x <-  textGrob("Exposure to Genocide (Children)", gp=gpar(fontsize=16))
# mo_p.x <-  textGrob("Exposure to Genocide (Mothers)", gp=gpar(fontsize=16))
b.x = textGrob("Exposure to Genocide", gp=gpar(fontsize=16))

# ch_p <- grid.arrange(arrangeGrob(com_ch_plot$PRDM8$cg05059566,
#                                  top = ch_p.x))
# mo_p <- grid.arrange(arrangeGrob(com_mo_plot$PRDM8$cg05059566,
                                 # top = mo_p.x))

prdm8_p1 <- plot_grid(com_ch_plot$PRDM8$cg05059566,
                      com_mo_plot$PRDM8$cg05059566) # gene

prdm8_p1 <- grid.arrange(arrangeGrob(prdm8_p1, bottom = b.x))

y.grob <- textGrob(expression(atop("Methylation (%)", ~"cg05059566" ~ italic(" (PRDM8)"))),
                   gp=gpar(fontsize=16), rot=90)

prdm8_p1 <- grid.arrange(arrangeGrob(prdm8_p1, left = y.grob))

plot_grid(bcor_p1, prdm8_p1, nrow = 2, labels = 'AUTO')

png("E:/Rawanda EWAS/Plots/Analysis_Plots/bcor&prmd8_T-test1_up_dt.png", width = 550, height = 450)
plot_grid(bcor_p1, prdm8_p1, nrow = 2, labels = 'AUTO')
dev.off()


# VWDE -----------------
vwde_p1 <- grid.arrange(arrangeGrob(com_ch_plot$VWDE$cg03579179, top = t1))
vwde_p2 <- grid.arrange(arrangeGrob(com_mo_plot$VWDE$cg03579179, top = t2))

vwde_p <- plot_grid(vwde_p1, vwde_p2, rel_heights=c(0.2, 2))
# p1 <- plot_grid(com_ch_plot$VWDE$cg03579179, com_mo_plot$VWDE$cg03579179) # gene VWDE
y.grob <- textGrob(expression(atop("Methylation (%)", ~"cg03579179" ~ italic(" (VWDE)"))),
                   gp=gpar(fontsize=16), rot=90)

p1 <- grid.arrange(arrangeGrob(vwde_p, left = y.grob))


p2 <- plot_grid(com_ch_plot$VWDE$cg06484146, com_mo_plot$VWDE$cg06484146)
y.grob <- textGrob(expression(atop("Methylation (%)", ~"cg06484146" ~ italic(" (VWDE)"))),
                   gp=gpar(fontsize=16), rot=90)

p2 <- grid.arrange(arrangeGrob(p2, left = y.grob))

plot_grid(com_ch_plot$VWDE$cg20607287, com_mo_plot$VWDE$cg20607287)

#
# ch_p <- grid.arrange(arrangeGrob(com_ch_plot$VWDE$cg20607287,
#                                  bottom = ch_p.x))
# mo_p <- grid.arrange(arrangeGrob(com_mo_plot$VWDE$cg20607287,
#                                  bottom = mo_p.x))

vwde_p3 <- plot_grid(com_ch_plot$VWDE$cg20607287,
                     com_mo_plot$VWDE$cg20607287)
vwde_p3 <- grid.arrange(arrangeGrob(vwde_p3, bottom = b.x))
y.grob <- textGrob(expression(atop("Methylation (%)", ~"cg20607287" ~ italic(" (VWDE)"))),
                   gp=gpar(fontsize=16), rot=90)

vwde_p1 <- grid.arrange(arrangeGrob(vwde_p3, left = y.grob))

plot_grid(p1, p2, vwde_p1, nrow = 3)

png("E:/Rawanda EWAS/Plots/Analysis_Plots/vwde_T-test1_up_dt.png", width = 550, height = 670)
plot_grid(p1, p2, vwde_p1, nrow = 3, labels = 'AUTO')
dev.off()



# plot other two cpgs
plot_grid(com_ch_plot$VWDE$cg06484146, com_mo_plot$VWDE$cg06484146,labels = LETTERS[3:4])
plot_grid(com_ch_plot$VWDE$cg20607287, com_mo_plot$VWDE$cg20607287,labels = LETTERS[3:4])

ch_p.x <-  textGrob("Exposure to Genocide (Children)", gp=gpar(fontsize=16))
mo_p.x <-  textGrob("Exposure to Genocide (Mothers)", gp=gpar(fontsize=16))

ch_p <- grid.arrange(arrangeGrob(com_ch_plot$VWDE$cg06484146,
                                 bottom = ch_p.x))
mo_p <- grid.arrange(arrangeGrob(com_mo_plot$VWDE$cg03579179,
                                 bottom = mo_p.x))

vwde_p1 <- plot_grid(ch_p, mo_p,
                     labels = LETTERS[5:6])
y.grob <- textGrob(expression(atop("Methylation (%)", ~"cg03579179" ~ italic(" (VWDE)"))),
                   gp=gpar(fontsize=16), rot=90)

vwde_p1 <- grid.arrange(arrangeGrob(vwde_p1, left = y.grob))

plot_grid(bcor_p1, prdm8_p1, vwde_p1, nrow = 3)



# ---------------------- END -------------------------



# -------------------Over-representation analysis-------------------------------
# get gene annotation
get_annotation <- function(gene_ids){
  df <- bitr(gene_ids, fromType = "SYMBOL",
             toType = c("SYMBOL", "ENTREZID"),
             OrgDb = 'org.Hs.eg.db')

  }


# Pathway analysis of common genes
# function for over representation with background gene list
enrich_wd_background_gns <- function(gene_ls, bkg_gene_lis, ontology){
    ego <- enrichGO(gene          = gene_ls$ENTREZID,
                    OrgDb         = org.Hs.eg.db,
                    ont           = ontology,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    minGSSize = 10,
                    #universe = bkg_gene_lis$ENTREZID,
                    readable      = TRUE)

  }


bkg_genes <- c(myResults_ch_sig[, 1], myResults_mo_sig[, 1])
genes <- get_annotation(gene_ids = bkg_genes)
ch_genes <- get_annotation(gene_ids = myResults_ch_sig[, 1])
mo_genes <- get_annotation(gene_ids = myResults_mo_sig[, 1])

com_genes <- get_annotation(gene_ids = common_ch_mo$Row.names)



ontology_names <- c('CC', 'MF', 'BP')

# common in children and mothers
Over_rep_wd_background <- lapply(ontology_names, function(x){
  message("Performing overrepresentation analysis : ", x)
  over_rep <- enrich_wd_background_gns(gene_ls = com_genes, bkg_gene_lis = genes,
                                       ontology = x)
})

# Mothers
Over_rep_wd_background_mo <- lapply(ontology_names, function(x){
  message("Performing overrepresentation analysis : ", x)
  over_rep <- enrich_wd_background_gns(gene_ls = mo_genes, bkg_gene_lis = genes,
                                       ontology = x)
})



# Children
Over_rep_wd_background_ch <- lapply(ontology_names, function(x){
  message("Performing overrepresentation analysis : ", x)
  over_rep <- enrich_wd_background_gns(gene_ls = ch_genes, bkg_gene_lis = genes,
                                       ontology = x)
})


# plot
dotplot(Over_rep_wd_background[[1]], showCategory = 30,
        title = "Cellular Components")

heatplot(Over_rep_wd_background[[1]]) + ggplot2::ggtitle("Cellular Components")



# -----------------------------------------------
# For NR3C1
# Now check the direction of effect between exposed vs unexposed
# This analysis was done previously in 2019 were 14 cpgs are in the DMR region in mothers
# Now recently jan 2022, we thought to look at the direction

# cpgs
lea_edge_cpgs <- c("cg06968181", "cg06521673", "cg21702128", "cg10847032", "cg26720913",
          "cg16335926", "cg23430507", "cg24026230", "cg18146873", "cg06952416",
          "cg18068240", "cg08818984", "cg14558428", "cg18849621")

# pull out 14 cpgs from beta
lea_edge_beta <- beta_vals[which(rownames(beta_vals) %in% lea_edge_cpgs), ]
dim(lea_edge_beta)

# mothers
lea_edge_m <- lea_edge_beta[, which(colnames(lea_edge_beta) %in% colnames(beta_mothers))]
dim(lea_edge_m)
lea_edge_m <- merge(t(lea_edge_m), pheno_mothers, by = 0)


# children
lea_edge_c <- lea_edge_beta[, which(colnames(lea_edge_beta) %in% colnames(beta_children))]
dim(lea_edge_c)
lea_edge_c <- merge(t(lea_edge_c), pheno_children, by = 0)

get_t.test <- function(b_vals, g_name=NULL, s_name ){
  res_l <- list()
  cpgs <- colnames(b_vals)[which(grepl('^cg', colnames(b_vals)))]
  y <- "Exposure"
  for(i in 1:length(cpgs)){
    f <- as.formula(paste0(cpgs[i], '~', y))
    res <- t.test(f, data = b_vals)
    df <- data.frame('t' = res$statistic,'p.value' =  res$p.value,
                     "Mean No group" =  res$estimate[1],
                     "Mean Yes group" = res$estimate[2],
                     "conf level 95" = paste0(res$conf.int[1], ",", res$conf.int[2]))
    res_l[[length(res_l)+1]] <- df
  }
  names(res_l) <- cpgs
  do.call(rbind, res_l)
  }

m_ttest <- get_t.test(b_vals =  lea_edge_m)
colnames(m_ttest) <- paste0(colnames(m_ttest),".Mothers")


c_ttest <- get_t.test(b_vals =  lea_edge_c)
colnames(c_ttest) <- paste0(colnames(c_ttest),".Children")

comb <- merge(m_ttest, c_ttest, by = 0)
View(comb)
write.csv(comb, "E:/Rawanda EWAS/Data/NR3C1_leading_edge_cpgs_effects.csv", row.names = F)
