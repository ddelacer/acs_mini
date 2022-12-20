################################################################################
# Broad Institute Mini project
# Compraisons consist of:
# Subset data of SNP data from individuals
# Will perform PCA and identify populations from individuals 
# 
# Created by David de la Cerda
# script name: mini_acs_PCA.R
# 
# 
# 
# input: VCFtools stats output, Plink PCA scores, Predicted ancestry PCA scores
#        using fraposa
# output: Histograms, PCA, predicted ancestry table, table with PCs per individual
# required software: R version 4.0.4
################################################################################
library(tidyverse)
library(egg)
###Raw vs Filtered vcf files
#raw inputs
#input variant quality
var_qual <- read_delim("mini_project.lqual", delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)
#missing variants
var_miss <- read_delim("mini_project.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
#input freq
var_freq <- read_delim("mini_project.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
#individ missing
ind_miss  <- read_delim("mini_project.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

#heterozygosity 
ind_het <- read_delim("mini_project.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
#change to filtered directory
#filtered inputs
var_qual_filt <- read_delim("mini_project.lqual", delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)
#missing variants
var_miss_filt <- read_delim("mini_project.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
#input freq
var_freq_filt <- read_delim("mini_project.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
#individ missing
ind_miss_filt  <- read_delim("mini_project.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

#heterozygosity 
ind_het_filt <- read_delim("mini_project.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)


#quality histograms
var_hist = ggplot(var_qual, aes(qual,fill="raw")) +
  geom_density(alpha = 0.3) +
  geom_density(aes(qual,fill="filtered"), alpha = 0.3,data=var_qual_filt) +
  scale_fill_manual(name = "VCF type", values = c(raw = "gray", filtered = "dodgerblue1")) +
  theme_minimal() +
  #xlim(0, 1e6)
  xlim(0, 1e5)

#missing
var_fmiss = ggplot(var_miss, aes(fmiss,fill="raw")) +
  geom_density(fill="gray",alpha = 0.3) +
  theme_minimal()

#summary(var_miss$fmiss)
#only zeros for filtered
#summary(var_miss_filt$fmiss)

#MAF
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
var_freq_filt$maf <- var_freq_filt %>% select(a1, a2) %>% apply(1, function(z) min(z))

var_maf = ggplot(var_freq, aes(maf,fill="raw")) +
  geom_density(alpha = 0.3) +
  geom_density(aes(maf,fill="filtered"), alpha = 0.5,data=var_freq_filt) +
  scale_fill_manual(name = "VCF type", values = c(raw = "gray", filtered = "dodgerblue1")) +
  theme_minimal()

#summary(var_freq$maf)
#summary(var_freq_filt$maf)

#missing per individuals, filtered has zero
ind_missing = ggplot(ind_miss, aes(fmiss)) +
  geom_density(fill = "gray", alpha = 0.3) +
  theme_minimal()
#summary(ind_miss_filt$fmiss)

#HZ
ind_hz = ggplot(ind_het, aes(f,fill="raw")) +
  geom_density(alpha = 0.3) +
  geom_density(aes(f,fill="filtered"), alpha = 0.3,data=ind_het_filt) +
  scale_fill_manual(name = "VCF type", values = c(raw = "gray", filtered = "dodgerblue1")) +
  theme_minimal()

#Stats comparision
#pdf("vcf_compare_stats.pdf")
ggarrange(var_hist,var_fmiss,var_maf,ind_missing,ind_hz)
#dev.off()

######PCA and prediction PCA
#read in PCA information from plink
pca_acs <- read_table2("acsminiproject.eigenvec", col_names = FALSE)
eigenval_acs <- scan("acsminiproject.eigenval")

pca_acs <- pca_acs[,-1]
# set names
names(pca_acs)[1] <- "ind"
names(pca_acs)[2:ncol(pca_acs)] <- paste0("PC", 1:(ncol(pca_acs)-1))

#import the individuals file provided
acs_ind = read_delim("acs_mini_project_labels.txt", delim = "\t",
                     col_names = c("ind", "ancestry"), skip = 1)
pca_merge <- as.tibble(data.frame(pca_acs, acs_ind))

#variance explained
pve_acs <- data.frame(PC = 1:20, pve = eigenval_acs/sum(eigenval_acs)*100)
#plot variance explained
ggplot(pve_acs, aes(PC, pve)) + geom_bar(stat = "identity") +
  ylab("Percentage variance explained") +
  theme_minimal()

#plot PCA with unknown individuals
pca_unknown1 = ggplot(pca_merge, aes(PC1, PC2, col = ancestry)) + geom_point(size = 2) +
  xlab(paste0("PC1 (", signif(pve_acs$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve_acs$pve[2], 3), "%)")) +
  coord_equal() + theme_minimal()
pca_unknown2 = ggplot(pca_merge, aes(PC1, PC3, col = ancestry)) + geom_point(size = 2) +
  xlab(paste0("PC1 (", signif(pve_acs$pve[1], 3), "%)")) + ylab(paste0("PC3 (", signif(pve_acs$pve[3], 3), "%)")) +
  coord_equal() + theme_minimal()
pca_unknown3 = ggplot(pca_merge, aes(PC1, PC4, col = ancestry)) + geom_point(size = 2) +
  xlab(paste0("PC1 (", signif(pve_acs$pve[1], 3), "%)")) + ylab(paste0("PC4 (", signif(pve_acs$pve[4], 3), "%)")) +
  coord_equal() + theme_minimal()

#output known and unknown ancestry IDs
unknown_ids = acs_ind[rowSums(is.na(acs_ind))>0,]
known_ids = subset(acs_ind,ancestry!="NA")
#need extra column for family
known_ids$fam = known_ids$ind
known_ids =known_ids[,c("ind","fam","ancestry")]
#export these ids for subsetting in bcftools
write.table(unknown_ids[,1],file="study_ids.txt",quote=F,sep="\t",row.names=F,col.names=F)
write.table(known_ids[,1],file="reference_ids.txt",quote=F,sep="\t",row.names=F,col.names=F)
write.table(known_ids,file="reference.popu",quote=F,sep="\t",row.names=F,col.names=F)
#import the predictions and see how good they are
pca_pred = read.table(file="study.popu",sep="\t",header = F)
#merge predictions with dataset
pca_merge_pred = merge(pca_merge,pca_pred,by.x="ind",by.y="V1",all.x=T)
pca_merge_pred$ancestry_all <- ifelse(is.na(pca_merge_pred$ancestry), pca_merge_pred$V3, pca_merge_pred$ancestry)


#new PC plot with predictions
pca_pred1 = ggplot(pca_merge_pred, aes(PC1, PC2, col = ancestry_all)) + geom_point(size = 2) +
  xlab(paste0("PC1 (", signif(pve_acs$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve_acs$pve[2], 3), "%)")) +
  coord_equal() + theme_minimal()

pca_pred2 = ggplot(pca_merge_pred, aes(PC1, PC3, col = ancestry_all)) + geom_point(size = 2) +
  xlab(paste0("PC1 (", signif(pve_acs$pve[1], 3), "%)")) + ylab(paste0("PC3 (", signif(pve_acs$pve[3], 3), "%)")) +
  coord_equal() + theme_minimal()

pca_pred3 = ggplot(pca_merge_pred, aes(PC1, PC4, col = ancestry_all)) + geom_point(size = 2) +
  xlab(paste0("PC1 (", signif(pve_acs$pve[1], 3), "%)")) + ylab(paste0("PC4 (", signif(pve_acs$pve[4], 3), "%)")) +
  coord_equal() + theme_minimal()

#plot together and export into pdf
pdf("pca_compare.pdf")
ggarrange(pca_unknown1, pca_pred1,
          pca_unknown2, pca_pred2,
          pca_unknown3, pca_pred3,nrow = 3,ncol=2)
dev.off()
#output PC scores
write.table(pca_acs,file="pruned_pca_scores.txt",quote=F,row.names=F,col.names=T)
#output the ancestry predictions
write.table(pca_merge_pred[,c("ind","ancestry","ancestry_all")],file="acs_mini_project_predicted_labels.txt",quote=F,row.names=F,col.names=T ) 


