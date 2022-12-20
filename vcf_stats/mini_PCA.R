################################################################################
# Broad Institute Mini project
# Compraisons consist of:
# Subset data of SNP data from individuals
# Will perform PCA and identify populations from individuals 
# 
# Created by David de la Cerda
# script name: mini_PCA.R
# 
# 
# 
# input: VCFtools output, subset VCF file
# output: Histograms, PCA
# required software: R version 4.0.4
################################################################################

library(tidyverse)
#input variant quality
var_qual <- read_delim("mini_project.lqual", delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)
#input depth
var_depth <- read_delim("mini_project.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
#missing variants
var_miss <- read_delim("mini_project.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
#input freq
var_freq <- read_delim("mini_project.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
#input individual
ind_depth <- read_delim("mini_project.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)
#individ missing
ind_miss  <- read_delim("mini_project.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

#heterozygosity 
ind_het <- read_delim("mini_project.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)

#quality histogram
var_hist = ggplot(var_qual, aes(qual)) +
   geom_density(fill = "dodgerblue1", color = "black", alpha = 0.3) +
   theme_minimal() +
   #xlim(0, 1e6)
  xlim(0, 1e5)

#depth histogram
#not working
summary(var_depth$mean_depth)
dep_hist = ggplot(var_depth, aes(mean_depth)) +
  geom_density(fill = "dodgerblue1", color = "black", alpha = 0.3) +
  theme_minimal()
#consider excluding outliers
#dep_hist + xlim(LOWERLIM, UPPERLIM)


var_fmiss = ggplot(var_miss, aes(fmiss)) +
  geom_density(fill = "dodgerblue1", color = "black", alpha = 0.3) +
  theme_minimal()

summary(var_miss$fmiss)

#MAF
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
var_maf = ggplot(var_freq, aes(maf)) +
  geom_density(fill = "dodgerblue1", color = "black", alpha = 0.3) +
  theme_minimal()
summary(var_freq$maf)

#missing per individuals
ggplot(ind_miss, aes(fmiss)) +
  geom_density(fill = "dodgerblue1", color = "black", alpha = 0.3) +
  theme_minimal()


#HZ
ggplot(ind_het, aes(f)) +
  geom_density(fill = "dodgerblue1", color = "black", alpha = 0.3) +
  theme_minimal()


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

#plot PCA
ggplot(pca_merge, aes(PC1, PC2, col = ancestry)) + geom_point(size = 2) +
xlab(paste0("PC1 (", signif(pve_acs$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve_acs$pve[2], 3), "%)")) +
coord_equal() + theme_minimal()


ggplot(pca_merge, aes(PC1, PC3, col = ancestry)) + geom_point(size = 2) +
  xlab(paste0("PC1 (", signif(pve_acs$pve[1], 3), "%)")) + ylab(paste0("PC3 (", signif(pve_acs$pve[3], 3), "%)")) +
  coord_equal() + theme_minimal()

#already have the "predictions" included, just need to define colors to match


ggplot(subset(pca_merge, ancestry != "na"), aes(PC1, PC2, col = ancestry)) + geom_point(size = 2) +
  geom_point(data=subset(pca_merge, ancestry=="na") , color="gray",size=2,shape=1) +
  
  xlab(paste0("PC1 (", signif(pve_acs$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve_acs$pve[2], 3), "%)")) +
  coord_equal() + theme_minimal()

pca_known = subset(pca_merge, ancestry != "na")

#adjusting code from Meyer Lab
#make function based on different ancestries provided
ancestry_labels = unique(pca_known$ancestry)

pred_ancestry = function(df,ac_labels){
  
  sub = subset(df,ancestry==ac_labels)
  pc1_mean = mean(sub$PC1)
  pc2_mean = mean(sub$PC2)
  #calculate euclidean distance
  sub$euclid_dist <- sqrt((sub$PC1 - pc1_mean)^2 +
                             (sub$PC2 - pc2_mean)^2)
  #calculate max euclidean distance
  max_euclid_dist <- max(sub$euclid_dist)
  
  #check against all samples
  df$euclid_dist <- sqrt((df$PC1 - pc1_mean)^2 +
                                 (df$PC2 - pc1_mean)^2)

  #df$assignment = ifelse(df$euclid_dist > (max_euclid_dist*1.5), paste("non-",ac_labels,sep=""),paste(ac_labels) )
  #binary instead of string
  df$assignment = ifelse(df$euclid_dist > (max_euclid_dist*1.5), 0,1 )
  estimate_df = df
  return(estimate_df)
}

#estimate_an_list_old = estimate_an_list

#list of data frames
estimate_an_list = list("list",length(ancestry_labels))
for(group in c(1:length(ancestry_labels))){
  estimate_an_list[[group]] = pred_ancestry(pca_merge,ancestry_labels[group])
}

names(estimate_an_list) = ancestry_labels
library(UpSetR)
#mutations <- read.csv( system.file("extdata", "mutations.csv", package = "UpSetR"), header=T, sep = ",")
#pull out unkown ancestry individuals from the list of estimates
unkown_estimate_an_list = lapply(estimate_an_list, subset, is.na(ancestry))
#multimerge function
multimerge <- function (mylist) {
  ## mimics a recursive merge or full outer join
  
  unames <- unique(unlist(lapply(mylist, rownames)))
  
  n <- length(unames)
  
  out <- lapply(mylist, function(df) {
    
    tmp <- matrix(nr = n, nc = ncol(df), dimnames = list(unames,colnames(df)))
    tmp[rownames(df), ] <- as.matrix(df)
    rm(df); gc()
    
    return(tmp)
  })
  
  stopifnot( all( sapply(out, function(x) identical(rownames(x), unames)) ) )
  
  bigout <- do.call(cbind, out)
  colnames(bigout) <- paste(rep(names(mylist), sapply(mylist, ncol)), unlist(sapply(mylist, colnames)), sep = "_")
  return(bigout)
}

unknown_merge = as.data.frame(multimerge(unkown_estimate_an_list))
#keep columns of interest
unknown_sets = unknown_merge[,c("eas_ind",grep("assignment",names(unknown_merge),value=T))]

unknown_sets_num = unknown_sets
#unknown_sets_num = as.data.frame(apply(unknown_sets[,2:length(names(unknown_sets))], 2, as.numeric)) 
unknown_sets_num = as.data.frame(apply(unknown_sets, 2, as.numeric))

unknown_sets_num[,2:length(names(unknown_sets))] <- sapply(unknown_sets_num[,2:length(names(unknown_sets))], as.numeric)


upset(unknown_sets_num,order.by = "freq",
      sets=c("afr_assignment","amr_assignment","eas_assignment","fin_assignment",
             "nfe_assignment","sas_assignment"),empty.intersections = "on")


upset(unknown_sets[,2:length(names(unknown_sets))] )

#output known and unknown ancestry IDs
unknown_ids = acs_ind[rowSums(is.na(acs_ind))>0,]
known_ids = subset(acs_ind,ancestry!="NA")
#need extra column for family
known_ids$fam = known_ids$ind
known_ids =known_ids[,c("ind","fam","ancestry")]

#export these ids for subsetting in other program
write.table(unknown_ids[,1],file="study_ids.txt",quote=F,sep="\t",row.names=F,col.names=F)
write.table(known_ids[,1],file="reference_ids.txt",quote=F,sep="\t",row.names=F,col.names=F)
write.table(known_ids,file="reference.popu",quote=F,sep="\t",row.names=F,col.names=F)

#import the predictions and see how good they are
pca_pred = read.table(file="study.popu",sep="\t",header = F)
#merge predictions with dataset
pca_merge_pred = merge(pca_merge,pca_pred,by.x="ind",by.y="V1",all.x=T)
pca_merge_pred$ancestry_all <- ifelse(is.na(pca_merge_pred$ancestry), pca_merge_pred$V3, pca_merge_pred$ancestry)


#new PC plot with predictions
ggplot(pca_merge_pred, aes(PC1, PC2, col = ancestry_all)) + geom_point(size = 2) +
  #geom_point(data=pca_merge_pred,aes(PC1, PC2, col = ancestry), size=2) +
  xlab(paste0("PC1 (", signif(pve_acs$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve_acs$pve[2], 3), "%)")) +
  coord_equal() + theme_minimal()

ggplot(pca_merge_pred, aes(PC1, PC3, col = ancestry_all)) + geom_point(size = 2) +
  #geom_point(data=pca_merge_pred,aes(PC1, PC2, col = ancestry), size=2) +
  xlab(paste0("PC1 (", signif(pve_acs$pve[1], 3), "%)")) + ylab(paste0("PC3 (", signif(pve_acs$pve[3], 3), "%)")) +
  coord_equal() + theme_minimal()

ggplot(pca_merge_pred, aes(PC1, PC4, col = ancestry_all)) + geom_point(size = 2) +
  #geom_point(data=pca_merge_pred,aes(PC1, PC2, col = ancestry), size=2) +
  xlab(paste0("PC1 (", signif(pve_acs$pve[1], 3), "%)")) + ylab(paste0("PC3 (", signif(pve_acs$pve[4], 3), "%)")) +
  coord_equal() + theme_minimal()


#functions for pulling ancestry comparisons
Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

Union <- function (x) {  
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}
Setdiff <- function (x, y) {
  # Remove the union of the y's from the common x's. 
  # x and y are lists of characters.
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}

vennlist = function(comparelist){
  combs <- 
    unlist(lapply(1:length(comparelist), 
                  function(j) combn(names(comparelist), j, simplify = FALSE)),
           recursive = FALSE)
  #make sure list is correct
  names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
  #get the names of them?
  #str(combs)
  #put the names in elements
  elements <- 
    lapply(combs, function(i) Setdiff(comparelist[i], comparelist[setdiff(names(comparelist), i)]))
  return(elements)
}



prac = subset(pca_known,ancestry == ancestry_labels[1])
pc1_mean = mean(prac$PC1)
pc2_mean = mean(prac$PC2)
prac$euclid_dist <- sqrt((prac$PC1 - pc1_mean)^2 +
                                   (prac$PC2 - pc2_mean)^2)
max_euclid_dist <- max(prac$euclid_dist)
#check against all samples
data_all = pca_merge
#calculate distances 
data_all$euclid_dist <- sqrt((data_all$PC1 - pc1_mean)^2 +
                                (data_all$PC2 - pc1_mean)^2)

data_all$assignment = ifelse(data_all$euclid_dist > (max_euclid_dist*1.5), paste("non-",ancestry_labels[1],sep=""),paste(ancestry_labels[1]) )

paste(ancestry_labels[1])
paste("non-",ancestry_labels[1],sep="")

df1$Category<-ifelse(df1$Frequency>25,"Good","Bad")


non_ancestry = dplyr::filter(data_all, .data$euclid_dist >
                               (max_euclid_dist * 1.5))

non_europeans <- dplyr::filter(data_name, .data$euclid_dist >
                                 (max_euclid_dist * europeanTh))

#all_european <- dplyr::filter(data_all, .data$Pop %in% refPopulation)
#euro_pc1_mean <- mean(all_european$PC1)
#euro_pc2_mean <- mean(all_european$PC2)

#all_european$euclid_dist <- sqrt((all_european$PC1 - euro_pc1_mean)^2 +
#                                   (all_european$PC2 - euro_pc2_mean)^2)

#max_euclid_dist <- max(all_european$euclid_dist)

## Find samples' distances to reference Europeans ####
data_name <- dplyr::filter(data_all, .data$Pop == name)
data_name$euclid_dist <- sqrt((data_name$PC1 - euro_pc1_mean)^2 +
                                (data_name$PC2 - euro_pc2_mean)^2)
non_europeans <- dplyr::filter(data_name, .data$euclid_dist >
                                 (max_euclid_dist * europeanTh))
fail_ancestry <- dplyr::select(non_europeans, .data$FID, .data$IID)
legend_rows <- round(nrow(colors)/legend_labels_per_row)







ggplot(df,aes(logFC,-log10(P.Value),color=-log10(P.Value)))+
  geom_point(shape=16,size=4, show.legend=F,color="#FF2600") + 
  #geom_point(data=subset(df,-log10(P.Value)<3) , color="gray",size=2) +
  geom_point(data=subset(df,logFC<=2 & logFC>0) , color="gray",size=4) +
  geom_point(data=subset(df,logFC>=-2 & logFC<0) , color="gray",size=4)



op <- par(mar=c(4,4,1,1), ps=10)
plot(pca$x[,pc], col=COLOR[iris.train$Species], cex=PCH[1], 
     xlab=paste0("PC ", pc[1], " (", expl.var[pc[1]], "%)"), 
     ylab=paste0("PC ", pc[2], " (", expl.var[pc[2]], "%)")
)
points(pred[,pc], col=COLOR[iris.valid$Species], pch=PCH[2])
legend("topright", legend=levels(iris$Species), fill = COLOR, border=COLOR)
legend("topleft", legend=c("training data", "validation data"), col=1, pch=PCH)
par(op)



pred <- predict(pca, newdata=iris.valid[,1:4])




# split data into 2 parts for pca training (75%) and prediction (25%)
set.seed(1)
samp <- sample(nrow(iris), nrow(iris)*0.75)
iris.train <- iris[samp,]
iris.valid <- iris[-samp,]

# conduct PCA on training dataset
pca <- prcomp(iris.train[,2:5], retx=TRUE, center=TRUE, scale=TRUE)
expl.var <- round(pca$sdev^2/sum(pca$sdev^2)*100) # percent explained variance

# prediction of PCs for validation dataset
pred <- predict(pca, newdata=iris.valid[,2:5])


?predict
?prcomp

#plink QC?
install.packages("plinkQC")
library(plinkQC)

?evaluate_check_ancestry

indir <- system.file("extdata", package="plinkQC")

