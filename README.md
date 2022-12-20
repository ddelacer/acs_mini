# Broad ACS mini project

##Created by David de la Cerda

Files here were provided from Dr. Grace Tiao and were used to generate a PCA plot and estimate the ancestry of unknown individuals in the dataset provided.

##Navigating the files

The primary files of interest are:

`mini_final.sh`

This executable was used to filter, generate PCs, and estimate ancestry

`mini_acs_PCA.R`

This script was used to visualize stats from the raw and filtered VCF files as well as visualize the PCA of individuals from unknown and predicted ancestry. It was also used to export tables related to ancestry classification and PC scores.

The directories are:

`vcf_stats`

This contains the descriptive statistics used to visualize the data provided.

`vcf_stats_filtered`

This contains the descriptive statistics used to visualize the filtered data once the raw data were visualized.

`deliverables`

This contains necessary files asked to be provided. The descriptions are:

`acs_mini_project_predicted_labels.txt`

This contains the information where the columns are original individual ids, original ancestry assignment, and predicted ancestry assignment.

`pca_compare.pdf`

This contains visualization of the distribution of PC values with the unknown and then the predicted classifications side-by-side. PC1 versus PC2, PC3, and PC4 are visualized to get a sense of separation.

`pruned_pca_scores.txt`

This file contains PC scores of each individual after pruning.

#Questions or concerns?
Please contact David de la Cerda at daviddlc09@gmail.com