#.libPaths( "/users/yiwangz/lib" )
packages <- c("matrixStats","stats","cluster","ggplot2", "RVAideMemoire", 
              "e1071", "plot3D", "scatterplot3d", "reshape2", 
              "tidyr", "ggthemes","knockoff", "class", "bmrm", "purrr")
ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.r-project.org/")
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)


#####################################################################################
#####################################################################################
setwd("~/Google Drive/Yiwang_Zhou/Pilot_Projects/UKBB_TM_NLP/clustering_results/Manuscript_revision")
data <- read.csv("Whole_dataset_ukbb_pilot.csv", header = TRUE, 
                 na.strings = c("", " ", "NA", "NaN", "Inf", "."))
n.sub <- nrow(data)

set.seed(1234)
index_train <- sample(seq(1, n.sub, 1), floor(0.80*n.sub), replace = FALSE)
data.train <- data[index_train,]
data.test <- data[-index_train,]
n.train <- nrow(data.train)
n.test <- nrow(data.test)

#### biomarker and feature are both from the training dataset
biomarker <- data.train[,2:3298]
scaled.biomarker <- as.data.frame(lapply(biomarker, scale))
feature <- data.train[,3299:ncol(data.train)]

##################################################################################
######################### hierarchical clustering ################################
# Compute hierarchical clustering method in hclust (complete, average, single, ward.D)
# biomarker.hc.ward <- scaled.biomarker %>%
#   scale() %>%                    # Scale the data
#   dist(method = "euclidean") %>% # Compute dissimilarity matrix
#   hclust(method = "ward.D2")     # Compute hierachical clustering
# biomarker.hc.comp <- scaled.biomarker %>%
#   scale() %>%                    # Scale the data
#   dist(method = "euclidean") %>% # Compute dissimilarity matrix
#   hclust(method = "complete")     # Compute hierachical clustering
# biomarker.hc.sing <- scaled.biomarker %>%
#   scale() %>%                    # Scale the data
#   dist(method = "euclidean") %>% # Compute dissimilarity matrix
#   hclust(method = "single")     # Compute hierachical clustering
# biomarker.hc.aver <- scaled.biomarker %>%
#   scale() %>%                    # Scale the data
#   dist(method = "euclidean") %>% # Compute dissimilarity matrix
#   hclust(method = "average")     # Compute hierachical clustering
# 
# # visualize the hierarchical clustering
# pdf("hierarchical_clustering_ward.pdf")
# plot(biomarker.hc.ward)
# pdf("hierarchical_clustering_comp.pdf")
# plot(biomarker.hc.comp)
# pdf("hierarchical_clustering_sing.pdf")
# plot(biomarker.hc.sing)
# pdf("hierarchical_clustering_aver.pdf")
# plot(biomarker.hc.aver)
# dev.off()

m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# function to compute coefficient
ac <- function(x) {
  agnes(scaled.biomarker, method = x)$ac
}

capture.output(print(map_dbl(m, ac)), file = "agglomerative_coefficient.txt")
