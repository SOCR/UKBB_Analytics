#.libPaths( "/users/yiwangz/lib" )
packages <- c("matrixStats","stats","cluster","ggplot2", "RVAideMemoire", 
              "e1071", "plot3D", "scatterplot3d", "reshape2", 
              "tidyr", "ggthemes","knockoff", "class", "bmrm")
ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.r-project.org/")
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)


#####################################################################################
#####################################################################################
data <- read.csv("~/Google Drive/Yiwang_Zhou/Pilot_Projects/UKBB_TM_NLP/clustering_results/clustering_exploration/Whole_dataset_ukbb_pilot.csv", header = TRUE, 
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

#####################################################################################
#################### 1000 times repeat of k-means clustering ########################
kpp_init = function(dat, K, i) {
  set.seed(1234+i)
  
  x = as.matrix(dat)
  n = nrow(x)
  # Randomly choose a first center
  centers = matrix(NA, nrow=K, ncol=ncol(x))
  
  centers[1,] = as.matrix(x[sample(1:n, 1),])
  for (k in 2:K) {
    # Calculate dist^2 to closest center for each point
    dists = matrix(NA, nrow=n, ncol=k-1)
    for (j in 1:(k-1)) {
      temp = sweep(x, 2, centers[j,], '-')
      dists[,j] = rowSums(temp^2)
    }
    dists = rowMins(dists)
    # Draw next center with probability proportional to dist^2
    cumdists = cumsum(dists)
    prop = runif(1, min=0, max=cumdists[n])
    centers[k,] = as.matrix(x[min(which(cumdists > prop)),])
  }
  return(centers)
}

iter <- 1000
cls2 <- matrix(0,n.train,iter)
for(i in 1:iter){
  cluster_2 <- kmeans(scaled.biomarker, kpp_init(scaled.biomarker, 2, i), iter.max=100, algorithm='Lloyd')
  cls2[,i] <- cluster_2$cluster
}

cls2 <- as.data.frame(cls2)

write.csv(cls2, file = "kmeans_cluster2_1000.txt")







