#.libPaths( "/users/yiwangz/lib" )
packages <- c("matrixStats","stats","cluster","ggplot2", "RVAideMemoire", 
              "e1071", "plot3D", "scatterplot3d", "reshape2", 
              "tidyr", "ggthemes","knockoff", "class", "bmrm", "factoextra", "plotly",
              "plyr", "c060", "randomForest","caret", "partykit", "Rtsne", 
              "neuralnet", "meanShiftR")
ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.r-project.org/")
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)


#####################################################################################
#####################################################################################
################################ read in the entire dataset #########################
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

#####################################################################################
#####################################################################################
## check the distribution of the categorical variables in the training and testing dataset
index_sen_feel <- which(colnames(feature)=="X1950.0.0")
index_depress <- which(colnames(feature)=="X4598.2.0")
index_wor_feel <- which(colnames(feature)=="X1980.0.0")
index_miserable <- which(colnames(feature)=="X1930.2.0")

table(feature[,index_sen_feel])/n.train
table(feature[,index_depress])/n.train
table(feature[,index_wor_feel])/n.train
table(feature[,index_miserable])/n.train

table(data.test[,"X1950.0.0"])/n.test
table(data.test[,"X4598.2.0"])/n.test
table(data.test[,"X1980.0.0"])/n.test
table(data.test[,"X1930.2.0"])/n.test

#####################################################################################
#####################################################################################
## get the optimal cluster number
## 1) Silhouette value 2) within-cluster sum of squares 
## use the R scripts clustering_optimization_kmeans.R & clustering_optimization_hie.R
pdf("kmeans_tuning_silhouette_value.pdf")
fviz_nbclust(scaled.biomarker, kmeans, method = "silhouette",k.max=10)+
  geom_vline(xintercept=2,linetype = 2)+
  ylim(0,0.08)+
  scale_y_continuous(breaks=c(0.00,0.02,0.04,0.06,0.08))+
  ggtitle("")
dev.off()
labs(subtitle = "Silhouette method")+
  
  pdf("kmeans_tuning_wss.pdf")
fviz_nbclust(scaled.biomarker, kmeans, method = "wss") +
  geom_vline(xintercept=2,linetype = 3)+
  ggtitle("")
dev.off()

pdf("hierarchical_tuning_silhouette_value.pdf")
fviz_nbclust(scaled.biomarker, hcut, method = "silhouette")+
  geom_vline(xintercept=2,linetype = 2)+
  ylim(0,0.08)+
  scale_y_continuous(breaks=c(0.00,0.02,0.04,0.06,0.08))+
  ggtitle("")
dev.off()
labs(subtitle = "Silhouette method")+
  
  pdf("hierarchical_tuning_wss.pdf")
fviz_nbclust(scaled.biomarker, hcut, method = "wss") +
  geom_vline(xintercept=2,linetype = 3)+
  ggtitle("")
dev.off()
#####################################################################################
################################ k-means clustering #################################
################################## function define ##################################
############## obtain an optimal initialization for k-means clustering ##############
kpp_init = function(dat, K) {
  x = as.matrix(dat)
  n = nrow(x)
  # Randomly choose a first center
  centers = matrix(NA, nrow=K, ncol=ncol(x))
  set.seed(123)
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

cluster_2 <- kmeans(scaled.biomarker, kpp_init(scaled.biomarker, 2), iter.max=100, algorithm='Lloyd')
scaled.biomarker$cluster_2 <- cluster_2$cluster

#####################################################################################
############################### silhouette values ###################################
dis <- dist(scaled.biomarker[,-3298])
sil_2 <- silhouette(cluster_2$cluster, dis)

summary(sil_2)
# pdf("sil_kpp_init_kmeans_k2.pdf")
# plot(sil_2)
# dev.off()
#####################################################################################
#################### 1000 times repeat of k-means clustering ########################
## repeat kmeans clustering for 1000 times 
## use the R script kmeans_1000.R
iter <- 1000
cls2 <- matrix(0,n.train,iter)
for(i in 1:iter){
  cluster_2 <- kmeans(scaled.biomarker, kpp_init(scaled.biomarker, 2, i), iter.max=100, algorithm='Lloyd')
  cls2[,i] <- cluster_2$cluster
}

cls2 <- as.data.frame(cls2)

write.csv(cls2, file = "kmeans_cluster2_1000.txt")
#####################################################################################
########### agglomerative coefficient values for hierarchical clustering ############
## choose the linkage for hierarchical clustering
## use the R script hierarchical_agglomerative_coefficient.R
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# function to compute coefficient
ac <- function(x) {
  agnes(scaled.biomarker, method = x)$ac
}

capture.output(print(map_dbl(m, ac)), file = "agglomerative_coefficient.txt")
###################################################################################
###################################################################################
###################################################################################
## 8559 subjects have concordant clusters by kmeans and hierarchical clustering (2 clusters)
biomarker.hc.ward <- scaled.biomarker[,-3298] %>%
  scale() %>%                    # Scale the data
  dist(method = "euclidean") %>% # Compute dissimilarity matrix
  hclust(method = "ward.D2")     # Compute hierachical clustering

##save the hierarchical clustering results
hc_cluster_2 <- cutree(biomarker.hc.ward, 2)
kmeans_cluster_2 <- cluster_2$cluster

table(kmeans_cluster_2, hc_cluster_2)/n.train

kmeans_hc_data <- data.frame("kmeans"=kmeans_cluster_2, "hc"=hc_cluster_2)
write.csv(kmeans_hc_data, file = "kmeans_hc_cluster.csv")


#####################################################################################
###################### biomarker summary for 2 clusters #############################
#####################################################################################
cluster1 <- scaled.biomarker[which(scaled.biomarker$cluster_2 == 1),]
cluster2 <- scaled.biomarker[which(scaled.biomarker$cluster_2 == 2),]

cluster1 <- cluster1[,1:(dim(cluster1)[2]-1)]
n1 <- nrow(cluster1)  ##4242
cluster2 <- cluster2[,1:(dim(cluster2)[2]-1)]
n2 <- nrow(cluster2)  ##3689

mean1 <- apply(cluster1, 2, mean)
sd1 <- apply(cluster1, 2, sd)
median1 <- apply(cluster1, 2, median)
range1 <- apply(cluster1, 2, function(x) max(x)-min(x))
skew1 <- apply(cluster1, 2, function(x) skewness(x))
kurtosis1 <- apply(cluster1, 2, function(x) kurtosis(x))

mean2 <- apply(cluster2, 2, mean)
sd2 <- apply(cluster2, 2, sd)
median2 <- apply(cluster2, 2, median)
range2 <- apply(cluster2, 2, function(x) max(x)-min(x))
skew2 <- apply(cluster2, 2, function(x) skewness(x))
kurtosis2 <- apply(cluster2, 2, function(x) kurtosis(x))

## t-test
ttest <- function(i){
  model <- t.test(cluster1[,i], cluster2[,i], alternative = "two.sided", var.equal = FALSE)
  t <- model$statistic
  p <- model$p.value
  
  result <- list(t, p)
  names(result) <- c("t","p")
  return(result)
}

ttest_result <- lapply(1:ncol(biomarker), function(i) ttest(i))
t_test_value <- do.call(c, lapply(1:ncol(biomarker), function(i) ttest_result[[i]]$t))
t_test_p <- do.call(c, lapply(1:ncol(biomarker), function(i) ttest_result[[i]]$p))


## KS test
kstest <- function(i){
  model <- ks.test(cluster1[,i], cluster2[,i], alternative = "two.sided", exact = NULL)
  ks <- model$statistic
  p <- model$p.value
  
  result <- list(ks, p)
  names(result) <- c("ks","p")
  return(result)
}

kstest_result <- lapply(1:ncol(biomarker), function(i) kstest(i))
ks_test_value <- do.call(c, lapply(1:ncol(biomarker), function(i) kstest_result[[i]]$ks))
ks_test_p <- do.call(c, lapply(1:ncol(biomarker), function(i) kstest_result[[i]]$p))


## Wilcoxon Test 
wiltest <- function(i){
  model <- wilcox.test(cluster1[,i], cluster2[,i], alternative = "two.sided")
  wil <- model$statistic
  p <- model$p.value
  
  result <- list(wil, p)
  names(result) <- c("wil","p")
  return(result)
}

wiltest_result <- lapply(1:ncol(biomarker), function(i) wiltest(i))
wil_test_value <- do.call(c, lapply(1:ncol(biomarker), function(i) wiltest_result[[i]]$wil))
wil_test_p <- do.call(c, lapply(1:ncol(biomarker), function(i) wiltest_result[[i]]$p))


## unordered data.summary.cluster.2
data.summary.cluster.2 <- data.frame("mean_cluster1"=mean1,
                                     "sd_cluster1"=sd1,
                                     "median_cluster1"=median1,
                                     "range_cluster1"=range1,
                                     "skew_cluster1"=skew1,
                                     "kurtosis_cluster1"=kurtosis1,
                                     "mean_cluster2"=mean2,
                                     "sd_cluster2"=sd2,
                                     "median_cluster2"=median2,
                                     "range_cluster2"=range2,
                                     "skew_cluster2"=skew2,
                                     "kurtosis_cluster2"=kurtosis2,
                                     "t_value"=t_test_value,
                                     "t_p"=t_test_p,
                                     "ks_value"=ks_test_value,
                                     "ks_p"=ks_test_p,
                                     "wil_value"=wil_test_value,
                                     "wil_p"=wil_test_p)

#####################################################################################
########## Bonferroni correction of the threshold of the p-values  ##################
p.thresh <- 0.05/ncol(biomarker)
length(which(data.summary.cluster.2$t_p < p.thresh))
index.t <- which(data.summary.cluster.2$t_p < p.thresh)
sig.t <- numeric(nrow(data.summary.cluster.2))
sig.t[index.t] <- 1
data.summary.cluster.2$sig.t <- sig.t

length(which(data.summary.cluster.2$ks_p < p.thresh))
index.ks <- which(data.summary.cluster.2$ks_p < p.thresh)
sig.ks <- numeric(nrow(data.summary.cluster.2))
sig.ks[index.ks] <- 1
data.summary.cluster.2$sig.ks <- sig.ks

length(which(data.summary.cluster.2$wil_p < p.thresh))
index.wil <- which(data.summary.cluster.2$wil_p < p.thresh)
sig.wil <- numeric(nrow(data.summary.cluster.2))
sig.wil[index.wil] <- 1
data.summary.cluster.2$sig.wil <- sig.wil

length(which(data.summary.cluster.2$sig.t == 1 
             & data.summary.cluster.2$sig.ks == 1
             & data.summary.cluster.2$sig.wil == 1))

data.summary.cluster.2 <- data.summary.cluster.2[order(abs(data.summary.cluster.2$t_value), decreasing = TRUE),]
data.summary.cluster.2$t.order <- seq(1,nrow(data.summary.cluster.2),1)
data.summary.cluster.2 <- data.summary.cluster.2[order(data.summary.cluster.2$ks_value, decreasing = TRUE),]
data.summary.cluster.2$ks.order <- seq(1,nrow(data.summary.cluster.2),1)
data.summary.cluster.2 <- data.summary.cluster.2[order(data.summary.cluster.2$wil_p, decreasing = FALSE),]
data.summary.cluster.2$wil.order <- seq(1,nrow(data.summary.cluster.2),1)

##ordered data.summary.cluster.2
data.summary.cluster.2$sum.order <- data.summary.cluster.2$t.order+data.summary.cluster.2$ks.order+data.summary.cluster.2$wil.order
data.summary.cluster.2 <- data.summary.cluster.2[order(data.summary.cluster.2$sum.order),]

setwd("~/Google Drive/Yiwang_Zhou/Pilot_Projects/UKBB_TM_NLP/clustering_results/Manuscript_revision/result")
write.csv(data.summary.cluster.2, file = "data_summary_cluster_2.csv")

####################################################################################
####################################################################################
###################### density plot of the top 20 biomarkers #######################
setwd("~/Google Drive/Yiwang_Zhou/Pilot_Projects/UKBB_TM_NLP/clustering_results/Manuscript_revision/result")
data.summary.cluster.2 <- read.csv("data_summary_cluster_2.csv", header = TRUE)
top.selected <- as.character(data.summary.cluster.2[1:20,1])

scaled.biomarker$cluster_2 <- as.factor(scaled.biomarker$cluster_2)

bar.data <- scaled.biomarker[, top.selected]
colnames(bar.data) <- c(paste("rh_BA_exvivo_area_", "_rh_WhiteSurfArea_area", sep = "\n"),
                        paste("lh_BA_exvivo_area_", "_lh_WhiteSurfArea_area", sep = "\n"),
                        paste("rh_aparc_area_", "_lh_WhiteSurfArea_area", sep = "\n"),
                        paste("rh_aparc.a2009s_area_", "_rh_WhiteSurfArea_area", sep = "\n"),
                        paste("lh_aparc_area__lh_", "WhiteSurfArea_area", sep = "\n"),
                        paste("lh_aparc.a2009s_area_", "_lh_WhiteSurfArea_area", sep = "\n"),
                        "aseg__SupraTentorialVol",
                        "aseg__SupraTentorialVolNotVent",
                        "aseg__SupraTentorialVolNotVentVox",
                        "aseg__BrainSegVol",
                        "aseg__BrainSegVolNotVent",
                        "aseg__BrainSegVolNotVentSurf",
                        "aseg__CortexVol",
                        "aseg__rhCortexVol",
                        "aseg__MaskVol",
                        "aseg__lhCortexVol",
                        "aseg__TotalGrayVol",
                        paste("rh_aparc.DKTatlas_area_", "_rh_superiortemporal_area", sep = "\n"),
                        paste("rh_aparc.DKTatlas_area_", "_rh_superiorfrontal_area", sep = "\n"),
                        paste("lh_aparc.DKTatlas_area_", "_lh_superiortemporal_area", sep = "\n"))

bar.data$id <- 1:nrow(bar.data)
bar.data.m <- melt(bar.data, "id")
bar.data.m$cluster <- rep(scaled.biomarker$cluster_2, times=20)

png("Figure_4.png", height = 8, width = 15, units = "in", res = 300)
ggplot(data = bar.data.m, aes(x=value))+
  geom_density(aes(fill=cluster),alpha = 0.5,linetype="solid", kernel="epanechnikov", bw=0.3)+
  facet_wrap( ~ variable, ncol=5)+
  xlim(-5,5)+
  ylim(0,0.8)+
  scale_fill_manual(name=NULL,values=c("red2", "royalblue2"),labels=c("cluster1", "cluster2"))+
  theme(legend.position = "bottom",
        text = element_text(size=13),
        legend.text=element_text(size=14),
        panel.background = element_rect(fill = NA),
        axis.text = element_text(size=12),
        axis.title=element_text(size=14))
dev.off()


middle.selected <- as.character(data.summary.cluster.2[c(1640,1645,1649,
                                                         1651:1658,1661, 
                                                         1663,1664:1670),1])
bar.data <- scaled.biomarker[, middle.selected]
colnames(bar.data) <- c(paste("rh_aparc.a2009s_thickness_", "_rh_Pole_temporal_thickness", sep = "\n"),
                        paste("rh_aparc.DKTatlas_meancurv_", "_rh_rostralanteriorcingulate_meancurv", sep = "\n"),
                        paste("lh_BA_exvivo_thicknessstd_", "_lh_BA4a_exvivo_thicknessstd", sep = "\n"),
                        paste("rh_aparc.a2009s_thicknessstd_", "_rh_S_oc.temp_med.Lingual_thicknessstd", sep = "\n"),
                        paste("lh_aparc.DKTatlas_gauscurv_", "_lh_caudalanteriorcingulate_gauscurv", sep = "\n"),
                        paste("rh_aparc.pial_meancurv_", "_rh_parsopercularis_meancurv", sep = "\n"),
                        paste("rh_aparc_meancurv_", "_rh_rostralmiddlefrontal_meancurv", sep = "\n"),
                        paste("rh_aparc.a2009s_thicknessstd_", "_rh_G_cingul.Post.dorsal_thicknessstd", sep = "\n"),
                        paste("rh_aparc.a2009s_meancurv_", "_rh_G_temp_sup.Plan_polar_meancurv", sep = "\n"),
                        paste("lh_aparc_meancurv_", "_lh_supramarginal_meancurv", sep = "\n"),
                        paste("rh_aparc_meancurv_", "_rh_supramarginal_meancurv", sep = "\n"),
                        paste("lh_aparc.pial_meancurv_", "_lh_rostralanteriorcingulate_meancurv", sep = "\n"),
                        paste("lh_aparc.a2009s_meancurv_", "_lh_S_temporal_inf_meancurv", sep = "\n"),
                        paste("rh_aparc.DKTatlas_meancurv_", "_rh_posteriorcingulate_meancurv", sep = "\n"),
                        paste("rh_aparc.a2009s_foldind_", "_rh_G_oc.temp_med.Parahip_foldind", sep = "\n"),
                        paste("rh_aparc_thickness_", "_rh_parahippocampal_thickness", sep = "\n"),
                        paste("lh_aparc.pial_meancurv_", "_lh_medialorbitofrontal_meancurv", sep = "\n"),
                        paste("rh_aparc.a2009s_meancurv_", "_rh_G_rectus_meancurv", sep = "\n"),
                        paste("rh_aparc_gauscurv_", "_rh_caudalanteriorcingulate_gauscurv", sep = "\n"),
                        paste("rh_BA_exvivo_thicknessstd_", "_rh_BA4a_exvivo_thicknessstd", sep = "\n"))

bar.data$id <- 1:nrow(bar.data)
bar.data.m <- melt(bar.data, "id")
bar.data.m$cluster <- rep(scaled.biomarker$cluster_2, times=20)

png("Figure_S4a.png", height = 8, width = 15, units = "in", res = 300)
ggplot(data = bar.data.m, aes(x=value))+
  geom_density(aes(fill=cluster),alpha = 0.5,linetype="solid", kernel="epanechnikov", bw=0.3)+
  facet_wrap( ~ variable, ncol=5)+
  xlim(-5,5)+
  ylim(0,0.8)+
  scale_fill_manual(name=NULL,values=c("red2", "royalblue2"),labels=c("cluster1", "cluster2"))+
  theme(legend.position = "bottom",
        text = element_text(size=11),
        legend.text=element_text(size=14),
        panel.background = element_rect(fill = NA),
        axis.text = element_text(size=12),
        axis.title=element_text(size=14))
dev.off()



bottom.selected <- as.character(data.summary.cluster.2[c(3275:3278,3279,3282,
                                                         3283,3285:3297),1])
bar.data <- scaled.biomarker[, bottom.selected]
colnames(bar.data) <- c(paste("lh_aparc.a2009s_gauscurv_", "_lh_G.S_occipital_inf_gauscurv", sep = "\n"),
                        paste("rh_aparc.a2009s_thickness_", "_rh_S_occipital_ant_thickness", sep = "\n"),
                        paste("lh_aparc.a2009s_meancurv_", "_lh_G_temporal_middle_meancurv", sep = "\n"),
                        paste("rh.w.g.pct.mean_", "_transversetemporal", sep = "\n"),
                        paste("rh_aparc_thickness_", "_rh_precuneus_thickness", sep = "\n"),
                        paste("lh_aparc.DKTatlas_thicknessstd_", "_lh_parsorbitalis_thicknessstd", sep = "\n"),
                        paste("rh_aparc.pial_meancurv_", "_rh_bankssts_meancurv", sep = "\n"),
                        paste("rh_aparc.DKTatlas_thickness_", "_rh_precuneus_thickness", sep = "\n"),
                        paste("rh_aparc.DKTatlas_thickness_", "_rh_transversetemporal_thickness", sep = "\n"),
                        paste("rh_BA_exvivo_thickness_", "_rh_V2_exvivo_thickness", sep = "\n"),
                        paste("lh_aparc_meancurv_", "_lh_pericalcarine_meancurv", sep = "\n"),
                        paste("rh_BA_exvivo_thickness_", "_rh_MT_exvivo_thickness", sep = "\n"),
                        paste("rh_aparc.a2009s_thickness_", "_rh_S_temporal_sup_thickness", sep = "\n"),
                        paste("rh_BA_exvivo.thresh_meancurv_", "_rh_entorhinal_exvivo_meancurv", sep = "\n"),
                        paste("lh_BA_exvivo_thickness_", "_lh_V2_exvivo_thickness", sep = "\n"),
                        paste("lh_BA_exvivo_meancurv_", "_lh_V2_exvivo_meancurv", sep = "\n"),
                        paste("rh_aparc_thicknessstd_", "_rh_caudalmiddlefrontal_thicknessstd", sep = "\n"),
                        paste("rh_aparc_thickness_", "_rh_transversetemporal_thickness", sep = "\n"),
                        paste("lh_aparc.a2009s_thicknessstd_", "_lh_Lat_Fis.ant.Vertical_thicknessstd", sep = "\n"),
                        paste("rh_aparc.a2009s_thickness_", "_rh_G_temp_sup.G_T_transv_thickness", sep = "\n"))

bar.data$id <- 1:nrow(bar.data)
bar.data.m <- melt(bar.data, "id")
bar.data.m$cluster <- rep(scaled.biomarker$cluster_2, times=20)

png("Figure_S4b.png", height = 8, width = 15, units = "in", res = 300)
ggplot(data = bar.data.m, aes(x=value))+
  geom_density(aes(fill=cluster),alpha = 0.5,linetype="solid", kernel="epanechnikov", bw=0.3)+
  facet_wrap( ~ variable, ncol=5)+
  xlim(-5,5)+
  ylim(0,0.8)+
  scale_fill_manual(name=NULL,values=c("red2", "royalblue2"),labels=c("cluster1", "cluster2"))+
  theme(legend.position = "bottom",
        text = element_text(size=11),
        legend.text=element_text(size=14),
        panel.background = element_rect(fill = NA),
        axis.text = element_text(size=12),
        axis.title=element_text(size=14))
dev.off()

################################################################################
#################################################################################
## summary statistics 
data_summary <- biomarker[,top.selected]
data_summary_1 <- data_summary[cluster_2$cluster==1,]
data_summary_2 <- data_summary[cluster_2$cluster==2,]

summary_table <- data.frame("names"=colnames(data_summary),
                            "mean1"=apply(data_summary_1, 2, mean),
                            "median1"=apply(data_summary_1, 2, median),
                            "sd1"=apply(data_summary_1, 2, sd),
                            "mean2"=apply(data_summary_2, 2, mean),
                            "median2"=apply(data_summary_2, 2, median),
                            "sd2"=apply(data_summary_2, 2, sd))

rownames(summary_table) <- NULL
summary_table[,2:7] <- sapply(2:7, function(i) round(summary_table[,i], digits = 0))
write.csv(summary_table, file = "Table1.csv")

############################## supervised clustering #############################
##################################################################################
######################### K-nearest neighbor clustering ##########################
## we need to use cross validation to validate the prediction of the cluster labels
set.seed(1234)
knn_cv_index <- balanced.cv.fold(scaled.biomarker$cluster_2, num.cv = 5)

knn_pred <- numeric(5)
for(i in 1:5){
  knn_train <- scaled.biomarker[which(knn_cv_index!=i),]
  knn_valid <- scaled.biomarker[which(knn_cv_index==i),]
  
  model <- knn(knn_train[,-3298], knn_valid[,-3298], knn_train$cluster_2, k = 10, prob = FALSE, use.all = TRUE)
  knn_pred[i] <- sum(model==knn_valid$cluster_2)/nrow(knn_valid)
}
mean(knn_pred)  ##0.937 (a very high consistency)
##################################################################################
######################## artificial neural networks  #############################
## we need to use cross validation to validate the prediction of the cluster labels
set.seed(1234)
ann_cv_index <- balanced.cv.fold(scaled.biomarker$cluster_2, num.cv = 5)

pred <- function(nn, dat) {
  yhat = compute(nn, dat)$net.result
  yhat = apply(yhat, 1, which.max)-1
  return(yhat)
}

ann_pred <- numeric(5)
for(i in 1:5){
  ann_train_x <- scaled.biomarker[which(ann_cv_index!=i),-3298]
  ann_train_y <- scaled.biomarker[which(ann_cv_index!=i),3298]
  ann_train_y_ind <- model.matrix(~factor(ann_train_y)-1)
  colnames(ann_train_y_ind) = c("cluster1","cluster2")
  ann_train = cbind(ann_train_x, ann_train_y_ind)
  ann_valid <- scaled.biomarker[which(ann_cv_index==i),]
  
  ann_char <- "cluster1+cluster2~lh_aparc_area__lh_bankssts_area"
  for(j in 2:3297){
    ann_char <- paste(ann_char, colnames(ann_train)[j], sep = "+")
  }
  ann_formula <- as.formula(ann_char)
  
  model <- neuralnet(ann_formula,data=ann_train,linear.output=FALSE,lifesign='full')
  
  ann_pred[i] <- mean(pred(model, ann_valid[,-3298]) == as.factor(ann_valid[,3298]-1))
  #table(pred(model, ann_valid[,-11]), as.factor(ann_valid[,11]-1))
}
mean(ann_pred) ##0.973

##################################################################################
############################### mixture resolving  ###############################
## mixture resolving is a model based clustering method
## mixture resolving is an unsupervised clustering method

##################################################################################
########################### mode seeking algorithms ##############################
# set.seed(1234)
# 
# msc_train <- as.matrix(scaled.biomarker)
# 
# h <- rep(0.5,ncol(msc_train)-1)
# model <- meanShift(msc_train[,-3298],msc_train[,-3298],alpha=0,
#                    bandwidth=h,iterations=1)
# msc_cluster <- model$assignment
# 
# length(which(msc_cluster == msc_train[,3298]))/nrow(msc_train)
###################################################################################
###################################################################################
#################### clustering consistency 1000 kmeans ###########################
temp.cluster2 <- read.csv("kmeans_cluster2_1000.txt")
temp.cluster2 <- temp.cluster2[,-1]

for(i in 1:1000){
  if(temp.cluster2[1,i] == 1){
    next
  }
  else{
    for(j in 1:n.train){
      if(temp.cluster2[j,i] == 2){
        temp.cluster2[j,i] <- 1
      }
      else{
        temp.cluster2[j,i] <- 2
      }
    }
  }
}

mis.cluster2 <- numeric(n.train)
for(i in 1:n.train){
  mis.cluster2[i] <- min(length(which(temp.cluster2[i,] == 1)), length(which(temp.cluster2[i,] == 2)))
}
mis.cluster2.table <- as.data.frame(table(mis.cluster2))  ## 57 is the cutoff (we can set it as 60)

c1 <- numeric(n.train)
for(i in 1:n.train){
  c1[i] <- length(which(temp.cluster2[i,] == 1))
}
true.cluster2 <- rep(1, n.train)
for(i in 1:n.train){
  if(c1[i] < 500){
    true.cluster2[i] <- 2
  }
}

index.1 <- which(mis.cluster2 <= 60 & true.cluster2 == 1)  # cluster 1
index.2 <- which(mis.cluster2 <= 60 & true.cluster2 == 2)  # cluster 2
index.3 <- which(mis.cluster2 > 60)

label <- numeric(n.train)
label[index.1] <- 1
label[index.2] <- 2
label[index.3] <- 3
label <- as.factor(label)

dist_eu <- dist(scaled.biomarker[,-3298], method = "euclidean")
save(dist_eu, file = "dist_eu.Rdata")
mds_eu <- cmdscale(dist_eu, eig=FALSE, k=2)

colnames(mds_eu) <- c("Coordinate1", "Coordinate2")
mds_eu <- as.data.frame(mds_eu)
mds_eu$label <- label[1:n.train]

png("Figure_2a.png", width = 16, height = 8, units = 'in', res=300)
ggplot(mds_eu, aes(x=Coordinate1, y=Coordinate2, color=label)) +
  geom_point(size = 2)+
  scale_color_manual(name=NULL,
                     breaks = c("1", "2", "3"),
                     values=c("red2", "royalblue2", "springgreen3"),
                     labels=c("MCR<=0.06 Cluster=1", "MCR<=0.06 Cluster=2", "MCR>0.06"))+
  theme_bw()+
  theme(text = element_text(size=20), 
        legend.text=element_text(size=18), legend.position = "bottom")
dev.off()

########## calculating kmeans clustering consistency
Niter <- 999
cluster0_1 <- which(temp.cluster2[,1]==1)
cluster0_2 <- which(temp.cluster2[,1]==2)

stat <- c()

for (i in 1:999){
  
  c1_new <- table(temp.cluster2[cluster0_1,i+1]) #clustering result among original c1
  c2_new <- table(temp.cluster2[cluster0_2,i+1]) #clustering result among original c2
  
  stat <-cbind(stat,c(max(prop.table(c1_new)),max(prop.table(c2_new))))
}

robust_result <- as.data.frame(cbind(rowMeans(stat),apply(stat,1,sd)))
names(robust_result) <- c("Mean_Proportions","SD")
write.csv(robust_result,"clustering(robustness).csv",row.names = F)

##################################################################################
##################################################################################
###################################### PCA #######################################
set.seed(1234)
pca_biomarker <- prcomp(as.matrix(biomarker), scale. = TRUE)
summary(pca_biomarker)

pca_x <- as.data.frame(pca_biomarker$x)

eigen <- get_eigenvalue(pca_biomarker)
eigen

qualit_vars <- as.factor(cluster_2$cluster)

png("Figure_2b.png", height = 6, width = 8, units = "in", res = 300)
fviz_pca_ind(pca_biomarker, axes = c(1, 2), label = "none", 
             habillage = qualit_vars, addEllipses = TRUE, 
             title = "")+
  scale_color_manual(values=c("red2", "royalblue2","green"))+
  theme(legend.position = "bottom",
        panel.background=element_rect(fill = NA, colour="gray25"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey90"),
        text = element_text(size=20),
        legend.text=element_text(size=18))+
  xlab("PC1 (11.5%)")+
  ylab("PC2 (9.1%)")
dev.off()



## 3D plot of PCA
m.data <- as.data.frame(pca_biomarker$x)
m.data <- m.data[,1:3]
colnames(m.data) <- c("PC1", "PC2", "PC3")
m.data$label <- cluster_2$cluster
m.data$label <- as.factor(m.data$label)
p <- plot_ly(m.data, x = ~PC1, y = ~PC2, z = ~PC3, 
             color = ~label, colors = c("red2", "royalblue2")) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = "PC1"),
                      yaxis = list(title = "PC2"),
                      zaxis = list(title = "PC3")))
p

###################################################################################
###################################################################################
################################# categorical test ################################
numeric.feature <- as.data.frame(matrix(0,nrow = n.train, ncol = 0))

## all the numeric (including integer) data
for(i in 1:ncol(feature)){
  if(class(feature[,i]) == "numeric" | class(feature[,i]) == "integer"){
    numeric.feature <- cbind(numeric.feature, feature[,i])
    colnames(numeric.feature)[ncol(numeric.feature)] <- colnames(feature)[i]
  }
}

numeric.feature <- numeric.feature[!(apply(numeric.feature, 2, function(x) all(is.na(x))))]
numeric.feature <- numeric.feature[!(apply(numeric.feature, 2, function(x) length(which(!is.na(x)))==1))]

###################################################################################
integer.feature <- as.data.frame(matrix(0,nrow = n.train, ncol = 0))

for(i in 1:ncol(numeric.feature)){
  if(class(numeric.feature[,i]) == "integer"){
    integer.feature <- cbind(integer.feature, numeric.feature[,i])
    colnames(integer.feature)[ncol(integer.feature)] <- colnames(numeric.feature)[i]
  }
}

integer.feature[integer.feature == "-3"] <- NA
integer.feature[integer.feature == "-1"] <- NA
integer.feature[integer.feature == "-6"] <- NA
integer.feature <- as.data.frame(lapply(integer.feature, factor))
integer.feature$cluster <- scaled.biomarker$cluster_2


chi_value <- rep(NA,ncol(integer.feature))
chi_p <- rep(NA,ncol(integer.feature))
fisher_p <- rep(NA,ncol(integer.feature))

for(i in 1:ncol(integer.feature)){
  if(length(levels(integer.feature[,i])) == 2){
    tt <- as.matrix(table(integer.feature[,i], integer.feature$cluster))
    tryCatch(
      if(tt[1,1]>=5 & tt[1,2]>=5 & tt[2,1]>=5 & tt[2,2]>=5){
        temp.test <- chisq.test(integer.feature[,i], integer.feature$cluster)
        chi_value[i] <- temp.test$statistic
        chi_p[i] <- temp.test$p.value
      }, error=function(e) NULL
    )
  }
  else if(length(levels(integer.feature[,i])) > 2 & length(levels(integer.feature[,i])) <= 4){
    tt <- as.matrix(table(integer.feature[,i], integer.feature$cluster))
    temp.test <- fisher.test(tt, alternative = "two.sided", workspace = 800000,hybrid = TRUE)
    fisher_p[i] <- temp.test$p.value
  }
}


chi_test <- data.frame("names"=colnames(integer.feature), "chi_value"=chi_value, "chi_p"=chi_p)
chi_test <- chi_test[complete.cases(chi_test),]
chi_test <- chi_test[order(chi_test$chi_p, decreasing = FALSE),]

p.thresh <- 0.05/nrow(chi_test)
index.chi <- which(chi_test$chi_p < p.thresh)
selected.chi <- chi_test[index.chi,]


fisher_test <- data.frame("names"=colnames(integer.feature), "fisher_p"=fisher_p)
fisher_test <- fisher_test[complete.cases(fisher_test),]
fisher_test <- fisher_test[order(fisher_test$fisher_p, decreasing = FALSE),]

p.thresh <- 0.05/nrow(fisher_test)
index.fisher <- which(fisher_test$fisher_p < p.thresh)
selected.fisher <- fisher_test[index.fisher,]

write.csv(selected.chi, "selected_categorical_feature_chi.csv")
write.csv(selected.fisher, "selected_categorical_feature_fisher.csv")

#####
table(integer.feature[,"X20117.0.0"], integer.feature$cluster)

table(integer.feature[integer.feature$cluster==1, "X20117.0.0"])/length(which(integer.feature$cluster==1))
table(integer.feature[integer.feature$cluster==2, "X20117.0.0"])/length(which(integer.feature$cluster==2))

###################################################################################
###################################################################################
###################################################################################
selected.chi <- read.csv("selected_categorical_feature_chi.csv", header=TRUE)
selected.chi$names <- as.character(selected.chi$names)
selected.chi <- selected.chi[which(selected.chi$names %in% c("X31.0.0", "X1950.0.0", "X1980.0.0",
                                                             "X2040.2.0", "X4598.2.0")),]


names <- as.character(selected.chi$names)
description <- as.character(selected.chi$description)

index <- which(colnames(integer.feature) %in% names)
int.feature <- integer.feature[,index]
int.feature$cluster <- integer.feature$cluster

int.feature$X31.0.0 <- revalue(int.feature$X31.0.0, c("0"="female", "1"="male"))
int.feature$X1950.0.0 <- revalue(int.feature$X1950.0.0, c("1"="yes", "0"="no"))
int.feature$X1980.0.0 <- revalue(int.feature$X1980.0.0, c("1"="yes", "0"="no"))
int.feature$X2040.2.0 <- revalue(int.feature$X2040.2.0, c("1"="yes", "0"="no"))
int.feature$X4598.2.0 <- revalue(int.feature$X4598.2.0, c("1"="yes", "0"="no"))

seq.formula <- list(nrow(selected.chi)) 
for(i in 1:nrow(selected.chi)){
  seq.formula[[i]] <-  as.formula(paste("~ cluster", names[i], sep = "+"))
}

i<-5
mosaicplot(seq.formula[[i]], data = int.feature, shade = TRUE,
           main = "", xlab="cluster", ylab=description[i], cex=1)
###################################################################################
###################################################################################
###################################################################################
selected.fisher <- read.csv("selected_categorical_feature_fisher.csv", header=TRUE)
selected.fisher$names <- as.character(selected.fisher$names)
selected.fisher <- selected.fisher[which(selected.fisher$names=="X1200.0.0"),]

names <- as.character(selected.fisher$names)
description <- as.character(selected.fisher$description)

index <- which(colnames(integer.feature) %in% names)
int.feature <- integer.feature[,index]
int.feature <- data.frame("X1200.0.0"=int.feature, "cluster"=integer.feature$cluster)  

int.feature$X1200.0.0 <- revalue(int.feature$X1200.0.0, c("1"="Never/rarely", "2"="Sometimes", "3"="Usually"))

seq.formula <-  as.formula(paste("~ cluster", names, sep = "+"))

mosaicplot(seq.formula, data = int.feature, shade = TRUE,
           main = "", xlab="cluster", ylab=description, cex=1)


################################################################################
################################################################################
############################## random forest ###################################
biomarker.20 <- read.csv("data_summary_cluster_2.csv", header=TRUE)
selected.biomarker <- biomarker.20[1:20,1]
selected.biomarker <- as.character(selected.biomarker)

#################################################################################
##### Sensitivity / hurt feelings ###############################################
feature.20 <- c("X31.0.0","X1980.0.0","X2040.2.0","X2030.0.0",
                "X1618.2.0","X2090.0.0","X2000.0.0","X1930.0.0",
                "X1210.0.0","X1970.2.0","X4653.2.0","X4598.2.0",
                "X4631.2.0","X1200.0.0","X1170.0.0",
                "X1190.2.0","X2080.0.0","X20117.0.0","X2877.0.0",
                "X1249.2.0","X2050.0.0")

num.miss <- numeric(length(feature.20))
for(i in 1:length(feature.20)){
  num.miss[i] <- length(which(is.na(feature[,feature.20[i]])))
}
index <- which(num.miss > 1000)
feature.20 <- feature.20[-index]


rf.data <- cbind(scaled.biomarker[,selected.biomarker],feature[,feature.20])
rf.data[rf.data==-3] <- NA
rf.data[rf.data==-1] <- NA
rf.data[rf.data==-6] <- NA
biomarker.names <- c("rh_BA_exvivo_area__rh_WhiteSurfArea_area",
                     "lh_BA_exvivo_area__lh_WhiteSurfArea_area",
                     "rh_aparc_area__rh_WhiteSurfArea_area",
                     "rh_aparc.a2009s_area__rh_WhiteSurfArea_area",
                     "lh_aparc_area__lh_WhiteSurfArea_area",
                     "lh_aparc.a2009s_area__lh_WhiteSurfArea_area",
                     "aseg__SupraTentorialVol",
                     "aseg__SupraTentorialVolNotVent",
                     "aseg__SupraTentorialVolNotVentVox",
                     "aseg__BrainSegVol",
                     "aseg__BrainSegVolNotVent",
                     "aseg__BrainSegVolNotVentSurf",
                     "aseg__CortexVol",
                     "aseg__rhCortexVol",
                     "aseg__MaskVol",
                     "aseg__lhCortexVol",
                     "aseg__TotalGrayVol",
                     "rh_aparc.DKTatlas_area__rh_superiortemporal_area",
                     "rh_aparc.DKTatlas_area__rh_superiorfrontal_area",
                     "lh_aparc.DKTatlas_area__lh_superiortemporal_area")  

colnames(rf.data) <- c(biomarker.names, 
                       "sex",
                       "Worrier_anxious_feelings",
                       "Risk_taking",
                       "Guilty_feelings",
                       "Alcohol_usually_taken_with_meals",
                       "Seen_doctor_for_nerves_anxiety_tension_depression",
                       "Worry_too_long_after_embarrassment",
                       "Miserableness",
                       "Snoring",
                       "Nervous_feelings",
                       "Ever_highly_irritable_argumentative_for_2_days",
                       "Ever_depressed_for_a_whole_week",
                       "Ever_unenthusiastic_disinterested_for_a_whole_week",
                       "Sleeplessness_insomnia",
                       "Getting_up_in_morning",
                       "Nap_during_day",
                       "Frequency_of_tiredness_lethargy_in_last_2_weeks",
                       "Alcohol_drinker_status",
                       "Past_tobacco_smoking",
                       "Frequency_of_depressed_mood_in_last_2_weeks")  
  

rf.data$label <- feature[,"X1950.0.0"]
rf.data$label[rf.data$label=="-3"] <- NA
rf.data$label[rf.data$label=="-1"] <- NA
rf.data$label <- as.factor(rf.data$label)

set.seed(1234)
cross.validation <- function(rf.data, nfolds, iter){
  
  pred.accuracy <- numeric(iter)
  pred.accuracy.upper <- numeric(iter)
  pred.accuracy.lower <- numeric(iter)
  mcnemar.p <- numeric(iter)
  sensitivity <- numeric(iter)
  specificity <- numeric(iter)
  var.importance <- NULL
  
  for(i in 1:iter){
    foldid <- balancedFolds(rf.data$label, cross.outer=nfolds)
    
    pred.accuracy.temp <- numeric(nfolds)
    pred.accuracy.upper.temp <- numeric(nfolds)
    pred.accuracy.lower.temp <- numeric(nfolds)
    mcnemar.p.temp <- numeric(nfolds)
    sensitivity.temp <- numeric(nfolds)
    specificity.temp <- numeric(nfolds)
    
    var.imp.temp <- NULL
    
    for(j in 1:nfolds){
      index <- which(foldid == j)
      rf.train <- rf.data[-index,]
      rf.test <- rf.data[index,]
      
      fit <- randomForest(label ~ ., data = rf.train, ntree=500, importance=T,
                          na.action=na.omit)
      #plot(fit)
      
      rf.test$pred.label <- predict(fit ,rf.test)
      temp <- confusionMatrix(data=rf.test$pred.label, reference=rf.test$label)
      
      pred.accuracy.temp[j] <- temp$overall[1]
      pred.accuracy.upper.temp[j] <- temp$overall[4]
      pred.accuracy.lower.temp[j] <- temp$overall[3]
      sensitivity.temp[j] <- temp$byClass[1]
      specificity.temp[j] <- temp$byClass[2]
      
      var.imp <- data.frame(importance(fit,type=2))
      # make row names as columns
      var.imp$Variables <- row.names(var.imp)
      rownames(var.imp) <- NULL
      var.imp <- var.imp[order(var.imp$MeanDecreaseGini,decreasing = T),]
      
      if(j == 1){
        var.imp.temp <- var.imp
      }
      else{
        var.imp.temp <- merge(var.imp.temp, var.imp, by="Variables")
      }
      
    }
    
    if(i == 1){
      var.importance <- var.imp.temp
    }
    else{
      var.importance <- merge(var.importance, var.imp.temp, by="Variables")
    }
    
    pred.accuracy[i] <- mean(pred.accuracy.temp)
    pred.accuracy.upper[i] <- mean(pred.accuracy.upper.temp)
    pred.accuracy.lower[i] <- mean(pred.accuracy.lower.temp)
    mcnemar.p[i] <- mean(mcnemar.p.temp)
    sensitivity[i] <- mean(sensitivity.temp)
    specificity[i] <- mean(specificity.temp)
  }
  
  rownames(var.importance) <- var.importance$Variables
  var.importance <- var.importance[,-1]
  var.importance.value <- as.matrix(apply(var.importance,1,mean))
  
  result <- list(mean(pred.accuracy),
                 mean(pred.accuracy.upper),
                 mean(pred.accuracy.lower),
                 mean(mcnemar.p),
                 mean(sensitivity),
                 mean(specificity),
                 var.importance.value)
  names(result) <- c("accuracy","accuracy.upper","accuracy.lower","mcnemar.p",
                     "sensitivity","specificity","var.imp")
  return(result)
}

cv.result <- cross.validation(rf.data, nfolds = 5, iter = 10)
accuracy.1 <- cv.result$accuracy # 0.7204027
accuracy.upper.1 <- cv.result$accuracy.upper # 0.7526738
accuracy.lower.1 <- cv.result$accuracy.lower # 0.6863406
mcnemar.p.1 <- cv.result$mcnemar.p # 0
sensitivity.1 <- cv.result$sensitivity # 0.6837484
specificity.1 <- cv.result$specificity # 0.7535395
var.imp.1 <- cv.result$var.imp

#varImpPlot(fit, sort = T, main="Variable Importance", n.var=20)

#var.imp <- data.frame(importance(fit,type=2))
# make row names as columns
#var.imp$Variables <- row.names(var.imp)
#rownames(var.imp) <- NULL
#var.imp <- var.imp[order(var.imp$MeanDecreaseGini,decreasing = T),]

var.imp.1 <- as.data.frame(var.imp.1)
var.imp.1$names <- row.names(var.imp.1)
colnames(var.imp.1) <- c("MeanDecreaseGini", "Variables")
var.imp.1 <- var.imp.1[order(var.imp.1$MeanDecreaseGini,decreasing = T),]

plot.data <- var.imp.1[1:20,]
plot.data$Variables <- factor(plot.data$Variables, levels = rev(plot.data$Variables))
png("Figure_6a.png", height = 8, width = 8, units = "in", res = 300)
p1 <- ggplot(data=plot.data, aes(x=Variables, y=MeanDecreaseGini))+
  geom_point(shape=1)+
  coord_flip()+
  ylim(0,150)+
  theme_bw()+
  ylab("Mean decrease Gini value")+
  xlab("")+
  ggtitle("Sensitivity/hurt feelings")+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) 
p1
dev.off()

fit <- randomForest(label ~ ., data = rf.data, ntree=500, importance=T,
                    na.action=na.omit)

data.test.scaled.biomarker <- as.data.frame(lapply(data.test[,2:3298], scale))
data.test.feature <- data.test[,3299:ncol(data.test)]
data.test.fit <- cbind(data.test.scaled.biomarker[,selected.biomarker], 
                       data.test.feature[,feature.20])
data.test.fit[data.test.fit==-3] <- NA
data.test.fit[data.test.fit==-1] <- NA
data.test.fit[data.test.fit==-6] <- NA
colnames(data.test.fit) <- c(biomarker.names, 
                             "sex",
                             "Worrier_anxious_feelings",
                             "Risk_taking",
                             "Guilty_feelings",
                             "Alcohol_usually_taken_with_meals",
                             "Seen_doctor_for_nerves_anxiety_tension_depression",
                             "Worry_too_long_after_embarrassment",
                             "Miserableness",
                             "Snoring",
                             "Nervous_feelings",
                             "Ever_highly_irritable_argumentative_for_2_days",
                             "Ever_depressed_for_a_whole_week",
                             "Ever_unenthusiastic_disinterested_for_a_whole_week",
                             "Sleeplessness_insomnia",
                             "Getting_up_in_morning",
                             "Nap_during_day",
                             "Frequency_of_tiredness_lethargy_in_last_2_weeks",
                             "Alcohol_drinker_status",
                             "Past_tobacco_smoking",
                             "Frequency_of_depressed_mood_in_last_2_weeks")  
pred.label <- predict(fit ,data.test.fit)
true.label <- data.test[,"X1950.0.0"]
true.label[true.label=="-3"] <- NA
true.label[true.label=="-1"] <- NA
true.label <- as.factor(true.label)
confusionMatrix(pred.label, true.label) 


#################################################################################
##### Ever depressed for a whole week 2 #########################################
feature.20 <- c("X31.0.0","X1950.0.0","X1980.0.0","X2040.2.0","X2030.0.0",
                "X1618.2.0","X2000.0.0","X1930.0.0",
                "X1210.0.0","X1970.2.0","X4653.2.0",
                "X4631.2.0","X1200.0.0","X1170.0.0",
                "X1190.2.0","X2080.0.0","X20117.0.0","X2877.0.0",
                "X1249.2.0")

num.miss <- numeric(length(feature.20))
for(i in 1:length(feature.20)){
  num.miss[i] <- length(which(is.na(feature[,feature.20[i]])))
}
index <- which(num.miss > 1000)
feature.20 <- feature.20[-index]


rf.data <- cbind(scaled.biomarker[,selected.biomarker],feature[,feature.20])
rf.data[rf.data==-3] <- NA
rf.data[rf.data==-1] <- NA
rf.data[rf.data==-6] <- NA
colnames(rf.data) <- c(biomarker.names, 
                       "sex",
                       "Sensitivity_hurt_feelings",
                       "Worrier_anxious_feelings",
                       "Risk_taking",
                       "Guilty_feelings",
                       "Alcohol_usually_taken_with_meals",
                       "Worry_too_long_after_embarrassment",
                       "Miserableness",
                       "Snoring",
                       "Nervous_feelings",
                       "Ever_highly_irritable_argumentative_for_2_days",
                       "Ever_unenthusiastic_disinterested_for_a_whole_week",
                       "Sleeplessness_insomnia",
                       "Getting_up_in_morning",
                       "Nap_during_day",
                       "Frequency_of_tiredness_lethargy_in_last_2_weeks",
                       "Alcohol_drinker_status",
                       "Past_tobacco_smoking")  

rf.data$label <- feature[,"X4598.2.0"]
rf.data$label[rf.data$label=="-3"] <- NA
rf.data$label[rf.data$label=="-1"] <- NA
rf.data$label <- as.factor(rf.data$label)

set.seed(1234)
cv.result <- cross.validation(rf.data, nfolds = 5, iter = 10)
accuracy.3 <- cv.result$accuracy # 0.777967
accuracy.upper.3 <- cv.result$accuracy.upper # 0.8073186
accuracy.lower.3 <- cv.result$accuracy.lower # 0.7463962
mcnemar.p.3 <- cv.result$mcnemar.p # 0
sensitivity.3 <- cv.result$sensitivity # 0.9115463
specificity.3 <- cv.result$specificity # 0.6396355
var.imp.3 <- cv.result$var.imp

#varImpPlot(fit, sort = T, main="Variable Importance", n.var=20)

#var.imp <- data.frame(importance(fit,type=2))
# make row names as columns
#var.imp$Variables <- row.names(var.imp)
#rownames(var.imp) <- NULL
#var.imp <- var.imp[order(var.imp$MeanDecreaseGini,decreasing = T),]

var.imp.3 <- as.data.frame(var.imp.3)
var.imp.3$names <- row.names(var.imp.3)
colnames(var.imp.3) <- c("MeanDecreaseGini", "Variables")
var.imp.3 <- var.imp.3[order(var.imp.3$MeanDecreaseGini,decreasing = T),]

plot.data <- var.imp.3[1:20,]
plot.data$Variables <- factor(plot.data$Variables, levels = rev(plot.data$Variables))
png("Figure_6b.png", height = 8, width = 8, units = "in", res = 300)
p2 <- ggplot(data=plot.data, aes(x=Variables, y=MeanDecreaseGini))+
  geom_point(shape=1)+
  coord_flip()+
  ylim(0,400)+
  theme_bw()+
  ylab("Mean decrease Gini value")+
  xlab("")+
  ggtitle("Ever depressed for a whole week")+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) 
p2
dev.off()

fit <- randomForest(label ~ ., data = rf.data, ntree=500, importance=T,
                    na.action=na.omit)

data.test.scaled.biomarker <- as.data.frame(lapply(data.test[,2:3298], scale))
data.test.feature <- data.test[,3299:ncol(data.test)]
data.test.fit <- cbind(data.test.scaled.biomarker[,selected.biomarker], 
                       data.test.feature[,feature.20])
data.test.fit[data.test.fit==-3] <- NA
data.test.fit[data.test.fit==-1] <- NA
data.test.fit[data.test.fit==-6] <- NA
colnames(data.test.fit) <- c(biomarker.names, 
                             "sex",
                             "Sensitivity_hurt_feelings",
                             "Worrier_anxious_feelings",
                             "Risk_taking",
                             "Guilty_feelings",
                             "Alcohol_usually_taken_with_meals",
                             "Worry_too_long_after_embarrassment",
                             "Miserableness",
                             "Snoring",
                             "Nervous_feelings",
                             "Ever_highly_irritable_argumentative_for_2_days",
                             "Ever_unenthusiastic_disinterested_for_a_whole_week",
                             "Sleeplessness_insomnia",
                             "Getting_up_in_morning",
                             "Nap_during_day",
                             "Frequency_of_tiredness_lethargy_in_last_2_weeks",
                             "Alcohol_drinker_status",
                             "Past_tobacco_smoking")  
pred.label <- predict(fit ,data.test.fit)
true.label <- data.test[,"X4598.2.0"]
true.label[true.label=="-3"] <- NA
true.label[true.label=="-1"] <- NA
true.label <- as.factor(true.label)
confusionMatrix(pred.label, true.label) 


###############################################################################
##### Worrier / anxious feelings ##############################################
feature.20 <- c("X31.0.0","X1950.0.0","X2040.2.0","X2030.0.0",
                "X1618.2.0","X2090.0.0","X2000.0.0","X1930.0.0",
                "X1210.0.0","X1970.2.0","X4653.2.0","X4598.2.0",
                "X4631.2.0","X1200.0.0","X1170.0.0",
                "X1190.2.0","X2080.0.0","X20117.0.0","X2877.0.0",
                "X1249.2.0","X2050.0.0")

num.miss <- numeric(length(feature.20))
for(i in 1:length(feature.20)){
  num.miss[i] <- length(which(is.na(feature[,feature.20[i]])))
}
index <- which(num.miss > 1000)
feature.20 <- feature.20[-index]


rf.data <- cbind(scaled.biomarker[,selected.biomarker],feature[,feature.20])
rf.data[rf.data==-3] <- NA
rf.data[rf.data==-1] <- NA
rf.data[rf.data==-6] <- NA
colnames(rf.data) <- c(biomarker.names, 
                       "sex",
                       "Sensitivity_hurt_feelings",
                       "Risk_taking",
                       "Guilty_feelings",
                       "Alcohol_usually_taken_with_meals",
                       "Seen_doctor_for_nerves_anxiety_tension_depression",
                       "Worry_too_long_after_embarrassment",
                       "Miserableness",
                       "Snoring",
                       "Nervous_feelings",
                       "Ever_highly_irritable_argumentative_for_2_days",
                       "Ever_depressed_for_a_whole_week",
                       "Ever_unenthusiastic_disinterested_for_a_whole_week",
                       "Sleeplessness_insomnia",
                       "Getting_up_in_morning",
                       "Nap_during_day",
                       "Frequency_of_tiredness_lethargy_in_last_2_weeks",
                       "Alcohol_drinker_status",
                       "Past_tobacco_smoking",
                       "Frequency_of_depressed_mood_in_last_2_weeks")  

rf.data$label <- feature[,"X1980.0.0"]
rf.data$label[rf.data$label=="-3"] <- NA
rf.data$label[rf.data$label=="-1"] <- NA
rf.data$label <- as.factor(rf.data$label)

set.seed(1234)
cv.result <- cross.validation(rf.data, nfolds = 5, iter = 10)
accuracy.4 <- cv.result$accuracy # 0.7393196
accuracy.upper.4 <- cv.result$accuracy.upper # 0.7708003
accuracy.lower.4 <- cv.result$accuracy.lower # 0.7058916
mcnemar.p.4 <- cv.result$mcnemar.p # 0
sensitivity.4 <- cv.result$sensitivity # 0.7230401
specificity.4 <- cv.result$specificity # 0.7548048
var.imp.4 <- cv.result$var.imp

#varImpPlot(fit, sort = T, main="Variable Importance", n.var=20)

#var.imp <- data.frame(importance(fit,type=2))
# make row names as columns
#var.imp$Variables <- row.names(var.imp)
#rownames(var.imp) <- NULL
#var.imp <- var.imp[order(var.imp$MeanDecreaseGini,decreasing = T),]

var.imp.4 <- as.data.frame(var.imp.4)
var.imp.4$names <- row.names(var.imp.4)
colnames(var.imp.4) <- c("MeanDecreaseGini", "Variables")
var.imp.4 <- var.imp.4[order(var.imp.4$MeanDecreaseGini,decreasing = T),]

plot.data <- var.imp.4[1:20,]
plot.data$Variables <- factor(plot.data$Variables, levels = rev(plot.data$Variables))
png("Figure_6c.png", height = 8, width = 8, units = "in", res = 300)
p3 <- ggplot(data=plot.data, aes(x=Variables, y=MeanDecreaseGini))+
  geom_point(shape=1)+
  coord_flip()+
  ylim(0,105)+
  theme_bw()+
  ylab("Mean decrease Gini value")+
  xlab("")+
  ggtitle("Worrier/anxious feelings")+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) 
p3
dev.off()

fit <- randomForest(label ~ ., data = rf.data, ntree=500, importance=T,
                    na.action=na.omit)

data.test.scaled.biomarker <- as.data.frame(lapply(data.test[,2:3298], scale))
data.test.feature <- data.test[,3299:ncol(data.test)]
data.test.fit <- cbind(data.test.scaled.biomarker[,selected.biomarker], 
                       data.test.feature[,feature.20])
data.test.fit[data.test.fit==-3] <- NA
data.test.fit[data.test.fit==-1] <- NA
data.test.fit[data.test.fit==-6] <- NA
colnames(data.test.fit) <- c(biomarker.names, 
                             "sex",
                             "Sensitivity_hurt_feelings",
                             "Risk_taking",
                             "Guilty_feelings",
                             "Alcohol_usually_taken_with_meals",
                             "Seen_doctor_for_nerves_anxiety_tension_depression",
                             "Worry_too_long_after_embarrassment",
                             "Miserableness",
                             "Snoring",
                             "Nervous_feelings",
                             "Ever_highly_irritable_argumentative_for_2_days",
                             "Ever_depressed_for_a_whole_week",
                             "Ever_unenthusiastic_disinterested_for_a_whole_week",
                             "Sleeplessness_insomnia",
                             "Getting_up_in_morning",
                             "Nap_during_day",
                             "Frequency_of_tiredness_lethargy_in_last_2_weeks",
                             "Alcohol_drinker_status",
                             "Past_tobacco_smoking",
                             "Frequency_of_depressed_mood_in_last_2_weeks")  
pred.label <- predict(fit ,data.test.fit)
true.label <- data.test[,"X1980.0.0"]
true.label[true.label=="-3"] <- NA
true.label[true.label=="-1"] <- NA
true.label <- as.factor(true.label)
confusionMatrix(pred.label, true.label) 

###############################################################################
##### Miserableness ###########################################################
feature.20 <- c("X31.0.0","X1950.0.0","X1980.0.0","X2040.2.0","X2030.0.0",
                "X1618.2.0","X2090.0.0","X2000.0.0",
                "X1210.0.0","X1970.2.0","X4653.2.0","X4598.2.0",
                "X4631.2.0","X1200.0.0","X1170.0.0",
                "X1190.2.0","X2080.0.0","X20117.0.0","X2877.0.0",
                "X1249.2.0","X2050.0.0")

num.miss <- numeric(length(feature.20))
for(i in 1:length(feature.20)){
  num.miss[i] <- length(which(is.na(feature[,feature.20[i]])))
}
index <- which(num.miss > 1000)
feature.20 <- feature.20[-index]


rf.data <- cbind(scaled.biomarker[,selected.biomarker],feature[,feature.20])
rf.data[rf.data==-3] <- NA
rf.data[rf.data==-1] <- NA
rf.data[rf.data==-6] <- NA
colnames(rf.data) <- c(biomarker.names, 
                       "sex",
                       "Sensitivity_hurt_feelings",
                       "Worrier_anxious_feelings",
                       "Risk_taking",
                       "Guilty_feelings",
                       "Alcohol_usually_taken_with_meals",
                       "Seen_doctor_for_nerves_anxiety_tension_depression",
                       "Worry_too_long_after_embarrassment",
                       "Snoring",
                       "Nervous_feelings",
                       "Ever_highly_irritable_argumentative_for_2_days",
                       "Ever_depressed_for_a_whole_week",
                       "Ever_unenthusiastic_disinterested_for_a_whole_week",
                       "Sleeplessness_insomnia",
                       "Getting_up_in_morning",
                       "Nap_during_day",
                       "Frequency_of_tiredness_lethargy_in_last_2_weeks",
                       "Alcohol_drinker_status",
                       "Past_tobacco_smoking",
                       "Frequency_of_depressed_mood_in_last_2_weeks")  

rf.data$label <- feature[,"X1930.0.0"]
rf.data$label[rf.data$label=="-3"] <- NA
rf.data$label[rf.data$label=="-1"] <- NA
rf.data$label <- as.factor(rf.data$label)

set.seed(1234)
cv.result <- cross.validation(rf.data, nfolds = 5, iter = 10)
accuracy.5 <- cv.result$accuracy # 0.7432627
accuracy.upper.5 <- cv.result$accuracy.upper # 0.7745725
accuracy.lower.5 <- cv.result$accuracy.lower # 0.7099731
mcnemar.p.5 <- cv.result$mcnemar.p # 0
sensitivity.5 <- cv.result$sensitivity # 0.8666349
specificity.5 <- cv.result$specificity # 0.5498084
var.imp.5 <- cv.result$var.imp

#varImpPlot(fit, sort = T, main="Variable Importance", n.var=20)

#var.imp <- data.frame(importance(fit,type=2))
# make row names as columns
#var.imp$Variables <- row.names(var.imp)
#rownames(var.imp) <- NULL
#var.imp <- var.imp[order(var.imp$MeanDecreaseGini,decreasing = T),]

var.imp.5 <- as.data.frame(var.imp.5)
var.imp.5$names <- row.names(var.imp.5)
colnames(var.imp.5) <- c("MeanDecreaseGini", "Variables")
var.imp.5 <- var.imp.5[order(var.imp.5$MeanDecreaseGini,decreasing = T),]

plot.data <- var.imp.5[1:20,]
plot.data$Variables <- factor(plot.data$Variables, levels = rev(plot.data$Variables))
png("Figure_6d.png", height = 8, width = 8, units = "in", res = 300)
p4 <- ggplot(data=plot.data, aes(x=Variables, y=MeanDecreaseGini))+
  geom_point(shape=1)+
  coord_flip()+
  ylim(0,105)+
  theme_bw()+
  ylab("Mean decrease Gini value")+
  xlab("")+
  ggtitle("Miserableness")+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) 
p4
dev.off()

fit <- randomForest(label ~ ., data = rf.data, ntree=500, importance=T,
                    na.action=na.omit)

data.test.scaled.biomarker <- as.data.frame(lapply(data.test[,2:3298], scale))
data.test.feature <- data.test[,3299:ncol(data.test)]
data.test.fit <- cbind(data.test.scaled.biomarker[,selected.biomarker], 
                       data.test.feature[,feature.20])
data.test.fit[data.test.fit==-3] <- NA
data.test.fit[data.test.fit==-1] <- NA
data.test.fit[data.test.fit==-6] <- NA
colnames(data.test.fit) <- c(biomarker.names, 
                             "sex",
                             "Sensitivity_hurt_feelings",
                             "Worrier_anxious_feelings",
                             "Risk_taking",
                             "Guilty_feelings",
                             "Alcohol_usually_taken_with_meals",
                             "Seen_doctor_for_nerves_anxiety_tension_depression",
                             "Worry_too_long_after_embarrassment",
                             "Snoring",
                             "Nervous_feelings",
                             "Ever_highly_irritable_argumentative_for_2_days",
                             "Ever_depressed_for_a_whole_week",
                             "Ever_unenthusiastic_disinterested_for_a_whole_week",
                             "Sleeplessness_insomnia",
                             "Getting_up_in_morning",
                             "Nap_during_day",
                             "Frequency_of_tiredness_lethargy_in_last_2_weeks",
                             "Alcohol_drinker_status",
                             "Past_tobacco_smoking",
                             "Frequency_of_depressed_mood_in_last_2_weeks")
pred.label <- predict(fit ,data.test.fit)
true.label <- data.test[,"X1930.0.0"]
true.label[true.label=="-3"] <- NA
true.label[true.label=="-1"] <- NA
true.label <- as.factor(true.label)
confusionMatrix(pred.label, true.label) 

png("Figure_6.png", height = 10, width = 15, units = "in", res = 300)
plot_grid(p1,p2,p3,p4, labels = "auto", ncol = 2)
dev.off()




###############################################################################
###############################################################################
## 20 biomarkers, 20 features, prediction accuracy 0.90 (0.86,0.92)

options(repos='http://cran.rstudio.org')
have.packages <- installed.packages()
cran.packages <- c('devtools','plotrix','randomForest','tree')
to.install <- setdiff(cran.packages, have.packages[,1])
if(length(to.install)>0) install.packages(to.install)

library(devtools)
if(!('reprtree' %in% installed.packages())){
  install_github('araastat/reprtree')
}
for(p in c(cran.packages, 'reprtree')) eval(substitute(library(pkg), list(pkg=p)))

library(randomForest)
library(reprtree)
library(partykit)

#selected.biomarker <- biomarker.20[c(1,10:15,17:18),1]
biomarker.names <- c("aseg__MaskVol", #*
                     "aseg__TotalGrayVol", #*
                     "rh_aparc.DKTatlas_area__rh_superiortemporal_area", #*
                     "rh_aparc.DKTatlas_area__rh_superiorfrontal_area", #*
                     "lh_aparc.DKTatlas_area__lh_superiortemporal_area") #*
# biomarker.names <- c("rh_BA_exvivo_area__rh_WhiteSurfArea_area",
#                      "lh_BA_exvivo_area__lh_WhiteSurfArea_area",
#                      "rh_aparc_area__rh_WhiteSurfArea_area",
#                      "rh_aparc.a2009s_area__rh_WhiteSurfArea_area",
#                      "lh_aparc_area__lh_WhiteSurfArea_area",
#                      "lh_aparc.a2009s_area__lh_WhiteSurfArea_area",
#                      "aseg__SupraTentorialVol",
#                      "aseg__SupraTentorialVolNotVent",
#                      "aseg__SupraTentorialVolNotVentVox",
#                      "aseg__BrainSegVol",
#                      "aseg__BrainSegVolNotVent",
#                      "aseg__BrainSegVolNotVentSurf",
#                      "aseg__CortexVol",
#                      "aseg__rhCortexVol",
#                      "aseg__MaskVol", #*
#                      "aseg__lhCortexVol",
#                      "aseg__TotalGrayVol", #*
#                      "rh_aparc.DKTatlas_area__rh_superiortemporal_area", #*
#                      "rh_aparc.DKTatlas_area__rh_superiorfrontal_area", #*
#                      "lh_aparc.DKTatlas_area__lh_superiortemporal_area") #*

feature.one <- c("X2000.0.0","X1980.0.0")  ##unenthusiastic

rf.data <- cbind(scaled.biomarker[,selected.biomarker[1:5]],feature[,feature.one])
rf.data[rf.data==-3] <- NA
rf.data[rf.data==-1] <- NA
rf.data[rf.data==-6] <- NA
colnames(rf.data) <- c(biomarker.names, "Worry_too_long_after_embarrassment","Worrier_anxious_feelings")

rf.data$label <- feature[,"X1950.0.0"]
rf.data$label[rf.data$label==-3] <- NA
rf.data$label[rf.data$label==-1] <- NA
rf.data$label <- as.factor(rf.data$label)

set.seed(1234)
sample.ind <- sample(2,nrow(rf.data), replace = T, prob = c(0.7,0.3))
train.data <- rf.data[sample.ind==1,]
valid.data <- rf.data[sample.ind==2,]

# train.data.complete <- train.data[which(!is.na(train.data$label)),]
# x <- ctree(label ~ ., data = train.data.complete)
# png("forest_plot.png", height = 10, width = 20, units = "in", res = 300)
# plot(x, type="simple")
# dev.off()

library(party)
train.data.complete <- train.data[which(!is.na(train.data$label)),]
x <- ctree(label ~ ., data = train.data.complete)
png("forest_plot.png", height = 10, width = 20, units = "in", res = 300)
plot(x, type="simple")
dev.off()

#################################################################################
#################################################################################
############################### t-SNE plot (2D, 3D) #############################
labels <- cluster_2$cluster
labels <- as.factor(labels)
colors <- c("red", "blue")
names(colors) <- unique(labels)

## perplexity in the range 5-50
## steps 10-1000 (I will try 500 1000 5000)

p.seq <- seq(5,50,by=5)
iter.seq <- c(500,1000,5000)

for(i in 1:length(iter.seq)){
  for(j in 1:length(p.seq)){
    
    i <- 1
    j <- 5 ##1-5 1-7
    
    set.seed(1234)
    tsne <- Rtsne(scaled.biomarker[,-3298], dims = 2, perplexity=p.seq[j], max_iter = iter.seq[i], 
                  verbose=TRUE, theta=0.5, pca=FALSE)
    
    m.data <- as.data.frame(tsne$Y)
    colnames(m.data) <- c("Coordinate1", "Coordinate2")
    m.data$cluster <- labels 
    m.data$cluster <- as.factor(m.data$cluster)
    
    png("Figure_2c.png", height = 6, width = 8, units = "in", res = 300)
    ggplot(m.data, aes(x=Coordinate1, y=Coordinate2, shape=cluster, color=cluster))+
      geom_point()+
      scale_color_manual(values=c("red2", "royalblue2"))+
      theme_bw()+
      theme(legend.position = "bottom",
            text = element_text(size=20),
            legend.text=element_text(size=18))
    dev.off()
    
    #tsne$origD
    #str(tsne)
    
    ## 128 points
    index <- which((m.data$Coordinate1 > -8 & m.data$Coordinate1 < -3) &
                     (m.data$Coordinate2 > 10.5 & m.data$Coordinate2 < 15))
    
    m.data$cluster3 <- as.numeric(as.character(m.data$cluster))
    m.data$cluster3[index] <- 3
    m.data$cluster3 <- as.factor(m.data$cluster3)
    
    png("Figure_S5.png", height = 6, width = 8, units = "in", res = 300)
    ggplot(m.data, aes(x=Coordinate1, y=Coordinate2, shape=cluster3, color=cluster3))+
      geom_point()+
      scale_color_manual(values=c("red2", "royalblue2", "green"))+
      theme_bw()+
      theme(legend.position = "bottom",
            text = element_text(size=20),
            legend.text=element_text(size=18))
    dev.off()
    
    label <- m.data$cluster3
    
  }
}
#####################################################################################
#####################################################################################

set.seed(1234)
tsne <- Rtsne(scaled.biomarker[,-3298], dims = 3, perplexity=25, max_iter = 500, 
              verbose=TRUE, theta=0.5, pca=FALSE)

m.data <- as.data.frame(tsne$Y)
colnames(m.data) <- c("Coordinate1", "Coordinate2", "Coordinate3")
m.data$cluster <- labels
m.data$cluster <- as.factor(m.data$cluster)

ggplot(m.data.scaled, aes(x=Y1, y=Y2, shape=cluster, color=cluster))+
  geom_point()+
  theme_bw()

p <- plot_ly(m.data, x = ~Coordinate1, y = ~Coordinate2, z = ~Coordinate3, 
             color = ~cluster, colors = c("red2", "royalblue2")) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = "Coordinate1"),
                      yaxis = list(title = "Coordinate2"),
                      zaxis = list(title = "Coordinate3")))
p

#################################################################################
#################################################################################
##################### three clusters according to the t-SNE plot ################
##################### compare the blue cluster and green cluster ################
#################################################################################
## the top 20 selected biomarkers
setwd("~/Google Drive/Yiwang_Zhou/Pilot_Projects/UKBB_TM_NLP/clustering_results/Manuscript_revision/result")
data.summary.cluster.2 <- read.csv("data_summary_cluster_2.csv", header = TRUE)
top.selected <- as.character(data.summary.cluster.2[1:20,1])

#scaled.biomarker$cluster_2 <- as.factor(scaled.biomarker$cluster_2)
scaled.biomarker$cluster_3 <- as.numeric(as.character(label))  ##label is got from tSNE plot

bar.data <- scaled.biomarker[, top.selected]
colnames(bar.data) <- c(paste("rh_BA_exvivo_area_", "_rh_WhiteSurfArea_area", sep = "\n"),
                        paste("lh_BA_exvivo_area_", "_lh_WhiteSurfArea_area", sep = "\n"),
                        paste("rh_aparc_area_", "_lh_WhiteSurfArea_area", sep = "\n"),
                        paste("rh_aparc.a2009s_area_", "_rh_WhiteSurfArea_area", sep = "\n"),
                        paste("lh_aparc_area__lh_", "WhiteSurfArea_area", sep = "\n"),
                        paste("lh_aparc.a2009s_area_", "_lh_WhiteSurfArea_area", sep = "\n"),
                        "aseg__SupraTentorialVol",
                        "aseg__SupraTentorialVolNotVent",
                        "aseg__SupraTentorialVolNotVentVox",
                        "aseg__BrainSegVol",
                        "aseg__BrainSegVolNotVent",
                        "aseg__BrainSegVolNotVentSurf",
                        "aseg__CortexVol",
                        "aseg__rhCortexVol",
                        "aseg__MaskVol",
                        "aseg__lhCortexVol",
                        "aseg__TotalGrayVol",
                        paste("rh_aparc.DKTatlas_area_", "_rh_superiortemporal_area", sep = "\n"),
                        paste("rh_aparc.DKTatlas_area_", "_rh_superiorfrontal_area", sep = "\n"),
                        paste("lh_aparc.DKTatlas_area_", "_lh_superiortemporal_area", sep = "\n"))

bar.data$id <- 1:nrow(bar.data)
bar.data.m <- melt(bar.data, "id")
bar.data.m$cluster <- rep(scaled.biomarker$cluster_3, times=20)
bar.data.m$cluster <- as.factor(bar.data.m$cluster)

png("Figure_S6a.png", height = 8, width = 15, units = "in", res = 300)
ggplot(data = bar.data.m, aes(x=value))+
  geom_density(aes(fill=cluster),alpha = 0.5,linetype="solid", kernel="epanechnikov", bw=0.3)+
  facet_wrap( ~ variable, ncol=5)+
  xlim(-5,5)+
  ylim(0,0.8)+
  scale_fill_manual(name=NULL,values=c("red2", "royalblue2","green"),labels=c("cluster1", "cluster2","cluster3"))+
  theme(legend.position = "bottom",
        text = element_text(size=13),
        legend.text=element_text(size=14),
        panel.background = element_rect(fill = NA),
        axis.text = element_text(size=12),
        axis.title=element_text(size=14))
dev.off()
######################################################################################
## blue cluster
cluster2 <- scaled.biomarker[which(scaled.biomarker$cluster_3 == 2),]
## green cluster
cluster3 <- scaled.biomarker[which(scaled.biomarker$cluster_3 == 3),]

cluster2 <- cluster2[,1:(ncol(cluster2)-2)]
n2 <- nrow(cluster2)
cluster3 <- cluster3[,1:(ncol(cluster3)-2)]
n3 <- nrow(cluster3)

mean2 <- apply(cluster2, 2, mean)
sd2 <- apply(cluster2, 2, sd)
median2 <- apply(cluster2, 2, median)
range2 <- apply(cluster2, 2, function(x) max(x)-min(x))
skew2 <- apply(cluster2, 2, function(x) skewness(x))
kurtosis2 <- apply(cluster2, 2, function(x) kurtosis(x))

mean3 <- apply(cluster3, 2, mean)
sd3 <- apply(cluster3, 2, sd)
median3 <- apply(cluster3, 2, median)
range3 <- apply(cluster3, 2, function(x) max(x)-min(x))
skew3 <- apply(cluster3, 2, function(x) skewness(x))
kurtosis3 <- apply(cluster3, 2, function(x) kurtosis(x))

## t-test
ttest <- function(i){
  model <- t.test(cluster2[,i], cluster3[,i], alternative = "two.sided", var.equal = FALSE)
  t <- model$statistic
  p <- model$p.value
  
  result <- list(t, p)
  names(result) <- c("t","p")
  return(result)
}

ttest_result <- lapply(1:ncol(biomarker), function(i) ttest(i))
t_test_value <- do.call(c, lapply(1:ncol(biomarker), function(i) ttest_result[[i]]$t))
t_test_p <- do.call(c, lapply(1:ncol(biomarker), function(i) ttest_result[[i]]$p))


## KS test
kstest <- function(i){
  model <- ks.test(cluster2[,i], cluster3[,i], alternative = "two.sided", exact = NULL)
  ks <- model$statistic
  p <- model$p.value
  
  result <- list(ks, p)
  names(result) <- c("ks","p")
  return(result)
}

kstest_result <- lapply(1:ncol(biomarker), function(i) kstest(i))
ks_test_value <- do.call(c, lapply(1:ncol(biomarker), function(i) kstest_result[[i]]$ks))
ks_test_p <- do.call(c, lapply(1:ncol(biomarker), function(i) kstest_result[[i]]$p))


## Wilcoxon Test 
wiltest <- function(i){
  model <- wilcox.test(cluster2[,i], cluster3[,i], alternative = "two.sided")
  wil <- model$statistic
  p <- model$p.value
  
  result <- list(wil, p)
  names(result) <- c("wil","p")
  return(result)
}

wiltest_result <- lapply(1:ncol(biomarker), function(i) wiltest(i))
wil_test_value <- do.call(c, lapply(1:ncol(biomarker), function(i) wiltest_result[[i]]$wil))
wil_test_p <- do.call(c, lapply(1:ncol(biomarker), function(i) wiltest_result[[i]]$p))

data.summary.cluster.3 <- data.frame("mean_cluster2"=mean2,
                                     "sd_cluster2"=sd2,
                                     "median_cluster2"=median2,
                                     "range_cluster2"=range2,
                                     "skew_cluster2"=skew2,
                                     "kurtosis_cluster2"=kurtosis2,
                                     "mean_cluster3"=mean3,
                                     "sd_cluster3"=sd3,
                                     "median_cluster3"=median3,
                                     "range_cluster3"=range3,
                                     "skew_cluster3"=skew3,
                                     "kurtosis_cluster3"=kurtosis3,
                                     "t_value"=t_test_value,
                                     "t_p"=t_test_p,
                                     "ks_value"=ks_test_value,
                                     "ks_p"=ks_test_p,
                                     "wil_value"=wil_test_value,
                                     "wil_p"=wil_test_p)
#####################################################################################
#####################################################################################
## Bonferroni correction of the threshold of the p-values
p.thresh <- 0.05/dim(biomarker)[2]
length(which(data.summary.cluster.3$t_p < p.thresh))
index.t <- which(data.summary.cluster.3$t_p < p.thresh)
sig.t <- numeric(nrow(data.summary.cluster.3))
sig.t[index.t] <- 1
data.summary.cluster.3$sig.t <- sig.t


length(which(data.summary.cluster.3$ks_p < p.thresh))
index.ks <- which(data.summary.cluster.3$ks_p < p.thresh)
sig.ks <- numeric(nrow(data.summary.cluster.3))
sig.ks[index.ks] <- 1
data.summary.cluster.3$sig.ks <- sig.ks

length(which(data.summary.cluster.3$wil_p < p.thresh))
index.wil <- which(data.summary.cluster.3$wil_p < p.thresh)
sig.wil <- numeric(nrow(data.summary.cluster.3))
sig.wil[index.wil] <- 1
data.summary.cluster.3$sig.wil <- sig.wil


length(which(data.summary.cluster.3$sig.t == 1 
             & data.summary.cluster.3$sig.ks == 1
             & data.summary.cluster.3$sig.wil == 1))

data.summary.cluster.3 <- data.summary.cluster.3[order(abs(data.summary.cluster.3$t_value), decreasing = TRUE),]
data.summary.cluster.3$t.order <- seq(1,nrow(data.summary.cluster.3),1)
data.summary.cluster.3 <- data.summary.cluster.3[order(data.summary.cluster.3$ks_value, decreasing = TRUE),]
data.summary.cluster.3$ks.order <- seq(1,nrow(data.summary.cluster.3),1)
data.summary.cluster.3 <- data.summary.cluster.3[order(data.summary.cluster.3$wil_p, decreasing = FALSE),]
data.summary.cluster.3$wil.order <- seq(1,nrow(data.summary.cluster.3),1)

data.summary.cluster.3$sum.order <- data.summary.cluster.3$t.order+data.summary.cluster.3$ks.order+data.summary.cluster.3$wil.order
data.summary.cluster.3 <- data.summary.cluster.3[order(data.summary.cluster.3$sum.order),]
write.csv(data.summary.cluster.3, file = "data_summary_cluster_3_kmeans_blue_vs_green.csv")

######################################################################################
data.summary.cluster.3 <- read.csv("data_summary_cluster_3_kmeans_blue_vs_green.csv", header = TRUE)
top.10.selected <- data.summary.cluster.3[c(1,5,6,10,12,13,14,18,26,28),1]
top.10.selected <- as.character(top.10.selected)

bar.data <- scaled.biomarker[which(scaled.biomarker$cluster_3 %in% c(2,3)), top.10.selected]
colnames(bar.data) <- c(paste("rh_aparc.DKTatlas_area_","_rh_postcentral_area", sep = "\n"),
                        paste("rh_aparc.a2009s_area_","_rh_S_postcentral_area", sep = "\n"),
                        paste("rh_aparc_area_","_rh_postcentral_area", sep = "\n"),
                        paste("rh_aparc.a2009s_volume_","_rh_S_postcentral_volume", sep = "\n"),
                        paste("rh_aparc.pial_area_","_rh_postcentral_area", sep = "\n"),
                        paste("rh_aparc.a2009s_thickness_","_rh_S_precentral.sup.part_thickness", sep = "\n"),
                        paste("rh_aparc.a2009s_area_","_rh_G_postcentral_area", sep = "\n"),
                        paste("rh_aparc.a2009s_thickness_","_rh_G_precentral_thickness", sep = "\n"),
                        paste("rh_aparc_thickness_","_rh_precentral_thickness", sep = "\n"),
                        paste("rh_aparc.a2009s_volume_","_rh_G_precentral_volume", sep = "\n"))


bar.data$id <- 1:nrow(bar.data)
bar.data.m <- melt(bar.data, "id")
bar.data.m$cluster <- rep(scaled.biomarker$cluster_3[which(scaled.biomarker$cluster_3 %in% c(2,3))], times=10)
bar.data.m$cluster <- as.factor(bar.data.m$cluster)

png("Figure_S6b.png", height = 4, width = 15, units = "in", res = 300)
ggplot(data = bar.data.m, aes(x=value))+
  geom_density(aes(fill=cluster),alpha = 0.5,linetype="solid", kernel="epanechnikov", bw=0.4)+
  facet_wrap( ~ variable, ncol=5)+
  xlim(-5,5)+
  ylim(0,0.8)+
  scale_fill_manual(name=NULL,values=c("royalblue2", "green2"),labels=c("cluster2", "cluster3"))+
  theme(legend.position = "bottom",
        text = element_text(size=13),
        legend.text=element_text(size=14),
        panel.background = element_rect(fill = NA),
        axis.text = element_text(size=12),
        axis.title=element_text(size=14))
dev.off()


## summary statistics
data_summary <- biomarker[,top.10.selected]
data_summary_2 <- data_summary[which(scaled.biomarker$cluster_3 == 2),]
data_summary_3 <- data_summary[which(scaled.biomarker$cluster_3 == 3),]

summary_table <- data.frame("names"=colnames(data_summary),
                            "mean1"=apply(data_summary_2, 2, mean),
                            "median1"=apply(data_summary_2, 2, median),
                            "sd1"=apply(data_summary_2, 2, sd),
                            "mean2"=apply(data_summary_3, 2, mean),
                            "median2"=apply(data_summary_3, 2, median),
                            "sd2"=apply(data_summary_3, 2, sd))

rownames(summary_table) <- NULL
summary_table[c(1:5,7,10),2:7] <- sapply(2:7, function(i) round(summary_table[c(1:5,7,10),i], digits = 0))
summary_table[c(6,8,9),2:7] <- sapply(2:7, function(i) round(summary_table[c(6,8,9),i], digits = 2))
write.csv(summary_table, file = "Table2.csv")

