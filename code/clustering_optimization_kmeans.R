#####################################################################################
#####################################################################################
#####################################################################################
.libPaths("/home/yiwangz/UKBB/lib")
packages <- c("matrixStats","stats","cluster","ggplot2", "factoextra", 
              "NbClust", "cowplot")
ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.r-project.org/")
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)

#####################################################################################
data <- read.csv("Whole_dataset_ukbb_pilot.csv", header = TRUE, 
                 na.strings = c("", " ", "NA", "NaN", "Inf", "."))
n.sub <- nrow(data)

set.seed(1234)
index_train <- sample(seq(1, n.sub, 1), floor(0.80*n.sub), replace = FALSE)
data.train <- data[index_train,]
n.train <- nrow(data.train)

#### biomarker and feature are both from the training dataset
biomarker <- data.train[,2:3298]
scaled.biomarker <- as.data.frame(lapply(biomarker, scale))

############################# optimization by kmeans ##################################
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
# labs(subtitle = "Elbow method")+

#pdf("kmeans_tuning_gap_statistics.pdf")
#set.seed(1234)
#fviz_nbclust(scaled.biomarker, kmeans, nstart = 25,  method = "gap_stat", nboot = 500)+
  #labs(subtitle = "Gap statistic method")
#dev.off()

#####################################################################################
#####################################################################################
#####################################################################################
nb_kmeans <- NbClust(scaled.biomarker, distance = "euclidean", min.nc = 2,
              max.nc = 10, method = "kmeans")
#nb_hc.ward.D <- NbClust(scaled.biomarker, distance = "euclidean", min.nc = 2,
                     #max.nc = 5, method = "ward.D")
# nb_hc.ward.D2 <- NbClust(scaled.biomarker, distance = "euclidean", min.nc = 2,
#                  max.nc = 10, method = "ward.D2")

capture.output(print(fviz_nbclust(nb_kmeans)), file = "output_nb_kmeans.txt")
#capture.output(print(fviz_nbclust(nb_hc.ward.D2)), file = "output_nb_hc.ward.D2.txt")

#fviz_nbclust(nb_kmeans)
#fviz_nbclust(nb_hc.ward.D)
#fviz_nbclust(nb_hc.ward.D2)

#####################################################################################
#####################################################################################
#####################################################################################
## make the optimization plot by ggplot2
sil.value <- data.frame(matrix(0,ncol = 3, nrow = 20))
colnames(sil.value) <- c("k", "silhouette_value", "clustering_method")
sil.value$k <- rep(seq(1,10,by=1), times=2)
sil.value$clustering_method <- rep(c("k_means", "hierarchical"), each=10)
sil.value$clustering_method <- factor(sil.value$clustering_method, levels = c("k_means", "hierarchical"))
sil.value$silhouette_value[1:10] <- c(0, 0.074, 0.067, 0.031, 0.030, 0.026,
                                      0.020, 0.017, 0.013, 0.011)
sil.value$silhouette_value[11:20] <- c(0, 0.056, 0.055, 0.016, 0.010, 0.006,
                                       0.008, 0.009, 0.006, -0.004)

png("Figure_1.png", height = 4, width = 6, units = "in", res = 300)
ggplot(data=sil.value, aes(x=k, y=silhouette_value)) +
  geom_line(color="black")+
  geom_point(color="black")+
  facet_wrap( ~ clustering_method, ncol = 2)+
  theme_bw()+
  xlab("number of cluster")+
  ylab("average Silhouette value")+
  scale_x_continuous(breaks = seq(1,10,1))
dev.off()


#####################################################################################
#####################################################################################
wss.value <- data.frame(matrix(0,ncol = 3, nrow = 20))
colnames(wss.value) <- c("k", "wss", "clustering_method")
wss.value$k <- rep(seq(1,10,by=1), times=2)
wss.value$clustering_method <- rep(c("k_means", "hierarchical"), each=10)
wss.value$clustering_method <- factor(wss.value$clustering_method, levels = c("k_means", "hierarchical"))
wss.value$wss[1:10] <- c(2.615e7, 2.425e7, 2.312e7, 2.264e7, 2.226e7, 
                         2.202e7, 2.182e7, 2.165e7, 2.155e7, 2.144e7)
wss.value$wss[11:20] <- c(2.614e7, 2.469e7, 2.370e7, 2.318e7, 2.289e7, 
                          2.268e7, 2.251e7, 2.234e7, 2.223e7, 2.213e7)

png("Figure_S2.png", height = 4, width = 6, units = "in", res = 300)
ggplot(data=wss.value, aes(x=k, y=wss)) +
  geom_line(color="black")+
  geom_point(color="black")+
  facet_wrap( ~ clustering_method, ncol = 2)+
  theme_bw()+
  xlab("number of cluster")+
  ylab("within-cluster sum of squares")+
  scale_x_continuous(breaks = seq(1,10,1))
dev.off()




