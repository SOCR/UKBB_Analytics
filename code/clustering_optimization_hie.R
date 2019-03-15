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

############################# optimization by hierarchical ############################
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
# labs(subtitle = "Elbow method")+

#pdf("hierarchical_tuning_gap_statistics.pdf")
#set.seed(1234)
#fviz_nbclust(scaled.biomarker, hcut, nstart = 25,  method = "gap_stat", nboot = 500)+
  #labs(subtitle = "Gap statistic method")
#dev.off()
#####################################################################################
#####################################################################################
#####################################################################################
#nb_kmeans <- NbClust(scaled.biomarker, distance = "euclidean", min.nc = 2,
#              max.nc = 10, method = "kmeans")
#nb_hc.ward.D <- NbClust(scaled.biomarker, distance = "euclidean", min.nc = 2,
                     #max.nc = 5, method = "ward.D")
nb_hc.ward.D2 <- NbClust(scaled.biomarker, distance = "euclidean", min.nc = 2,
                 max.nc = 10, method = "ward.D2")

#capture.output(print(fviz_nbclust(nb_kmeans)), file = "output_nb_kmeans.txt")
capture.output(print(fviz_nbclust(nb_hc.ward.D2)), file = "output_nb_hc.ward.D2.txt")

#fviz_nbclust(nb_kmeans)
#fviz_nbclust(nb_hc.ward.D)
#fviz_nbclust(nb_hc.ward.D2)

#####################################################################################
#####################################################################################
#####################################################################################
## make the optimization plot by ggplot2
# sil.value <- data.frame(matrix(0,ncol = 3, nrow = 20))
# colnames(sil.value) <- c("k", "silhouette_value", "clustering_method")
# sil.value$k <- rep(seq(1,10,by=1), times=2)
# sil.value$clustering_method <- rep(c("k_means", "hierarchical"), each=10)
# sil.value$clustering_method <- factor(sil.value$clustering_method, levels = c("k_means", "hierarchical"))
# sil.value$silhouette_value[1:10] <- c(0, 0.074, 0.068, 0.031, 0.030, 0.024,
#                                       0.018, 0.017, 0.018, 0.013)
# sil.value$silhouette_value[11:20] <- c(0, 0.057, 0.043, 0.038, 0.024, 0.014,
#                                        0.003, 0.004, -0.003, -0.006)
# 
# ggplot(data=sil.value, aes(x=k, y=silhouette_value)) +
#   geom_line(color="steelblue")+
#   geom_point(color="steelblue")+
#   facet_wrap( ~ clustering_method, ncol = 2)+
#   theme_bw()+
#   xlab("number of cluster k")+
#   ylab("average silhouette value")


# sil.value.kmeans <- data.frame(matrix(0, ncol=2, nrow=10))
# colnames(sil.value.kmeans) <- c("k", "silhouette_value")
# sil.value.kmeans$k <- seq(1,10,by=1)
# sil.value.kmeans$silhouette_value <- c(0, 0.074, 0.068, 0.031, 0.030, 0.024,
#                                       0.018, 0.017, 0.018, 0.013)
# 
# 
# sil.value.hierarchical <- data.frame(matrix(0, ncol=2, nrow=10))
# colnames(sil.value.hierarchical) <- c("k", "silhouette_value")
# sil.value.hierarchical$k <- seq(1,10,by=1)
# sil.value.hierarchical$silhouette_value <- c(0, 0.057, 0.043, 0.038, 0.024, 0.014,
#                                              0.003, 0.004, -0.003, -0.006)
# 
# 
# p1 <- ggplot(data=sil.value.kmeans, aes(x=k, y=silhouette_value))+
#   geom_line(color="steelblue")+
#   geom_point(color="steelblue")+
#   theme_bw()+
#   xlab("number of cluster k")+
#   ylab("average silhouette value")+
#   ggtitle("k-means clustering")+
#   ylim(-0.01,0.08)+
#   theme(plot.title = element_text(size=10),
#         axis.text.x = element_text(size=8),
#         axis.text.y = element_text(size=8))+
#   scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10))
# 
# p2 <- ggplot(data=sil.value.hierarchical, aes(x=k, y=silhouette_value))+
#   geom_line(color="steelblue")+
#   geom_point(color="steelblue")+
#   theme_bw()+
#   xlab("number of cluster k")+
#   ylab("average silhouette value")+
#   ggtitle("hierarchical clustering")+
#   ylim(-0.01,0.08)+
#   theme(plot.title = element_text(size=10),
#         axis.text.x = element_text(size=8),
#         axis.text.y = element_text(size=8))+
#   scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10))
# 
# plot_grid(p1, p2, labels = c("a", "b"))
# 
# 
# #####################################################################################
# #####################################################################################
# wss.value.kmeans <- data.frame(matrix(0, ncol=2, nrow=10))
# colnames(sil.value.kmeans) <- c("k", "wss")
# wss.value.kmeans$k <- seq(1,10,by=1)
# wss.value.kmeans$wss <- c(3.265e7, 3.024e7, 2.892e7, 2.833e7, 2.788e7,
#                           2.770e7, 2.730e7, 2.709e7, 2.693e7, 2.675e7)
# 
# 
# wss.value.hierarchical <- data.frame(matrix(0, ncol=2, nrow=10))
# colnames(wss.value.hierarchical) <- c("k", "wss")
# wss.value.hierarchical$k <- seq(1,10,by=1)
# wss.value.hierarchical$wss <- c(3.266e7, 3.072e7, 2.957e7, 2.902e7, 2.872e7, 
#                                 2.846e7, 2.819e7, 2.798e7, 2.783e7, 2.770e7)
# 
# 
# p1 <- ggplot(data=wss.value.kmeans, aes(x=k, y=wss))+
#   geom_line(color="steelblue")+
#   geom_point(color="steelblue")+
#   theme_bw()+
#   xlab("number of cluster k")+
#   ylab("total within sum of square")+
#   ggtitle("k-means clustering")+
#   ylim(2.6e7,3.3e7)+
#   theme(plot.title = element_text(size=10),
#         axis.text.x = element_text(size=8),
#         axis.text.y = element_text(size=8))+
#   scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10))
# 
# p2 <- ggplot(data=wss.value.hierarchical, aes(x=k, y=wss))+
#   geom_line(color="steelblue")+
#   geom_point(color="steelblue")+
#   theme_bw()+
#   xlab("number of cluster k")+
#   ylab("total within sum of square")+
#   ggtitle("hierarchical clustering")+
#   ylim(2.6e7,3.3e7)+
#   theme(plot.title = element_text(size=10),
#         axis.text.x = element_text(size=8),
#         axis.text.y = element_text(size=8))+
#   scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10))
# 
# plot_grid(p1, p2, labels = c("a", "b"))
# 
