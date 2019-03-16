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
#####################################################################################
#####################################################################################
#####################################################################################
nb_hc.ward.D2 <- NbClust(scaled.biomarker, distance = "euclidean", min.nc = 2,
                 max.nc = 10, method = "ward.D2")

capture.output(print(fviz_nbclust(nb_hc.ward.D2)), file = "output_nb_hc.ward.D2.txt")
