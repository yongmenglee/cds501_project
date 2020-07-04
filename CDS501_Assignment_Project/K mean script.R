###K-means Clustering

##Import package
library(ggplot2)
library(animation)


##Data-preprocessing
#Import the dataset into R environment
corona_raw <- read.csv(file="MN997409.1-4NY0T82X016-Alignment-HitTable.csv", header=TRUE, sep=",") 
#Drop unnecessary information
corona <- corona_raw[-c(1,2,11)]
#re-scale the dataset
corona_scale <- scale(corona)
#view data summary information
summary(corona)
summary(corona_scale)


##PCA Features Extraction
nComp <- 2
pca <- predict(prcomp(corona), newdata=corona)[,1:nComp] 


##View optimal K (unscaled data set)
set.seed(123)
kmean_withinss <- function(k) {
  cluster <- kmeans(corona, k, iter.max = 100, nstart = 10)
  return (cluster$tot.withinss)
}
max_k <-10
wss <- sapply(2:max_k, kmean_withinss)
elbow <-data.frame(2:max_k, wss)

# Plot the graph with gglop
ggplot(elbow, aes(x = X2.max_k, y = wss)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = seq(1, 20, by = 1))


##View optimal K (scaled data set)
set.seed(123)
kmean_withinss <- function(k) {
  cluster <- kmeans(corona_scale, k, iter.max = 100, nstart = 10)
  return (cluster$tot.withinss)
}
max_k <-10
wss <- sapply(2:max_k, kmean_withinss)
elbow <-data.frame(2:max_k, wss)

# Plot the graph with gglop
ggplot(elbow, aes(x = X2.max_k, y = wss)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = seq(1, 20, by = 1))

#K-mean Clustering
corona_cluster_1 <-kmeans(corona, 5)
corona_cluster_1$cluster
corona_cluster_1$centers
corona_cluster_1$size	

kmeans.ani(corona, 5)






