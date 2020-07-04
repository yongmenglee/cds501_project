##Import package
library(ggplot2)
library(animation)
library(e1071) #library for Cmean Fuzzy Clustering

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
set.seed(123)

##Fuzzy Clustering

##View optimal K (unscaled data set) Eulidean
cmean_withinerror_fuzzy <- function(k) {
  cluster <- cmeans(corona, k, iter.max = 100, dist = "euclidean")
  return (cluster$withinerror)
}
max_k <-10
wse <- sapply(2:max_k, cmean_withinerror_fuzzy)
elbowfuzzy <-data.frame(2:max_k, wse)
# Plot the graph with gglop
ggplot(elbowfuzzy, aes(x = X2.max_k, y = wse)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = seq(1, 20, by = 1))

##View optimal K (unscaled data set) Manhattan
cmean_withinerror_fuzzy <- function(k) {
  cluster <- cmeans(corona, k, iter.max = 100, dist = "manhattan")
  return (cluster$withinerror)
}
max_k <-10
wse <- sapply(2:max_k, cmean_withinerror_fuzzy)
elbowfuzzy <-data.frame(2:max_k, wse)
# Plot the graph with gglop
ggplot(elbowfuzzy, aes(x = X2.max_k, y = wse)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = seq(1, 20, by = 1))

##View optimal K (scaled data set) Eulidean
cmean_withinerror_fuzzy <- function(k) {
  cluster <- cmeans(corona_scale, k, iter.max = 100, dist = "euclidean")
  return (cluster$withinerror)
}
max_k <-10
wse <- sapply(2:max_k, cmean_withinerror_fuzzy)
elbowfuzzy <-data.frame(2:max_k, wse)
# Plot the graph with gglop
ggplot(elbowfuzzy, aes(x = X2.max_k, y = wse)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = seq(1, 20, by = 1))

##View optimal K (scaled data set) Manhattan
cmean_withinerror_fuzzy <- function(k) {
  cluster <- cmeans(corona_scale, k, iter.max = 100, dist = "manhattan")
  return (cluster$withinerror)
}
max_k <-10
wse <- sapply(2:max_k, cmean_withinerror_fuzzy)
elbowfuzzy <-data.frame(2:max_k, wse)
# Plot the graph with gglop
ggplot(elbowfuzzy, aes(x = X2.max_k, y = wse)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = seq(1, 20, by = 1))

#Fuzzy Clustering
corona_cluster_fuzzy <-cmeans(corona_scale, 5, dist = "euclidean")
corona_cluster_fuzzy$size
corona_cluster_fuzzy$cluster

clusterlabelscale <- as.data.frame(corona_cluster_fuzzy$cluster)
corona_label <- cbind(corona, clusterlabelscale)

cluster1 <- subset(corona_label, corona_cluster_fuzzy$cluster=='1')
cluster2 <- subset(corona_label, corona_cluster_fuzzy$cluster=='2')
cluster3 <- subset(corona_label, corona_cluster_fuzzy$cluster=='3')
cluster4 <- subset(corona_label, corona_cluster_fuzzy$cluster=='4')
cluster5 <- subset(corona_label, corona_cluster_fuzzy$cluster=='5')
as.data.frame(colMeans(cluster1[sapply(cluster1, is.numeric)]))
as.data.frame(colMeans(cluster2[sapply(cluster2, is.numeric)]))
as.data.frame(colMeans(cluster3[sapply(cluster3, is.numeric)]))
as.data.frame(colMeans(cluster4[sapply(cluster4, is.numeric)]))
as.data.frame(colMeans(cluster5[sapply(cluster5, is.numeric)]))


##PCA Features Extraction
princ <- prcomp(corona_scale)
nComp <- 2
corona_pca <- predict(princ, newdata=corona_scale)[,1:nComp]
as.data.frame(corona_pca)
corona_cluster_fuzzypca <-cmeans(corona_pca, 5, dist = "euclidean")
library(cluster)
clusplot(corona_pca, corona_cluster_fuzzypca$cluster, color=TRUE, shade=TRUE, labels=2, lines=0)
