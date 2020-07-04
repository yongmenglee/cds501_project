# CDS501 Project Final Report
# Group 3 - Coronavirus

# Read data from CSV file
corona <- read.table('MN997409.1-4NY0T82X016-Alignment-HitTable.csv', header = TRUE, sep=',')
summary(corona)

# Remove columns with only one unique value:-
corona <- corona[-c(1, 11)]
summary(corona)



# Hierarchical clustering



# === Data Preparation for Hierarchical Clustering
# ---
# Define columns to be used for hierarchical clustering
vars.to.use <- colnames(corona)[-1]   # now we only have 9 columns

# Transform data frame into a scaled matrix
corona_scale <- scale(corona[,vars.to.use])


# Create distance matrix
d_euclidean <- dist(corona_scale, method = "euclidean")



# === Perform hierarchical clustering + dendrogram
# ---
library(ggplot2)

# Helper function #1: Print selected columns of the data for each cluster
print_clusters <- function(labels, k) { 
  for(i in 1:k) { 
    print(paste("cluster", i)) 
    print(corona[labels==i, c("X.subject.acc.ver.", 
                              "X.bit.score.", 
                              "X.alignment.length.", 
                              "X.mismatches.", 
                              "X.gap.opens.")]) 
  }
}

# Helper function #2: Generate dendrogram
plot_hclust <- function(pfitHClust, k = NA, label = FALSE) {
  plot(pfitHClust, labels = label)
  if (is.na(k) == FALSE) {
    rect.hclust(pfitHClust, k = k)
  }
}

# ---
# Different clustering methods
plot_hclust(hclust(d_euclidean, method = "ward.D"))
plot_hclust(hclust(d_euclidean, method = "ward.D2"))
plot_hclust(hclust(d_euclidean, method = "complete"))
plot_hclust(hclust(d_euclidean, method = "average"))

# Different k using ward.D
plot_hclust(hclust(d_euclidean, method = "ward.D"), k = 4)
plot_hclust(hclust(d_euclidean, method = "ward.D"), k = 5)



# === Running clusterboot
# ---
# Import fpc library for clusterboot() function
library(fpc)

set.seed(123)

# Helper function to display the result of clusterboot
display.cboot.result <- function(cboot.result) {
  summary(cboot.result)
  
  groups <- cboot.result$result$partition
  print_clusters(groups, kbest.p)
  
  print("Vector of cluster stabilities")
  print(cboot.result$bootmean)
  
  print("How many times each cluster was dissolved?")
  print(cboot.result$bootbrd)
}

# 4 clusters
kbest.p <- 4

# Note: invisible(capture.output()) is used to suppress output
# ward.D
invisible(capture.output(cboot.hclust <- clusterboot(corona_scale, 
                                                     clustermethod = hclustCBI, 
                                                     method = "ward.D", 
                                                     k = kbest.p) ))
display.cboot.result(cboot.hclust)

# ward.D2
invisible(capture.output(cboot.hclust <- clusterboot(corona_scale, 
                                                     clustermethod = hclustCBI, 
                                                     method = "ward.D2", 
                                                     k = kbest.p) ))
display.cboot.result(cboot.hclust)

# complete
invisible(capture.output(cboot.hclust <- clusterboot(corona_scale, 
                                                     clustermethod = hclustCBI, 
                                                     method = "complete", 
                                                     k = kbest.p) ))
display.cboot.result(cboot.hclust)

# average
invisible(capture.output(cboot.hclust <- clusterboot(corona_scale, 
                                                     clustermethod = hclustCBI, 
                                                     method = "average", 
                                                     k = kbest.p) ))
display.cboot.result(cboot.hclust)


# 5 clusters
kbest.p <- 5

# Note: invisible(capture.output()) is used to suppress output

# ward.D
invisible(capture.output(cboot.hclust <- clusterboot(corona_scale, 
                                                     clustermethod = hclustCBI, 
                                                     method = "ward.D", 
                                                     k = kbest.p) ))
display.cboot.result(cboot.hclust)

# ward.D2
invisible(capture.output(cboot.hclust <- clusterboot(corona_scale, 
                                                     clustermethod = hclustCBI, 
                                                     method = "ward.D2", 
                                                     k = kbest.p) ))
display.cboot.result(cboot.hclust)

# complete
invisible(capture.output(cboot.hclust <- clusterboot(corona_scale, 
                                                     clustermethod = hclustCBI, 
                                                     method = "complete", 
                                                     k = kbest.p) ))
display.cboot.result(cboot.hclust)

# average
invisible(capture.output(cboot.hclust <- clusterboot(corona_scale, 
                                                     clustermethod = hclustCBI, 
                                                     method = "average", 
                                                     k = kbest.p) ))
display.cboot.result(cboot.hclust)



# === Optimal K
# ---
# Package for melt() function (in helper function #6)
library(reshape2)

# Helper function #1: Calculate Euclidean distance between two vectors
sqr_edist <- function(x, y) { sum((x-y)^2) }

# Helper function #2: Calculate WSS of a cluster
wss.cluster <- function(clustermat) {
  # Calculate the centroid of the cluster
  c0 <- apply(clustermat, 2, FUN=mean)   
  # Calculate the sum of Euclidean distances between pairs of data points
  # and the cluster, c0
  sum(apply(clustermat, 1, FUN = function(row) { sqr_edist(row,c0) } ))  
}

# Helper function #3: Calculate total WSS for all clusters
wss.total <- function(dmatrix, labels) {                          
  wsstot <- 0    
  k <- length(unique(labels))
  
  for(i in 1:k)
    # Extract each cluster, calculate the cluster's WSS
    # and add up all the values
    wsstot <- wsstot + wss.cluster(subset(dmatrix, labels == i))
  
  wsstot
}

# Helper function #4: Calculate total sum of squares
totss <- function(dmatrix) {                       
  grandmean <- apply(dmatrix, 2, FUN=mean)    
  sum(apply(dmatrix, 1, FUN = function(row) { sqr_edist(row, grandmean) } ))
}

# Helper function #5: Calculate CH criterion for "kmeans" or "hclust".
ch_criterion <- function(dmatrix, kmax, method = "kmeans", ch.method = "ward.D2") {           
  if (!(method %in% c("kmeans", "hclust"))) {
    stop("method must be one of c('kmeans', 'hclust')")
  }
  
  npts <- dim(dmatrix)[1]  # number of rows.      
  
  # Total sum of squares is independent of the clustering
  totss <- totss(dmatrix)
  
  wss <- numeric(kmax)    
  crit <- numeric(kmax)    
  
  # Calculate WSS for k = 1 (which is the total sums of squares)
  wss[1] <- (npts - 1) * sum(apply(dmatrix, 2, var))
  
  # Calculate WSS for k from 2 to kmax
  for (k in 2:kmax) {
    if (method=="kmeans") {
      # kmeans() returns the total WSS as one of its outputs
      clustering <- kmeans(dmatrix, k, nstart = 10, iter.max = 100)        
      wss[k] <- clustering$tot.withinss
    }
    else {
      # for hclust(), calculate total wss by hand
      d <- dist(dmatrix, method = "euclidean")        
      pfit <- hclust(d, method = ch.method)
      # pfit <- hclust(d, method = "single")
      labels <- cutree(pfit, k = k)        
      wss[k] <- wss.total(dmatrix, labels)      
    }    
  }
  
  # Calculate BSS for k from 1 to kmax
  bss <- totss - wss   
  
  # NOrmalize BSS by k-1
  crit.num <- bss/(0:(kmax-1))
  
  # Normalize WSS by npts-k
  crit.denom <- wss/(npts - 1:kmax)   
  
  # Return a vector of CH indices and of WSS for k from 1 to kmax
  # Also return total sum of squares
  list(crit = crit.num/crit.denom, wss = wss, totss = totss)
}

# Helper function #6: Plot both CH and WSS - find optimal k
plot_ch_criterion <- function(pmat, kmax = 10, method = "hclust", hClustMethod = "ward.D2") {
  clustcrit <- ch_criterion(pmat, kmax, method = method, ch.method = hClustMethod)
  
  print(clustcrit$crit)
  print(clustcrit$wss)
  
  # Create a data frame with number of clusters, the scaled CH criterion 
  # and the scaled WSS criterion (enable plotting on the same graph)
  critframe <- data.frame(k = 1:kmax, 
                          ch = scale(clustcrit$crit), 
                          wss = scale(clustcrit$wss))
  
  # melt() function: put the data frame in a shape suitable for ggplot
  critframe <- melt(critframe, id.vars = c("k"), 
                    variable.name="measure", 
                    value.name="score")
  
  plottitle <- paste0("CH vs WSS: ", hClustMethod)
  print(plottitle)
  
  # Plot both CH and WSS
  ggplot(critframe, aes(x = k, y = score, color = measure)) + 
    geom_point(aes(shape = measure)) + 
    geom_line(aes(linetype = measure)) + 
    scale_x_continuous(breaks = 1:kmax, labels = 1:kmax) +
    ggtitle(plottitle)
}

# Evaluate number of clusters
# ---
plot_ch_criterion(corona_scale, kmax = 10, hClustMethod = "ward.D")
plot_ch_criterion(corona_scale, kmax = 10, hClustMethod = "ward.D2")
plot_ch_criterion(corona_scale, kmax = 10, hClustMethod = "complete")
plot_ch_criterion(corona_scale, kmax = 10, hClustMethod = "average")



# === Plot distribution of using PCA
# ---
# Helper function: Visualize principal components
plot_hclust_pc <- function(pmat, nComp, groups, subject, method.name = NA) {
  princ <- prcomp(pmat)
  
  project <- predict(princ, newdata = pmat)[, 1:nComp]
  
  project.plus <- cbind(as.data.frame(project), 
                        cluster = as.factor(groups), 
                        subject = subject)
  
  plottitle <- ifelse(!is.na(method.name), paste0("Clusters: ", method.name), "Clusters")
  
  plotobj <- ggplot(project.plus, aes(x = PC1, y = PC2, color = cluster)) + 
    geom_point(aes(shape = cluster)) + 
    geom_text(aes(label = subject), 
              hjust = 0, 
              vjust = 1, 
              check_overlap = TRUE,
              show.legend = FALSE) +
    ggtitle(plottitle)
  
  list(plotobj = plotobj, project = project, princ = princ)
}


# k = 4
num.k <- 4
pMethod <- "ward.D"

pfit2 <- hclust(d_euclidean, method = pMethod)
groups <- cutree(pfit2, k = num.k)

# Print list of virus for each cluster
print_clusters(groups, num.k)

# Plot data with respect to each cluster
corona_orig <- corona[,vars.to.use]

# unscaled
pc <- plot_hclust_pc(corona_orig, 2, groups, corona$X.subject.acc.ver., method.name = pMethod)
pc$plotobj

# scaled
pc <- plot_hclust_pc(corona_scale, 2, groups, corona$X.subject.acc.ver., method.name = pMethod)
pc$plotobj


# k = 5
pMethod <- "ward.D"

pfit2 <- hclust(d_euclidean, method = pMethod)
groups <- cutree(pfit2, k = 5)

# Print list of virus for each cluster
print_clusters(groups, num.k)

# Plot data with respect to each cluster
corona_orig <- corona[,vars.to.use]

pc <- plot_hclust_pc(corona_orig, 2, groups, corona$X.subject.acc.ver., method.name = pMethod)
pc$plotobj

# scaled
pc <- plot_hclust_pc(corona_scale, 2, groups, corona$X.subject.acc.ver., method.name = pMethod)
pc$plotobj



# project <- pc$project
# princ <- pc$princ
# 
# # install.packages("factoextra")
# library(factoextra)
# 
# # Visualization: Clustering
# fviz_cluster(list(data = project, cluster = groups))
# 
# # Visualization (PCA): Show the percentage of variances explained by
# # each principal component.
# fviz_eig(princ)
# 
# # Visualization (PCA): Graph of individuals
# fviz_pca_ind(princ,
#              col.ind = "cos2", # Color by the quality of representation
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              # repel = TRUE,     # Avoid text overlapping
# )
# 
# # 
# # # Visualization (PCA): graph of variables
# # # - Highly (positive) correlated variables point to the same direction
# # # - Highly (negative) correlated variables point to the opposite direction
# # fviz_pca_var(princ,
# #              col.var = "contrib", # Color by contributions to the PC
# #              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
# #              repel = TRUE     # Avoid text overlapping
# # )
# # 
# # # WSS
# # fviz_nbclust(pmatrix, FUN = hcut, k.max = 20, method = "wss")
# # 
# # # CH criterion
# # fviz_nbclust(pmatrix, FUN = hcut, k.max = 20, method = "silhouette")











# Testing

library(FactoMineR)


HCPC(corona_scale, nb.clust = 4, min = 3, max = NULL, graph = TRUE)

# Compute PCA with ncp = 3
res.pca <- PCA(corona_scale, ncp = 4, graph = FALSE)
# Compute hierarchical clustering on principal components
res.hcpc <- HCPC(res.pca, graph = FALSE)

fviz_dend(res.hcpc, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
)

fviz_cluster(res.hcpc,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
)

# Principal components + tree
plot(res.hcpc, choice = "3D.map")

res.hcpc$desc.var$quanti

# #K-means
# pclusters <- kmeans(pmatrix, kbest.p, nstart=100, iter.max=100)
# summary(pclusters)
# pclusters$centers
# pclusters$size
# groups <- pclusters$cluster
# print_clusters(groups, kbest.p)
# clustering.ch <- kmeansruns(pmatrix, krange=1:10, criterion="ch") 
# clustering.ch$bestk 
# clustering.asw <- kmeansruns(pmatrix, krange=1:10, criterion="asw")
# clustering.asw$bestk
# clustering.ch$crit
# clustcrit$crit 
# critframe <- data.frame(k=1:10, ch=scale(clustering.ch$crit), asw=scale(clustering.asw$crit)) 
# critframe <- melt(critframe, id.vars=c("k"), variable.name="measure", value.name="score") 
# ggplot(critframe, aes(x=k, y=score, color=measure)) + geom_point(aes(shape=measure)) + geom_line(aes(linetype=measure)) + scale_x_continuous(breaks=1:10, labels=1:10)
# 
# 
# #Running clusterboot() with k-means 
# summary(clustering.ch)
# kbest.p<-5 
# cboot<-clusterboot(pmatrix, clustermethod=kmeansCBI,              
#                    runs=100,iter.max=100,              
#                    krange=kbest.p, seed=15555)
# groups <- cboot$result$partition 
# print_clusters(cboot$result$partition, kbest.p)
# cboot$bootmean
# cboot$bootbrd