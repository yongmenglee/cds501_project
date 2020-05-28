#Dataset import and drop column
corona <- read.csv(file="MN997409.1-4NY0T82X016-Alignment-HitTable.csv", header=TRUE, sep=",") 
corona <- corona[ -c(1,11) ]


#Hierarchical clustering
#Hierarchical clustering (Pre-process)
vars.to.use <- colnames(corona)[-1] 
pmatrix <- scale(corona[,vars.to.use])
pcenter <- attr(pmatrix, "scaled:center") 
pscale <- attr(pmatrix, "scaled:scale") 
d <- dist(pmatrix, method="euclidean") 
pfit <- hclust(d, method="ward.D2") 
plot(pfit, labels=corona$X.subject.acc.ver.) 
rect.hclust(pfit, k=5) 

#Hierarchical clustering
print_clusters <- function(labels, k) { 
  for(i in 1:k) { 
  print(paste("cluster", i)) 
  print(corona[labels==i,c("X.subject.acc.ver.","X.bit.score.","X.alignment.length.","X.mismatches.","X.gap.opens.")]) 
  }}
print_clusters(groups, 5) 

#Visualizing Hierarchical clustering
library(ggplot2) 
princ <- prcomp(pmatrix) 
nComp <- 2
project <- predict(princ, newdata=pmatrix)[,1:nComp] 
project.plus <- cbind(as.data.frame(project), cluster=as.factor(groups), subject=corona$X.subject.acc.ver.)
ggplot(project.plus, aes(x=PC1, y=PC2)) + geom_point(aes(shape=cluster)) + geom_text(aes(label=country), hjust=0, vjust=1)
ggplot(project.plus, aes(x=PC1, y=PC2)) + geom_point(aes(shape=cluster)) + geom_text(aes(label=subject), hjust=0, vjust=1)

#Running clusterboot
library(fpc)
kbest.p<-5 
cboot.hclust <- clusterboot(pmatrix,clustermethod=hclustCBI, method="ward.D2", k=kbest.p) 
summary(cboot.hclust$result) 
groups<-cboot.hclust$result$partition 
print_clusters(groups, kbest.p) 
cboot.hclust$bootmean 
cboot.hclust$bootbrd

#Optimal K
sqr_edist <- function(x, y) {sum((x-y)^2)}
wss.cluster <- function(clustermat) {           
  c0 <- apply(clustermat, 2, FUN=mean)    
  sum(apply(clustermat, 1, FUN=function(row){sqr_edist(row,c0)}))  
  }
wss.total <- function(dmatrix, labels) {                          
  wsstot <- 0    
  k <- length(unique(labels))    
  for(i in 1:k)      
    wsstot <- wsstot + wss.cluster(subset(dmatrix, labels==i)) 
  wsstot
}
totss <- function(dmatrix) {                       
  grandmean <- apply(dmatrix, 2, FUN=mean)    
  sum(apply(dmatrix, 1, FUN=function(row){sqr_edist(row, grandmean)}))} 
ch_criterion <- function(dmatrix, kmax, method="kmeans") {           
  if(!(method %in% c("kmeans", "hclust"))) {
    stop("method must be one of c('kmeans', 'hclust')")}    
  npts <- dim(dmatrix)[1]  # number of rows.      
  totss <- totss(dmatrix)    
  wss <- numeric(kmax)    
  crit <- numeric(kmax)    
  wss[1] <- (npts-1)*sum(apply(dmatrix, 2, var))   
  for(k in 2:kmax) {     
    if(method=="kmeans") {        
      clustering<-kmeans(dmatrix, k, nstart=10, iter.max=100)        
      wss[k] <- clustering$tot.withinss      
    }else {
      d <- dist(dmatrix, method="euclidean")        
      pfit <- hclust(d, method="ward.D")        
      labels <- cutree(pfit, k=k)        
      wss[k] <- wss.total(dmatrix, labels)      
    }    
  }    
  bss <- totss - wss   
  crit.num <- bss/(0:(kmax-1))   
  crit.denom <- wss/(npts - 1:kmax)   
  list(crit = crit.num/crit.denom, wss = wss, totss = totss)
} 
library(reshape2)
clustcrit <- ch_criterion(pmatrix, 10, method="hclust")
critframe <- data.frame(k=1:10, ch=scale(clustcrit$crit), wss=scale(clustcrit$wss))
critframe <- melt(critframe, id.vars=c("k"), variable.name="measure", value.name="score")
ggplot(critframe, aes(x=k, y=score, color=measure)) + geom_point(aes(shape=measure)) + geom_line(aes(linetype=measure)) + scale_x_continuous(breaks=1:10, labels=1:10) 


#K-means
pclusters <- kmeans(pmatrix, kbest.p, nstart=100, iter.max=100)
summary(pclusters)
pclusters$centers
pclusters$size
groups <- pclusters$cluster
print_clusters(groups, kbest.p)
clustering.ch <- kmeansruns(pmatrix, krange=1:10, criterion="ch") 
clustering.ch$bestk 
clustering.asw <- kmeansruns(pmatrix, krange=1:10, criterion="asw")
clustering.asw$bestk
clustering.ch$crit
clustcrit$crit 
critframe <- data.frame(k=1:10, ch=scale(clustering.ch$crit), asw=scale(clustering.asw$crit)) 
critframe <- melt(critframe, id.vars=c("k"), variable.name="measure", value.name="score") 
ggplot(critframe, aes(x=k, y=score, color=measure)) + geom_point(aes(shape=measure)) + geom_line(aes(linetype=measure)) + scale_x_continuous(breaks=1:10, labels=1:10)


#Running clusterboot() with k-means 
summary(clustering.ch)
kbest.p<-5 
cboot<-clusterboot(pmatrix, clustermethod=kmeansCBI,              
                   runs=100,iter.max=100,              
                   krange=kbest.p, seed=15555)
groups <- cboot$result$partition 
print_clusters(cboot$result$partition, kbest.p)
cboot$bootmean
cboot$bootbrd

# Title: A function to assign points to a cluster 
assign_cluster <- function(newpt, centers, xcenter=0, xscale=1) {
  xpt <- (newpt - xcenter)/xscale
  dists <- apply(centers, 1, FUN=function(c0){sqr_edist(c0, xpt)})
  which.min(dists)
}

# Title: An example of assigning points to cluster 
rnorm.multidim <- function(n, mean, sd, colstr="x") {           
  ndim <- length(mean)     
  data <- NULL     
  for(i in 1:ndim) {       
    col <- rnorm(n, mean=mean[[i]], sd=sd[[i]])       
    data<-cbind(data, col)     
  }     
  cnames <- paste(colstr, 1:ndim, sep='')     
  colnames(data) <- cnames     
  data   
}
mean1 <- c(1, 1, 1)                       
sd1 <- c(1, 2, 1)
mean2 <- c(10, -3, 5) 
sd2 <- c(2, 1, 2)
mean3 <- c(-5, -5, -5) 
sd3 <- c(1.5, 2, 1)
clust1 <- rnorm.multidim(100, mean1, sd1)
clust2 <- rnorm.multidim(100, mean2, sd2)
clust3 <- rnorm.multidim(100, mean3, sd3)
toydata <- rbind(clust3, rbind(clust1, clust2))
tmatrix <- scale(toydata)
tcenter <- attr(tmatrix, "scaled:center")
tscale<-attr(tmatrix, "scaled:scale")
kbest.t <- 3
tclusters <- kmeans(tmatrix, kbest.t, nstart=100, iter.max=100) 
tclusters$size
unscale <- function(scaledpt, centervec, scalevec) {
    scaledpt*scalevec + centervec
}
unscale(tclusters$centers[1,], tcenter, tscale)
mean2
unscale(tclusters$centers[2,], tcenter, tscale)
mean3
unscale(tclusters$centers[3,], tcenter, tscale) 
mean1
assign_cluster(rnorm.multidim(1, mean1, sd1),
               tclusters$centers, 
               tcenter, tscale)
assign_cluster(rnorm.multidim(1, mean2, sd1),      
               tclusters$centers, 
               tcenter, tscale)
assign_cluster(rnorm.multidim(1, mean3, sd1),        
               tclusters$centers, 
               tcenter, tscale) 
