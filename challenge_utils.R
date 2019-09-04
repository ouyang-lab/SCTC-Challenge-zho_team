library(NbClust)

esBinDist <- function(top10,MCC, geometry){
  maxBins <- apply(MCC, 2, which.max)
  maxBins <- geometry[maxBins,]
  
  ds <- rep(0, nrow(top10))
  for(i in 1:nrow(top10)){  
    rs = proxy::dist(geometry[top10[i,],], maxBins[i,])
    ds[i] <- mean(rs)
  }
  #plot(ds)
  print(mean(ds))
  ds
}  

FindGene <- function(cluster, mat){
  # mat: original data matrix, cell * gene
  output <- NULL
  genename <- colnames(mat)
  for (i in 1:max(cluster)){
    ## find all the elements belong to this cluster
    tmp.mat <- as.matrix(mat[, cluster == i])
    tmp.name <- genename[cluster == i]
    tmp.center <- as.matrix(rowMeans(tmp.mat))
    ## find the nearest pos to this center
    tmp.loc <- which.min(colSums((tmp.mat - 
                                    tmp.center %*% 
                                    t(rep(1, ncol(tmp.mat)))) ^ 2))
    output <- c(output, tmp.name[tmp.loc])
  }
  return(output)
}

GetGene_Hc <- function(scRNA, method, k){
  distmat <- dist(t(scRNA))
  hc <- hclust(distmat, method = method)
  clu <- cutree(hc, k = k)
  return(FindGene(clu, scRNA))
}

calCor <- function(dm_sel){  
  len1 = dim(dm_sel@insitu.matrix)[1]
  len2 = dim(dm_sel@binarized.data)[2]
  
  mcor <- matrix(0, len1, len2)
  for(i in 1:len1){
    vc = dm_sel@insitu.matrix[i,]
    mcor[i,] <- apply(dm_sel@binarized.data, 2, function(x) { 
      a = abs(x-vc)
      a[a>0] = 1
      1 - sum(a)/length(vc) 
    })
  }
  #mcor <- cor(t(dm_sel@insitu.matrix),  dm_sel@binarized.data)
  dm_sel@mcc.scores <-  mcor 
  dm_sel
}

getTop10 <- function(mcc, geometry, th=0.74){
  top10 <- matrix(0, nrow = ncol(mcc), ncol = 10)
  ks <- rep(0, ncol(mcc))
  topcor <- list()
  for(i in 1:ncol(mcc)){
    rc <- mcc[,i]
    od <- order(rc, decreasing = TRUE)
    
    k = length(which(rc > th)) # 
    if(k==0){
      k = 1
    }else if(k>100){
      k=100
    } 
    topd <- geometry[od[1:k],]
    topc <- rc[od[1:k]] 
    topcor[[i]] <- topc
    
    if(k > 3){ 
      
      if(k<10){
        res<-NbClust(data =topd , distance = "euclidean", min.nc= 2, max.nc=k-2, 
                     method = "centroid", index = "silhouette")
      }else{
        res<-NbClust(data =topd , distance = "euclidean", min.nc= 2, max.nc=8, 
                     method = "centroid", index = "silhouette")
      }
      
      cluid <- res$Best.partition
      table(cluid)
      
      mrc <- rep(0, length(unique(cluid)))
      src <- rep(0, length(unique(cluid)))
      for(j in 1:length(mrc)){
        mrc[j] <- mean(topc[which(cluid==j)])
        src[j] <- sum(topc[which(cluid==j)])
      }
      mid <- which(src==max(src))
      if(length(mid)>1){
        kid <- which.max(mrc[mid])
        mid <- mid[kid]
      } 
      
      binid <- as.numeric(names(cluid[which(cluid==mid)]))
      center <-  as.matrix(colMeans(geometry[binid,]))
      neards <- proxy::dist(  t(center),geometry[1:3039,])
    }else{
      mid <- od[1]
      center <- as.matrix(geometry[mid,])
      neards <- proxy::dist(center, geometry[1:3039,])
    }  
    
    nod <- order(neards[1,])
    top10[i,] <- nod[1:10]   
  }   
  list(top10 =top10,   topcor=topcor)
}

