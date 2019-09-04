source("newUtils.R") 
library(DistMap)

dim(normalized.data) 

folds_train <- read.table("/Users/cheny/Documents/scRNA3D/postCh/folds_train.csv", sep = "\t", 
                          header= FALSE, row.names = NULL, stringsAsFactors = F,
                          quote = "")
train.data <- list()
for(i in 1:10){
  train.data[[i]] <- as.numeric(strsplit(folds_train[i,1], split = ",")[[1]])
}

folds_train = 0 

folds_test <- read.table("/Users/cheny/Documents/scRNA3D/postCh/folds_test.csv", sep = "\t", 
                         header= FALSE, row.names = NULL, stringsAsFactors = F,
                         quote = "")
test.data <- list()
for(i in 1:10){
  test.data[[i]] <- as.numeric(strsplit(folds_test[i,1], split = ",")[[1]])
}

sapply(test.data, length)
 
geneAll <- matrix(0, 10, 60)

binDist60 = list()
geneCor60 = list()
geneCor60Know = list()
for(i in 1:10){
  #####1.Training: Use the training data to select genes 
  insitu.genes <- colnames(insitu.matrix)
  scRNA <- t(normalized.data[insitu.genes, train.data[[i]] ])
  selGene <- GetGene_Hc(scRNA, "mcquitty", 60)# Hclust-mcquitty
  geneAll[i,] <- selGene
  
  dm = new("DistMap",
           raw.data = raw.data,
           data=normalized.data[,train.data[[i]]],
           insitu.matrix=insitu.matrix[,selGene],
           geometry=as.matrix(geometry)) 
  qth <- getBestBinTh(dm, quantiles=seq(0.15, 0.5, 0.01)) 
  
  #####2.Prediction: use the selected genes and prediction data to predict the cell position  
  
  #####2.1  
  dm = new("DistMap",
           raw.data = raw.data,
           data=normalized.data[,test.data[[i]]],
           insitu.matrix=insitu.matrix[,selGene],
           geometry=as.matrix(geometry))
  #####2.2 Binarization
  dm <- binFun(dm, quantiles=qth)
  #####2.3 Calculate MCC: 3039bins * test cells
  dm <- calCor(dm) 
  maxc60 <- apply(dm@mcc.scores, 2, max) 
  v = summary(maxc60)
  #####2.4 Prediction top 10 bins for the test data
  bins <- getTop10(dm@mcc.scores, geometry, v[5]) 
  dim(bins$top10)
  
  if(1==0){
    geneM  <- matrix(selGene, nrow=6,ncol=10)
    res <- rbind(geneM ,bins$top10)
    rownames(res) <- c(rep("",6),  test.data[[i]] )
    fname = paste("post_challenge_zho_team/60genes", i, ".csv", sep = "_")
    write.table(res, file = fname,quote = FALSE, sep=",", 
                col.names=FALSE, row.names = TRUE) 
    
  }else if(1==2){
    tid = test.data[[i]]
    binDist60[[i]] = esBinDistAll(bins$top10, MCC[,tid], geometry)
    
    coff = getCor(dm, insitu.matrix, insitu.genes)  
    geneCor60[[i]] = coff$corM
    
    coff = getCorKnow(dm, insitu.matrix, insitu.genes) 
    geneCor60Know[[i]] = coff$corM
  }else{
    tid = test.data[[i]]
    binDist60[[i]] = scorePost(selGene,tid, bins$top10, ground.truth[tid,], ambig.locations, rawdm)
  } 
} 
 

binDist40 = list()
geneCor40 = list()
geneCor40Know = list()
gnum = 40 
geneAll <- matrix(0, 10, gnum)
for(i in 1:10){
  #####1.Training: Use the training data to select genes 
  insitu.genes <- colnames(insitu.matrix)
  scRNA <- t(normalized.data[insitu.genes, train.data[[i]] ])
  selGene <- GetGene_Hc(scRNA, "mcquitty", gnum)# Hclust-mcquitty
  geneAll[i,] <- selGene
  
  dm = new("DistMap",
           raw.data = raw.data,
           data=normalized.data[,train.data[[i]]],
           insitu.matrix=insitu.matrix[,selGene],
           geometry=as.matrix(geometry)) 
  qth <- getBestBinTh(dm, quantiles=seq(0.15, 0.5, 0.01)) 
  
  #####2.Prediction: use the selected genes and prediction data to predict the cell position  
  
  #####2.1  
  dm = new("DistMap",
           raw.data = raw.data,
           data=normalized.data[,test.data[[i]]],
           insitu.matrix=insitu.matrix[,selGene],
           geometry=as.matrix(geometry))
  #####2.2 Binarization
  dm <- binFun(dm, quantiles=qth)
  #####2.3 Calculate MCC: 3039bins * test cells
  dm <- calCor(dm) 
  maxc40 <- apply(dm@mcc.scores, 2, max) 
  v = summary(maxc40)
  #####2.4 Prediction top 10 bins for the test data
  bins <- getTop10(dm@mcc.scores, geometry, v[2])  
  
  if(1==0){
    geneM  <- matrix(selGene, nrow=4,ncol=10)
    res <- rbind(geneM,bins$top10)
    rownames(res) <- c(rep("",4),  test.data[[i]] )
    fname = paste("post_challenge_zho_team/40genes", i, ".csv", sep = "_")
    write.table(res, file = fname,quote = FALSE, sep=",", 
                col.names=FALSE, row.names = TRUE) 
  }else if(1==2){
    tid = test.data[[i]]
    binDist40[[i]] = esBinDistAll(bins$top10, MCC[,tid], geometry)
    
    coff = getCor(dm, insitu.matrix, insitu.genes)  
    geneCor40[[i]] = coff$corM
    
    coff = getCorKnow(dm, insitu.matrix, insitu.genes) 
    geneCor40Know[[i]] = coff$corM
  }else{
    tid = test.data[[i]]
    binDist40[[i]] = scorePost(selGene,tid, bins$top10, ground.truth[tid,], ambig.locations, rawdm)
  } 
}
 

binDist20 = list()
geneCor20 = list()
geneCor20Know = list()
gnum = 20
estDist20 <- list()
geneAll <- matrix(0, 10, gnum)
for(i in 1:10){
  #####1.Training: Use the training data to select genes 
  insitu.genes <- colnames(insitu.matrix)
  scRNA <- t(normalized.data[insitu.genes, train.data[[i]] ])
  selGene <- GetGene_Hc(scRNA, "ward", gnum)# Hclust-mcquitty
  geneAll[i,] <- selGene
  
  dm = new("DistMap",
           raw.data = raw.data,
           data=normalized.data[,train.data[[i]]],
           insitu.matrix=insitu.matrix[,selGene],
           geometry=as.matrix(geometry)) 
  qth <- getBestBinTh(dm, quantiles=seq(0.15, 0.5, 0.01)) 
  
  #####2.Prediction: use the selected genes and prediction data to predict the cell position  
  
  #####2.1  
  dm = new("DistMap",
           raw.data = raw.data,
           data=normalized.data[,test.data[[i]]],
           insitu.matrix=insitu.matrix[,selGene],
           geometry=as.matrix(geometry))
  #####2.2 Binarization
  dm <- binFun(dm, quantiles=qth) 
  #####2.3 Calculate MCC: 3039bins * test cells
  dm <- calCor(dm) 
  
  maxc20 <- apply(dm@mcc.scores, 2, max) 
  v = summary(maxc20)
  #####2.4 Prediction top 10 bins for the test data
  bins <- getTop10(dm@mcc.scores, geometry, v[1])  
  
  if(1==0){
    geneM  <- matrix(selGene, nrow=2,ncol=10)
    res  <- rbind(geneM, bins$top10)
    rownames(res) <- c(rep("",2),  test.data[[i]] )
    fname = paste("post_challenge_zho_team/20genes", i, ".csv", sep = "_")
    write.table(res, file = fname,quote = FALSE, sep=",", 
                col.names=FALSE, row.names = TRUE) 
  }else if(1==2){
    tid = test.data[[i]]
    binDist20[[i]] = esBinDistAll(bins$top10, MCC[,tid], geometry)
    
    coff = getCor(dm, insitu.matrix, insitu.genes)  
    geneCor20[[i]] = coff$corM
    
    coff = getCorKnow(dm, insitu.matrix, insitu.genes) 
    geneCor20Know[[i]] = coff$corM
  }else{
    tid = test.data[[i]]
    binDist20[[i]] = scorePost(selGene,tid, bins$top10, ground.truth[tid,], ambig.locations, rawdm)
  } 
} 