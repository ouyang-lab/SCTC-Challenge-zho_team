################################################### Score 1, 2, 3
## Modified from https://github.com/dream-sctc/Scoring/blob/master/dream_scoring_clean.R
########################## 
score <- function(genes,locations, ground.truth, ambig.locations, rawdm){
  
  if (!exists("rawdm")){
    
    submission <- read.csv(path,header=FALSE,stringsAsFactors = FALSE,na.strings = "")
    
    #separate the gene names from the location predictions
    gene.lines <- (4-sub)*2
    genes <- submission %>% slice(1:gene.lines)
    locations <- submission %>% slice(-1:-gene.lines) 
    
    #preprocess genes and locations, remove NAs, sort locations by cellid
    genes <- na.omit(genes %>% select(-1) %>% unlist %>% as.character)
    locations <- locations[order(locations[,1]),] %>% select(-1) %>% apply(2,as.numeric)
    
  } 
  
  #map every cell to its d84 value
  d84 <<- seq(nrow(ground.truth)) %>% map_dbl(function(j){
    #map every position to the norm of the difference in the geometry and calculate the mean
    ground.truth[j,] %>% map_dbl(~sqrt(sum((rawdm@geometry[.x,] - rawdm@geometry[ground.truth[j,1],])^2))) %>% mean
  })
  
  #do the same mapping for the submission as for d48
  dsub <- seq(nrow(locations)) %>% map_dbl(function(j){
    vals <- locations[j,] %>% map_dbl(~sqrt(sum((rawdm@geometry[.x,] - 
                                                   rawdm@geometry[ground.truth[j,1],])^2))) %>% mean
  })
  
  #calculate relative precision
  pk <- d84/dsub
  summary(pk)
  
  #s1
  
  #select fluorescence data only for the submitted subset of genes
  #reduced.insitu <- data.frame(insitu.matrix) %>% select(genes)
  reduced.insitu <- data.frame(rawdm@insitu.matrix[,genes])
  colnames(reduced.insitu) = genes
  #get binarized data from distmap
  reduced.ts <- data.frame(t(rawdm@binarized.data[genes,]))  
  colnames(reduced.ts) = genes
  
  # map every cell location prediction to the MCC between the ground truth location and 
  # the predicted most likely position, using the submitted subset of genes 
  mccrs <- seq(nrow(locations)) %>% map_dbl(~mccr(reduced.insitu[ground.truth[.x,1],],
                                                  reduced.insitu[locations[.x,1],]))
  
  #do not take into account the cells with ambiguous locations
  s1 <- sum(((pk/sum(pk)) * mccrs)[-ambig.locations])
  
  #s2
  #do not take into account the cells with ambiguous locations
  s2<- mean(pk[-ambig.locations])
  
  #s3
  
  #comparing rnaseq and fluorescence data using true locations
  true.mccs <- seq(ncol(reduced.ts)) %>% map_dbl(~mccr(reduced.insitu[ground.truth[-ambig.locations,1],.x],
                                                       reduced.ts[-ambig.locations,.x]))
  #.. using submitted locations
  competitor.mccs <- seq(ncol(reduced.ts)) %>% map_dbl(~mccr(reduced.insitu[locations[-ambig.locations,1],.x],
                                                             reduced.ts[-ambig.locations,.x]))
  
  #do not take into account the cells with ambiguous locations
  s3 <- sum(((true.mccs/sum(true.mccs)) * competitor.mccs))
  
  #select fluorescence data only for the submitted subset of genes
  #reduced.insitu <- data.frame(insitu.matrix) %>% select(genes)
  reduced.insitu <- data.frame(rawdm@insitu.matrix)
  colnames(reduced.insitu) = colnames(rawdm@insitu.matrix)
  #get binarized data from distmap
  reduced.ts <- data.frame(t(rawdm@binarized.data))  
  colnames(reduced.ts) = rownames(rawdm@binarized.data)
  
  true.mcc84 = seq(ncol(reduced.ts)) %>% map_dbl(~mccr(reduced.insitu[ground.truth[-ambig.locations,1],.x],
                                                       reduced.ts[-ambig.locations,.x]))
  
  competitor.mcc84 = seq(ncol(reduced.ts)) %>% map_dbl(~mccr(reduced.insitu[locations[-ambig.locations,1],.x],
                                                             reduced.ts[-ambig.locations,.x]))
  
  s4 = sum(((true.mcc84/sum(true.mcc84)) * competitor.mcc84))
  return(list(d84=d84, dsub=dsub, score= c(s1,s2,s3,s4), mccrs= mccrs, 
              true.mccs=true.mccs, competitor.mccs=competitor.mccs, 
              true.mcc84=true.mcc84,competitor.mcc84=competitor.mcc84 ))
}

##########For score 1,2,3 on CV results
scorePost <- function(genes,cellid,locations, ground.truth, ambig.locations, rawdm){
  
  
  #map every cell to its d84 value
  d84 <<- seq(nrow(ground.truth)) %>% map_dbl(function(j){
    #map every position to the norm of the difference in the geometry and calculate the mean
    ground.truth[j,] %>% map_dbl(~sqrt(sum((rawdm@geometry[.x,] - rawdm@geometry[ground.truth[j,1],])^2))) %>% mean
  })
  
  #do the same mapping for the submission as for d48
  dsub <- seq(nrow(locations)) %>% map_dbl(function(j){
    vals <- locations[j,] %>% map_dbl(~sqrt(sum((rawdm@geometry[.x,] - 
                                                   rawdm@geometry[ground.truth[j,1],])^2))) %>% mean
  })
  
  #calculate relative precision
  pk <- d84/dsub
  summary(pk)
  
  #s1
  
  #select fluorescence data only for the submitted subset of genes
  #reduced.insitu <- data.frame(insitu.matrix) %>% select(genes)
  reduced.insitu <- data.frame(rawdm@insitu.matrix[,genes])
  colnames(reduced.insitu) = genes
  #get binarized data from distmap
  reduced.ts <- data.frame(t(rawdm@binarized.data[genes,]))  
  colnames(reduced.ts) = genes
  
  # map every cell location prediction to the MCC between the ground truth location and 
  # the predicted most likely position, using the submitted subset of genes 
  mccrs <- seq(nrow(locations)) %>% map_dbl(~mccr(reduced.insitu[ground.truth[.x,1],],
                                                  reduced.insitu[locations[.x,1],]))
  
  ambig.indexes <- which(cellid %in% ambig.locations)
  
  #do not take into account the cells with ambiguous locations
  s1 <- sum(((pk/sum(pk)) * mccrs)[-ambig.indexes])
  
  #s2
  #do not take into account the cells with ambiguous locations
  s2<- mean(pk[-ambig.indexes])
  
  #s3
  
  #comparing rnaseq and fluorescence data using true locations
  true.mccs <- seq(ncol(reduced.ts)) %>% map_dbl(~mccr(reduced.insitu[ground.truth[-ambig.indexes,1],.x],
                                                       reduced.ts[cellid[-ambig.indexes],.x]))
  #.. using submitted locations
  competitor.mccs <- seq(ncol(reduced.ts)) %>% map_dbl(~mccr(reduced.insitu[locations[-ambig.indexes,1],.x],
                                                             reduced.ts[cellid[-ambig.indexes],.x]))
  
  #do not take into account the cells with ambiguous locations
  s3 <- sum(((true.mccs/sum(true.mccs)) * competitor.mccs))
  
  return(list(d84=d84, dsub=dsub, score= c(s1,s2,s3), mccrs= mccrs, 
              true.mccs=true.mccs, competitor.mccs=competitor.mccs))
}

####Figure 2
plotBinMaxvalue <- function(d1,d2,d3){
  comCorData <- data.frame(cbind(1:length(d1), d1, d2,d3)) 
  colnames(comCorData) = c("Cells","60","40","20")
  corMelt <- reshape2::melt(comCorData, id.var='Cells')
  colnames(corMelt) <- c("Cells", "Gene_number", "Probability")
  # 基函数
  ggplot(corMelt, aes(x = factor(Gene_number), y = Probability, fill = factor(Gene_number))) + 
    # 箱线图函数
    geom_boxplot(notch = TRUE) + 
    theme(legend.text=element_text(size=15), 
          axis.text=element_text(size=12), 
          axis.title=element_text(size=14,face="bold")) +  
    ggtitle("Max probability for each cell under selected genes")
}

### Fig.3
plotBinBoxs1 <- function(d1,d2,d3){
  comCorData <- data.frame(cbind(1:length(d1), d1, d2,d3)) 
  colnames(comCorData) = c("Bins","60","40","20")
  corMelt <- reshape2::melt(comCorData, id.var='Bins')
  colnames(corMelt) <- c("Bins", "Gene_number", "Frequency")
  # 基函数
  ggplot(corMelt, aes(x = factor(Gene_number), y = Frequency, fill = factor(Gene_number))) + 
    # 箱线图函数
    geom_boxplot(notch = TRUE) + 
    theme(legend.text=element_text(size=15), 
          axis.text=element_text(size=12), 
          axis.title=element_text(size=14,face="bold")) +  
    ggtitle("Frequency of gene expression for each bin under selected genes")
}


####Fig.5
plotSubplots2 <- function(S1, S2, S3 ){
  library(ggplot2)
  library(gridExtra)
  library(grid) 
  timevec = 1:length(S1)  
  TL = length(timevec)
  d1 = data.frame(cell=timevec, precision=S1)
  d2 = data.frame(cell=timevec, precision=S2)
  d3 = data.frame(cell=timevec, precision=S3) 
  cols = rainbow(8)[seq(1,8,2)]
  
  p1 = ggplot(d1, aes(x=cell, y=precision)) +  geom_point(size=2, shape=0,colour = cols[2]) + 
       geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.5) + ggtitle("Subchallenge 1: 60 genes")
  
  p2 = ggplot(d2, aes(x=cell, y=precision)) +  geom_point(size=2, shape=1,colour = cols[3])+ 
    geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.5)+ ggtitle("Subchallenge 2: 40 genes")
  
  p3 = ggplot(d3, aes(x=cell, y=precision)) +  geom_point(size=2, shape=2,colour = cols[4]) + 
    geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.5)+ ggtitle("Subchallenge 3: 20 genes")
  
  grid.arrange(p1, p2, p3, ncol = 1,
               top="Compare the relative precision for sub challenge 1, 2, 3")
  
}


###Fig.6
plotBarplot <- function(mcc1, mcc2, gene1, size1=20, size2= 6,
                        title="Compare the spatial gene expression in challenge 3"){
  # library
  library(ggplot2) 
   
  comCorData <- data.frame( x=gene1 , y=mcc1, z=mcc2) 
  colnames(comCorData) = c("Gene","84", as.character(length(gene1)))
  corMelt <- reshape2::melt(comCorData, id.var='Gene')
  colnames(corMelt) <- c("Gene", "Group", "MCC") 
  
  # Grouped
  ggplot(corMelt, aes(fill=Group, y=MCC, x=Gene)) + geom_bar(position="dodge", stat="identity")+  
    ggtitle(title)+
    theme(text = element_text(size=size1),
          axis.text.x = element_text(size=size2)) 
}

####Fig.7
predictSocre<- function(binDist, normalized.data,insitu.genes, train.data, raw.data, geometry, rawdm){
  scores = matrix(0, length(binDist), 3)
  for(i in 1:length(binDist)){
    scores[i,] = binDist[[i]]$score
  }
  colMeans(scores)
  sd(scores[,1])
  sd(scores[,2])
  sd(scores[,3])
  i = which.max(scores[,1])
  
  scRNA <- t(normalized.data[insitu.genes, train.data[[i]] ])
  
  if(sel==1){
    selGene <- GetGene_Hc(scRNA, "mcquitty", 60)# Hclust-mcquitty
  }else if(sel==2){
    selGene <- GetGene_Hc(scRNA, "mcquitty", 40)# Hclust-mcquitty
  }else{
    selGene <- GetGene_Hc(scRNA, "ward", 20)# Hclust-mcquitty
  } 
  
  dm = new("DistMap",
           raw.data = raw.data,
           data=normalized.data[,train.data[[i]]],
           insitu.matrix=insitu.matrix[,selGene],
           geometry=as.matrix(geometry)) 
  qth <- getBestBinTh(dm, quantiles=seq(0.15, 0.5, 0.01)) 
  
  
  scRNA <- t(normalized.data[insitu.genes, ])
  dms = new("DistMap",
             raw.data = raw.data,
             data= t(scRNA),
             insitu.matrix=insitu.matrix[,selGene],
             geometry=as.matrix(geometry))
  dms <- binarizeSingleCellData(dms, quantiles=qth) 
  dms <- calCor(dms)  
  bins <- getTop10(dms@mcc.scores, geometry, summary(maxc)[1]) 
  sc = score(selGene,bins$top10, ground.truth, ambig.locations, rawdm)
  sc$score 
}