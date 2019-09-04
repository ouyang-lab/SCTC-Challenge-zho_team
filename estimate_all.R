library(purrr)
library(mccr)
library(tidyverse) 
library(caret)
library(synapser)
library(DistMap)
source("estimate_util.R")
################################################## Raw DistMap
rawdm <- new("DistMap",
             raw.data = raw.data,
             data= normalized.data,
             insitu.matrix=insitu.matrix,
             geometry=as.matrix(geometry))
rawdm <- binarizeSingleCellData(rawdm, quantiles=seq(0.15, 0.5, 0.01))
rawdm <- mapCells(rawdm)
MCC <- rawdm@mcc.scores  
MCC <- log(MCC)
scB84 = rawdm@binarized.data

binarized.data = rawdm@binarized.data
insitu.matrix = rawdm@insitu.matrix

#GROUND TRUTH
ground.truth <<- t(apply(MCC,2,order,decreasing=TRUE))[,1:10]
ambig.locations <<- t(apply(MCC,2,sort,decreasing=TRUE))[,1:2]
ambig.locations <<- which(ambig.locations[,1] == ambig.locations[,2])
length(ambig.locations)

#map every cell to its d84 value
d84 <<- seq(nrow(ground.truth)) %>% map_dbl(function(j){
  #map every position to the norm of the difference in the geometry and calculate the mean
  ground.truth[j,] %>% map_dbl(~sqrt(sum((geometry[.x,] - geometry[ground.truth[j,1],])^2))) %>% mean
})


########################Fig.2
# maxc60, maxc40,  maxc20 are derived from challenge_code_1.R,challenge_code_2.R,challenge_code_3.R 
plotBinMaxvalue(maxc60, maxc40, maxc20)


########################Fig.3
# gene60, gene40,  gene20 are derived from challenge_code_1.R,challenge_code_2.R,challenge_code_3.R 
v1 = rowSums(insitu.matrix[,gene60])
v2 = rowSums(insitu.matrix[,gene40])
v3 = rowSums(insitu.matrix[,gene20])
plotBinBoxs1(v1,v2,v3)

library(VennDiagram)
venn.plot <- venn.diagram(x=list(gene60=gene60, gene40=gene40, gene20=gene20),
                          NULL, height = 250, width = 150, resolution =300, 
                          imagetype="png", col="white", 
                          fill=c(colors()[616], colors()[38], colors()[468]), 
                          alpha=c(0.9, 0.9, 0.9), lwd=c(1, 1, 1), 
                          cex= 1, cat.dist=c(-0.07, -0.07, -0.05), cat.pos=c(300, 60, 180),
                          cat.cex=1, main = "The overlaps among selected genes", main.cex=2)
# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
grid.draw(venn.plot)

########################Fig.4, we used the Excel to draw the bar plot of score 1,2,3
scoreraw = score(insitu.genes, ground.truth, ground.truth, ambig.locations, rawdm) 
score60 = score(gene60,bin60$top10, ground.truth, ambig.locations, rawdm)
score40 = score(gene40,bin40$top10, ground.truth, ambig.locations, rawdm)
score20 = score(gene20,bin20$top10, ground.truth, ambig.locations, rawdm) 

sl60=list()
for(i in 1:6){
  bin60 <- getTop10(dm60@mcc.scores, geometry, summary(maxc)[i])  
  sl60[[i]] = score(gene60,bin60$top10, ground.truth, ambig.locations, rawdm)
}

sl40=list()
for(i in 1:6){
  bin40 <- getTop10(dm40@mcc.scores, geometry, summary(maxc)[i])  
  sl40[[i]] = score(gene40,bin40$top10, ground.truth, ambig.locations, rawdm)
}

sl20=list()
for(i in 1:6){
  bin20 <- getTop10(dm20@mcc.scores, geometry, summary(maxc20)[i]) 
  sl20[[i]] = score(gene20,bin20$top10, ground.truth, ambig.locations, rawdm)
}

######################Fig.5(a). Fig.5(b) is drawn in Excel.
plotSubplots2(p60[-ambig.locations], p40[-ambig.locations], p20[-ambig.locations])

  
######################Fig.6

plotBarplot(score60$true.mccs, score60$competitor.mccs, gene60, size1=20, size2= 8,
            title="Compare predicted spatial gene expression accuracy in challenge 1")
plotBarplot(score40$true.mccs, score40$competitor.mccs, gene40, size1=20, size2= 8,
            title="Compare the predicted spatial gene expression accuracy in challenge 2")
plotBarplot(score20$true.mccs, score20$competitor.mccs, gene20, size1=20, size2= 10,
            title="Compare the predicted spatial gene expression accuracy in challenge 3")


######################Fig.7
#########Fig.7 are drawn by Excel. Users can get the scores for each fold from the binDist60, binDist40
########## and binDist20 list.
##########By use the score 1 from 10 fold CV, select the bese parameters to predicted the cell position for all cells.
binDist = binDist60

predictSocre(binDist60, normalized.data,insitu.genes, train.data, raw.data, geometry, rawdm)
predictSocre(binDist40, normalized.data,insitu.genes, train.data, raw.data, geometry, rawdm)
predictSocre(binDist20, normalized.data,insitu.genes, train.data, raw.data, geometry, rawdm)
  

