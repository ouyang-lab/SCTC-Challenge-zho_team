library(DistMap) 
library(NbClust)
source("challenge_utils.R")
##################################################################################
##1. Load raw data
#################################################################################
 
raw.data = read.table("dge_raw.txt", 
                      sep = "\t",
                      row.names = NULL,
                      stringsAsFactors = F,
                      quote = "")
raw.data.genes = raw.data$V1
raw.data$V1 = NULL
# Let’s fix the gene names that contains apostrophe – this would generate issues
print(grep("'",raw.data.genes,value = T,fixed = T))
raw.data.genes = gsub("'","",raw.data.genes,fixed = T)

raw.data = as.matrix(raw.data)
rownames(raw.data) = raw.data.genes

# Repeat for the normalised data
normalized.data = read.table("dge_normalized.txt",
                             sep = "\t",
                             row.names = NULL,
                             stringsAsFactors = F,
                             quote = "")
normalized.data.genes = normalized.data$row.names
normalized.data$row.names = NULL

# gene names with apostrophes
print(grep("'",normalized.data.genes,value = T,fixed = T))

normalized.data.genes = gsub("'","",normalized.data.genes,fixed = T)

normalized.data = as.matrix(normalized.data)
rownames(normalized.data) = normalized.data.genes
#Check that the gene names are identical in the raw and normalised dataset
stopifnot(all(normalized.data.genes == raw.data.genes))

#Import in situ datasets
insitu.matrix = read.table("binarized_bdtnp.csv", sep = ",",header = T)
insitu.genes_orig <- colnames(insitu.matrix)

# 2 gene names are not matched:
missingGenes = insitu.genes_orig[which(!insitu.genes_orig %in% normalized.data.genes)]
print(missingGenes)
# this was reported by Nikos
# lets fix this by changing the . characters in the gene names to -
insitu.genes = gsub(".","-",insitu.genes_orig,fixed = T)
# also replace .spl. --> (spl)
insitu.genes = gsub("-spl-","(spl)",insitu.genes,fixed = T)

id1 <- which(rownames(normalized.data) =="Blimp-1" )
id2 <- which(rownames(normalized.data) =="E(spl)m5-HLH" )
rownames(normalized.data)[id1] <- "Blimp-1" 
rownames(normalized.data)[id2] <- "E(spl)m5-HLH" 

table(insitu.genes %in% rownames(normalized.data))

# Now we can rename the genes in the institu.matrix with the correct names:
insitu.matrix = as.matrix(insitu.matrix)
colnames(insitu.matrix) = insitu.genes

# Read geometry data
geometry = read.csv("geometry.txt",sep = " ")
colnames(geometry) = c("x","y","z")
dim(geometry)

q2 <-  data.frame(x = geometry[,1], y = geometry[,2], z= geometry[,3]) 
plot_ly(q2, x=~x, y=~y,  z=~z,  mode="markers", # color = as.character(gid),
        marker = list(size = 6,opacity=1), scene='scene1')   


##################################################################################
##2 . Select 60 genes from 84 genes from scRNA data
#################################################################################

scRNA <- t(normalized.data[insitu.genes, ])
gene60 <- GetGene_Hc(scRNA, "mcquitty", 60)# Hclust-mcquitty


### remove duplicate names ‘betaCOP’
id <- which(rownames(normalized.data) =="betaCOP")
plot(normalized.data[id[1],])
plot(normalized.data[id[2],])
normalized.data <- normalized.data[-id[1],]

##################################################################################
##3 . Compute the probablity matrix of a cell belong to a bin position
#################################################################################
maxBins <- apply(MCC, 2, which.max)
maxBins <- geometry[maxBins,] 
 
selGene <- gene60
dm60 = new("DistMap",
         raw.data = raw.data,
         data= normalized.data,
         insitu.matrix=insitu.matrix[,selGene],
         geometry=as.matrix(geometry))
dm60 <- binarizeSingleCellData(dm60, quantiles=seq(0.15, 0.5, 0.01))

dm60 <- calCor(dm60) 
print(dim(dm60@insitu.matrix))
print(dim(dm60@binarized.data))  

##################################################################################
## 4. Predict the cell position
#################################################################################

maxc <- apply(dm60@mcc.scores, 2, max)
plot(maxc)
summary(maxc) 

bin60 <- getTop10(dm60@mcc.scores, geometry, summary(maxc)[5]) 

geneM60 <- matrix(selGene, nrow=6,ncol=10)
res60 <- rbind(geneM60,bin60$top10)
rownames(res60) <- c(rep("",6), 1:nrow(bin60$top10))
#write.table(res60, file = "60genes.csv",quote = FALSE, sep=",", 
#            col.names=FALSE, row.names = TRUE)



