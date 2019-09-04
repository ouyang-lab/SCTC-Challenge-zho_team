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

id <- which(rownames(normalized.data) =="betaCOP")
plot(normalized.data[id[1],])
plot(normalized.data[id[2],])
normalized.data <- normalized.data[-id[1],]

save(a=geometry, b= insitu.matrix, c=normalized.data, d=insitu.genes, file = "drosopila.RData")

##################################################################################
## 2. Select 40 genes from 84 genes from scRNA data
#################################################################################

scRNA <- t(normalized.data[insitu.genes, ])
gene40 <- GetGene_Hc(scRNA, "mcquitty", 40)

##################################################################################
##3 . Compute the probablity matrix of a cell belong to a bin position
#################################################################################
selGene <- gene40
dm40 = new("DistMap",
         raw.data = raw.data,
         data= normalized.data,
         insitu.matrix=insitu.matrix[,selGene],
         geometry=as.matrix(geometry))
dm40 <- binarizeSingleCellData(dm40, quantiles=seq(0.15, 0.5, 0.01))

dm40 <- calCor(dm40) # Compute the probablity matrix of a cell belong to a bin position
print(dim(dm40@insitu.matrix))
print(dim(dm40@binarized.data)) 

##################################################################################
## 4. Predict the cell position
#################################################################################
maxc40 <- apply(dm40@mcc.scores, 2, max)
plot(maxc40)
summary(maxc40)

#bin40 <- getTop10(dm40@mcc.scores, geometry, 0.78)
bin40 <- getTop10(dm40@mcc.scores, geometry,summary(maxc40)[2])
 
sl40=list()
for(i in 1:6){
  bin40 <- getTop10(dm40@mcc.scores, geometry,summary(maxc40)[i])
  sl40[[i]] = score(gene40,bin40$top10, ground.truth, ambig.locations, rawdm)
}

geneM40 <- matrix(selGene, nrow=4,ncol=10)
res40 <- rbind(geneM40,bin40$top10)
rownames(res40) <- c(rep("",4), 1:nrow(bin40$top10))
write.table(res40, file = "40genes.csv",quote = FALSE, sep=",", 
            col.names=FALSE, row.names = TRUE)

