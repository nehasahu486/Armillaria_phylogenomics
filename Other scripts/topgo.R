
library(biomaRt)
library(org.Mm.eg.db)
library(topGO)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tidyverse)
library(tibble)
library(GOsummaries, quietly=TRUE)
library(vegan, quietly=TRUE)
library(ggplot2, quietly=TRUE)


geneID2GO <- readMappings(file = "Gene_univ.txt")
View(geneID2GO)
geneUniverse <- names(geneID2GO) 



################UPREG##################
genesOfInterest <- read.csv("8vs3_upreg.csv",header=T)
head(genesOfInterest)
genesOfInterest <- as.character(genesOfInterest$ProtID) 
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse


####Biological Process
myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
resultTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
mysummary <- summary(attributes(resultTopgo)$score <= 0.05)
mysummary
numsignif <- as.integer(mysummary[[3]])
allRes <- GenTable(myGOdata, topgoFisher = resultTopgo, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif)
allRes
write.table(allRes, file="8vs3_upreg_BP.tsv", append="false", sep="\t")
# print a graph (to a pdf file) with the top 'numsignif' results:
showSigOfNodes(myGOdata, score(resultTopgo), firstSigNodes = 5, useInfo ='all')
printGraph(myGOdata, resultTopgo, firstSigNodes = 5, fn.prefix = "8vs3_upreg_BP", useInfo = "all", pdfSW = TRUE)



####Molecular function
myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
resultTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
mysummary <- summary(attributes(resultTopgo)$score <= 0.05)
mysummary
numsignif <- as.integer(mysummary[[3]])
allRes <- GenTable(myGOdata, topgoFisher = resultTopgo, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif)
allRes
write.table(allRes, file="8vs3_upreg_MF.tsv", append="false", sep="\t")
# print a graph (to a pdf file) with the top 'numsignif' results:
showSigOfNodes(myGOdata, score(resultTopgo), firstSigNodes = 5, useInfo ='all')
printGraph(myGOdata, resultTopgo, firstSigNodes = 5, fn.prefix = "8vs3_upreg_MF", useInfo = "all", pdfSW = TRUE)

###Cellular component
myGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
resultTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
mysummary <- summary(attributes(resultTopgo)$score <= 0.05)
mysummary
numsignif <- as.integer(mysummary[[3]])
allRes <- GenTable(myGOdata, topgoFisher = resultTopgo, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif)
allRes
write.table(allRes, file="8vs3_upreg_CC.tsv", append="false", sep="\t")
# print a graph (to a pdf file) with the top 'numsignif' results:
showSigOfNodes(myGOdata, score(resultTopgo), firstSigNodes = 5, useInfo ='all')
printGraph(myGOdata, resultTopgo, firstSigNodes = 5, fn.prefix = "8vs3_upreg_CC", useInfo = "all", pdfSW = TRUE)

#########################DOWNREG######################
genesOfInterest <- read.csv("8vs3_downreg.csv",header=T)
head(genesOfInterest)
genesOfInterest <- as.character(genesOfInterest$ProtID) 
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse


####Biological Process
myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
resultTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
mysummary <- summary(attributes(resultTopgo)$score <= 0.05)
mysummary
numsignif <- as.integer(mysummary[[3]])
allRes <- GenTable(myGOdata, topgoFisher = resultTopgo, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif)
allRes
write.table(allRes, file="8vs3_downreg_BP.tsv", append="false", sep="\t")
# print a graph (to a pdf file) with the top 'numsignif' results:
showSigOfNodes(myGOdata, score(resultTopgo), firstSigNodes = 5, useInfo ='all')
printGraph(myGOdata, resultTopgo, firstSigNodes = 5, fn.prefix = "8vs3_upreg_BP", useInfo = "all", pdfSW = TRUE)



####Molecular function
myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
resultTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
mysummary <- summary(attributes(resultTopgo)$score <= 0.05)
mysummary
numsignif <- as.integer(mysummary[[3]])
allRes <- GenTable(myGOdata, topgoFisher = resultTopgo, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif)
allRes
write.table(allRes, file="8vs3_downreg_MF.tsv", append="false", sep="\t")
# print a graph (to a pdf file) with the top 'numsignif' results:
showSigOfNodes(myGOdata, score(resultTopgo), firstSigNodes = 5, useInfo ='all')
printGraph(myGOdata, resultTopgo, firstSigNodes = 5, fn.prefix = "8vs3_upreg_MF", useInfo = "all", pdfSW = TRUE)

###Cellular component
myGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
resultTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
mysummary <- summary(attributes(resultTopgo)$score <= 0.05)
mysummary
numsignif <- as.integer(mysummary[[3]])
allRes <- GenTable(myGOdata, topgoFisher = resultTopgo, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif)
allRes
write.table(allRes, file="8vs3_downreg_CC.tsv", append="false", sep="\t")
# print a graph (to a pdf file) with the top 'numsignif' results:
showSigOfNodes(myGOdata, score(resultTopgo), firstSigNodes = 5, useInfo ='all')
printGraph(myGOdata, resultTopgo, firstSigNodes = 5, fn.prefix = "8vs3_upreg_CC", useInfo = "all", pdfSW = TRUE)


#####The following line gives all the p-values and then can be used further for p-value correction#########
#####Instead of specifying topNodes to some value N, you can get all available GO Terms using##########
allGO = usedGO(object = GOdata) 
# use it in GenTable as follows:
GenTable(GOdata, ... ,topNodes = length(allGO))

