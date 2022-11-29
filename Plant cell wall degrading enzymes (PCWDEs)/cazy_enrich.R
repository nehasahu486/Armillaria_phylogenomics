##############################################################################################################################################
#############################                                                                                    #############################
#############################     Shared Enriched CAZymes in Physalacriaceae and Ascomycota wrt WR/LD species    #############################
#############################                                                                                    #############################
##############################################################################################################################################

rm(list=ls())
list.files()
library(reshape2)
library(ggplot2)
library(phytools)
library(ape)
library(dplyr)
library(ggtree)
library(gridExtra)
library(ggtreeExtra)
library(tibble)
library(stringr)
library(patchwork)
library(ggnewscale)


###STEP1: Table for protein counts in each orthogroup fro different lifestyles
##Enriched CAZymes in Physalacriaceae with respect to LD+WR
pcwde1<-read.delim("pcwde_main_asco.tsv")
head(pcwde1)

##for our enrichment analysis, we need only 3 lifestyles - WR, LD and Physalacriaceae (Physac) species
pcwde2<-pcwde1[c(1,3,6,9)] #select which lifestyles you want
head(pcwde2) #check if they are there

pcwde<-pcwde2[rowSums(pcwde2[,-1])>0,]
head(pcwde)
length(pcwde2$OG_ID)
length(pcwde$OG_ID)

pcwde$A<-pcwde$Physac
pcwde$B<-pcwde$LD+pcwde$WR
pcwde$C<-sum(pcwde$Physac)-pcwde$A
pcwde$D<-sum(pcwde$LD,pcwde$WR)-pcwde$B
#pcwde$D<-sum(pcwde$Physac,pcwde$LD,pcwde$WR)-pcwde$A-pcwde$B-pcwde$C
pcwde$E<-sum(pcwde$A,pcwde$B,pcwde$C,pcwde$D)
head(pcwde)

EM<-pcwde[c(1,5:8)]
head(EM)
tail(EM)
Fstat<-sapply(1:nrow(EM), function(i2) fisher.test(matrix(unlist(EM[i2,2:5]),nrow = 2)))
EM$Pvalue<-unlist(t(Fstat)[,1])
EM$odds<-unlist(t(Fstat)[,3])
EM$FDR <- p.adjust(EM$Pvalue, method= "BH")
PT<-as_tibble(EM$Pvalue)
FDR<-as_tibble(EM$FDR)
FDRL<-abs(log10(FDR[,1]))
FDRL <- data.frame(lapply(FDRL, function(x) {gsub("Inf", "200", x)}))
szig<-FDR
szig[szig[1]<0.0001,]<-4
szig[szig[1]<0.001,]<-3
szig[szig[1]<0.05,]<-2
szig[szig[1]<=1,]<-0
szigm<-szig
szigm[szigm[1]>1,]<-2
szigm[szigm[1]<1,]<-1
EM$sig<-szig$value
EM$sigm<-szigm$value
EM$Dir<-EM$A/EM$B > as.numeric(EM$C)/as.numeric(EM$D)
write.table(EM, file='PhysacvsLDWR.tsv', quote=FALSE, sep='\t', row.names = F)
