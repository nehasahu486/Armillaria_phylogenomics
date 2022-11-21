rm(list=ls())

setwd("C:/Users/Dell/Documents/TopGO/Armi_compgen_updated/CAZy_PCWDEs/CAZy_enrichments/")
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


###MAKE the OG_ID_count file
og<-read.delim("ProtIDs with OG_IDs and IPRs.tsv")
og$Species<-sub('_[^_]*$','',og$ProtID)
head(og)

pcwde<-read.delim("PCWDE_clusters.tsv")
head(pcwde)
min(pcwde$cazy_perc)
length(unique(pcwde$OGID))

og1<-og[og$OrthoID %in% pcwde$OGID,]
head(og1)
length(unique(og1$OrthoID))

write.table(table(og1$Species,og1$OrthoID),"table.tsv",quote=F,row.names = T,sep="\t")

###main file is pcwde_main_asco

###################################################################################################################
#############Enrichments#####################

####taking only wood decay species in consideration _only LD_WR i.e. BaSR vs LD+WR
pcwde1<-read.delim("pcwde_main_asco.tsv")
head(pcwde1)

pcwde2<-pcwde1[c(1,3,6,9)] #select which lifestyles you want
head(pcwde2) #check if they are there

pcwde<-pcwde2[rowSums(pcwde2[,-1])>0,]
head(pcwde)
length(pcwde2$OG_ID)
length(pcwde$OG_ID)

pcwde$A<-pcwde$Ba_SR
pcwde$B<-pcwde$LD+pcwde$WR
pcwde$C<-sum(pcwde$Ba_SR)-pcwde$A
pcwde$D<-sum(pcwde$LD,pcwde$WR)-pcwde$B
#pcwde$D<-sum(pcwde$Ba_SR,pcwde$LD,pcwde$WR)-pcwde$A-pcwde$B-pcwde$C
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
write.table(EM, file='BaSRvsLDWR.tsv', quote=FALSE, sep='\t', row.names = F)


####Asco vs LD+WR
pcwde1<-read.delim("pcwde_main_asco.tsv")
head(pcwde1)

pcwde2<-pcwde1[c(1,2,6,9)] #select which lifestyles you want
head(pcwde2) #check if they are there

pcwde<-pcwde2[rowSums(pcwde2[,-1])>0,]
head(pcwde)
length(pcwde2$OG_ID)
length(pcwde$OG_ID)

pcwde$A<-pcwde$Asco
pcwde$B<-pcwde$LD+pcwde$WR
pcwde$C<-sum(pcwde$Asco)-pcwde$A
pcwde$D<-sum(pcwde$LD,pcwde$WR)-pcwde$B
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
write.table(EM, file='AscovsLDWR.tsv', quote=FALSE, sep='\t', row.names = F)



###################################################################################################################
#############Enrichments - cazy classes#####################
main<-read.delim("cazy_class_main.tsv")
head(main)


####taking only wood decay species in consideration _only LD_WR i.e. BaSR vs LD+WR
cazy1<-main
head(cazy1)

cazy2<-cazy1[c(1,2,5,9)] #select which lifestyles you want
head(cazy2) #check if they are there

cazy<-cazy2[rowSums(cazy2[,-1])>0,]
head(cazy)
length(rownames(cazy2))
length(rownames(cazy))

cazy$A<-cazy$Ba_SR
cazy$B<-cazy$LD+cazy$WR
cazy$C<-sum(cazy$Ba_SR)-cazy$A
cazy$D<-sum(cazy$Ba_SR,cazy$LD,cazy$WR)-cazy$A-cazy$B-cazy$C
cazy$E<-sum(cazy$A,cazy$B,cazy$C,cazy$D)
head(cazy)

EM<-cazy[c(1,5:8)]
head(EM)
tail(EM)
Fstat<-sapply(1:nrow(EM), function(i2) fisher.test(matrix(unlist(EM[i2,2:5]),nrow = 2)))
EM$Pvalue<-unlist(t(Fstat)[,1])
EM$FDR <- p.adjust(EM$Pvalue, method= "BH")
PT<-as.tibble(EM$Pvalue)
FDR<-as.tibble(EM$FDR)
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
write.table(EM, file='BaSRvsLDandWR_CAZYenriched.tsv', quote=FALSE, sep='\t', row.names = F)


####SR vs LD+WR
cazy1<-main
head(cazy1)

cazy2<-cazy1[c(1,5,7,9)] #select which lifestyles you want
head(cazy2) #check if they are there

cazy<-cazy2[rowSums(cazy2[,-1])>0,]
head(cazy)
length(rownames(cazy2))
length(rownames(cazy))

cazy$A<-cazy$SR
cazy$B<-cazy$LD+cazy$WR
cazy$C<-sum(cazy$SR)-cazy$A
cazy$D<-sum(cazy$SR,cazy$LD,cazy$WR)-cazy$A-cazy$B-cazy$C
cazy$E<-sum(cazy$A,cazy$B,cazy$C,cazy$D)
head(cazy)

EM<-cazy[c(1,5:8)]
head(EM)
tail(EM)
Fstat<-sapply(1:nrow(EM), function(i2) fisher.test(matrix(unlist(EM[i2,2:5]),nrow = 2)))
EM$Pvalue<-unlist(t(Fstat)[,1])
EM$FDR <- p.adjust(EM$Pvalue, method= "BH")
PT<-as.tibble(EM$Pvalue)
FDR<-as.tibble(EM$FDR)
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
write.table(EM, file='SRvsLDandWR_CAZYenriched.tsv', quote=FALSE, sep='\t', row.names = F)


##### enriched in physac wrt to other WR
head(main)
main2<-main[c(1,2,9)] #select which lifestyles you want
head(main2) #check if they are there

main<-main2[rowSums(main2[,-1])>0,]
head(main)
length(main2$OG_ID)
length(main$OG_ID)

main$A<-main$Ba_SR
main$B<-main$WR
main$C<-sum(main$Ba_SR)-main$A
main$D<-sum(main$Ba_SR,main$WR)-main$A-main$B-main$C
main$E<-sum(main$A,main$B,main$C,main$D)
head(main)

EM<-main[c(1,4:7)]
head(EM)
tail(EM)
Fstat<-sapply(1:nrow(EM), function(i2) fisher.test(matrix(unlist(EM[i2,2:5]),nrow = 2)))
EM$Pvalue<-unlist(t(Fstat)[,1])
EM$FDR <- p.adjust(EM$Pvalue, method= "BH")
PT<-as.tibble(EM$Pvalue)
FDR<-as.tibble(EM$FDR)
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
write.table(EM, file='PhysacvsWR_enriched.tsv', quote=FALSE, sep='\t', row.names = F)




pcwde<-read.delim("CAZy_clust.tsv")
head(pcwde)
rownames(pcwde)
pcwde$A<-pcwde$Ba_SR
pcwde$B<-pcwde$LD+pcwde$SR+pcwde$WR+pcwde$BR+pcwde$ECM+pcwde$PATH+pcwde$UR
pcwde$C<-sum(pcwde$Ba_SR)-pcwde$A
pcwde$D<-sum(pcwde$Ba_SR,pcwde$LD,pcwde$SR,pcwde$WR,pcwde$BR,pcwde$ECM,pcwde$PATH,pcwde$UR)-pcwde$A-pcwde$B-pcwde$C
pcwde$E<-sum(pcwde$A,pcwde$B,pcwde$C,pcwde$D)
head(pcwde)

EM<-pcwde[c(1,10,11,12,13)]
head(EM)
tail(EM)
Fstat<-sapply(1:nrow(EM), function(i2) fisher.test(matrix(unlist(EM[i2,2:5]),nrow = 2)))
EM$Pvalue<-unlist(t(Fstat)[,1])
EM$FDR <- p.adjust(EM$Pvalue, method= "BH")
library(tibble)
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
write.table(EM, file='BaSRvsallOthers_enriched.tsv', quote=FALSE, sep='\t', row.names = F)


####taking only wood decay species in consideration
pcwde<-read.delim("CAZy_clust.tsv")
head(pcwde)
pcwde$A<-pcwde$Ba_SR
pcwde$B<-pcwde$LD+pcwde$SR+pcwde$WR+pcwde$BR
pcwde$C<-sum(pcwde$Ba_SR)-pcwde$A
pcwde$D<-sum(pcwde$Ba_SR,pcwde$LD,pcwde$SR,pcwde$WR,pcwde$BR)-pcwde$A-pcwde$B-pcwde$C
pcwde$E<-sum(pcwde$A,pcwde$B,pcwde$C,pcwde$D)
head(pcwde)

EM<-pcwde[c(1,10,11,12,13)]
head(EM)
tail(EM)
Fstat<-sapply(1:nrow(EM), function(i2) fisher.test(matrix(unlist(EM[i2,2:5]),nrow = 2)))
EM$Pvalue<-unlist(t(Fstat)[,1])
EM$FDR <- p.adjust(EM$Pvalue, method= "BH")
PT<-as.tibble(EM$Pvalue)
FDR<-as.tibble(EM$FDR)
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
write.table(EM, file='BaSRvsOtheWD_enriched.tsv', quote=FALSE, sep='\t', row.names = F)


####taking only wood decay species in consideration _only LD_WR i.e. BaSR vs LD+WR
pcwde<-read.delim("CAZy_clust.tsv")
head(pcwde)
pcwde$A<-pcwde$Ba_SR
pcwde$B<-pcwde$LD+pcwde$WR
pcwde$C<-sum(pcwde$Ba_SR)-pcwde$A
pcwde$D<-sum(pcwde$Ba_SR,pcwde$LD,pcwde$WR)-pcwde$A-pcwde$B-pcwde$C
pcwde$E<-sum(pcwde$A,pcwde$B,pcwde$C,pcwde$D)
head(pcwde)

EM<-pcwde[c(1,10,11,12,13)]
head(EM)
tail(EM)
Fstat<-sapply(1:nrow(EM), function(i2) fisher.test(matrix(unlist(EM[i2,2:5]),nrow = 2)))
EM$Pvalue<-unlist(t(Fstat)[,1])
EM$FDR <- p.adjust(EM$Pvalue, method= "BH")
PT<-as.tibble(EM$Pvalue)
FDR<-as.tibble(EM$FDR)
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
write.table(EM, file='BaSRvsLDandWR_enriched.tsv', quote=FALSE, sep='\t', row.names = F)


####BaSR+SR vs LD+WR
pcwde<-read.delim("CAZy_clust.tsv")
head(pcwde)
pcwde$A<-pcwde$Ba_SR+pcwde$SR
pcwde$B<-pcwde$LD+pcwde$WR
pcwde$C<-sum(pcwde$Ba_SR)+sum(pcwde$SR)-pcwde$A
pcwde$D<-sum(pcwde$Ba_SR,pcwde$SR,pcwde$LD,pcwde$WR)-pcwde$A-pcwde$B-pcwde$C
pcwde$E<-sum(pcwde$A,pcwde$B,pcwde$C,pcwde$D)
head(pcwde)
EM<-pcwde[c(1,10,11,12,13)]
head(EM)
tail(EM)
Fstat<-sapply(1:nrow(EM), function(i2) fisher.test(matrix(unlist(EM[i2,2:5]),nrow = 2)))
EM$Pvalue<-unlist(t(Fstat)[,1])
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
write.table(EM, file='BaSR_SRvsLDandWR_enriched.tsv', quote=FALSE, sep='\t', row.names = F)

####SRvsLD+WR
pcwde<-read.delim("CAZy_clust.tsv")
head(pcwde)
pcwde$A<-pcwde$SR
pcwde$B<-pcwde$LD+pcwde$WR
pcwde$C<-sum(pcwde$SR)-pcwde$A
pcwde$D<-sum(pcwde$SR,pcwde$LD,pcwde$WR)-pcwde$A-pcwde$B-pcwde$C
pcwde$E<-sum(pcwde$A,pcwde$B,pcwde$C,pcwde$D)
head(pcwde)

EM<-pcwde[c(1,10,11,12,13)]
head(EM)
tail(EM)
Fstat<-sapply(1:nrow(EM), function(i2) fisher.test(matrix(unlist(EM[i2,2:5]),nrow = 2)))
EM$Pvalue<-unlist(t(Fstat)[,1])
EM$FDR <- p.adjust(EM$Pvalue, method= "BH")
PT<-as.tibble(EM$Pvalue)
FDR<-as.tibble(EM$FDR)
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
write.table(EM, file='SRvsLDandWR_enriched.tsv', quote=FALSE, sep='\t', row.names = F)

####BaSRvsWR
pcwde<-read.delim("CAZy_clust.tsv")
head(pcwde)
pcwde$A<-pcwde$Ba_SR
pcwde$B<-pcwde$WR
pcwde$C<-sum(pcwde$Ba_SR)-pcwde$A
pcwde$D<-sum(pcwde$Ba_SR,pcwde$WR)-pcwde$A-pcwde$B-pcwde$C
pcwde$E<-sum(pcwde$A,pcwde$B,pcwde$C,pcwde$D)
head(pcwde)

EM<-pcwde[c(1,10,11,12,13)]
head(EM)
tail(EM)
Fstat<-sapply(1:nrow(EM), function(i2) fisher.test(matrix(unlist(EM[i2,2:5]),nrow = 2)))
EM$Pvalue<-unlist(t(Fstat)[,1])
EM$FDR <- p.adjust(EM$Pvalue, method= "BH")
PT<-as.tibble(EM$Pvalue)
FDR<-as.tibble(EM$FDR)
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
write.table(EM, file='BaSRvsWR_enriched.tsv', quote=FALSE, sep='\t', row.names = F)


######significantly enriched pcwde plots

#read the tree first - so that we can sort the order of the bubble plot according to the order of tips in the tree
tree<-read.nexus("132sps_final_rooted.nex")
plotTree(tree)
splist<-data.frame(tree$tip.label)
head(splist)
colnames(splist)<-c("Species") #order of species to be used


#check first which OG_IDs to plot
p<-read.delim("overrep_pcwde.tsv")
head(p)
p1<-melt(p)
head(p1)
colnames(p1)<-c("Species","Lifestyle","OG_ID","CAZy")
head(p1)

p3<-left_join(splist,p1,by="Species")
head(p3)
p3$Species<-as.character(p3$Species)
p3$Species<-factor(p3$Species,levels=unique(p3$Species))
head(p3)
p2<-ggplot(p3,aes(x=Species,y=CAZy,fill=Lifestyle))+
  geom_bar(stat="identity",alpha=0.4,position="dodge")+
  theme_classic()+scale_fill_brewer(palette="Dark2")+
  ggtitle("All_sig_ovverep_OGIDs")+coord_flip()+
  facet_grid(~OG_ID, scales = "free")
p2


#read the CAZY_size file - this will be used to create bars
cazysize<-read.delim("CAZy_sizes.tsv")
head(cazysize)
sum(cazysize$CAZy_size)
cazysize$OG_ID<-c("CAZy sizes")
head(cazysize)
head(splist)
cazysize<-left_join(splist,cazysize,by="Species")
cazysize$Species<-as.character(cazysize$Species)
cazysize$Species<-factor(cazysize$Species,levels=unique(cazysize$Species))
head(cazysize)


cazysize_plot<-ggplot(cazysize,aes(x=OG_ID,y=Species,size=CAZy_size,fill=Lifestyle,alpha=0.2))+
  geom_point(aes(alpha=0.2,color=Lifestyle))+
  theme_minimal()+scale_fill_brewer(palette="Dark2",guide="none")+
  scale_color_brewer(palette = "Dark2",guide="none")+ggtitle("")+
  scale_x_discrete(guide=guide_axis(angle=90))+xlab("")+ylab("")+
  theme(panel.grid.major = element_blank(), legend.position = "right",axis.title.y=element_blank(),
        axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x =element_blank())

cazysize_plot

cazysize_plot<-ggplot(cazysize,aes(x=Species,y=CAZy_size,fill=Lifestyle))+
  geom_bar(stat="identity",alpha=0.4,position="dodge")+
  theme_classic()+scale_fill_brewer(palette="Dark2")+
  ggtitle("CAZy sizes")+coord_flip()

cazysize_plot<-ggplot(cazysize,aes(x=Species,y=CAZy_size,fill=Lifestyle))+
  geom_bar(stat="identity",alpha=0.4,position="dodge")+
  theme_minimal()+scale_fill_brewer(palette="Dark2")+
  ggtitle("CAZy sizes")+coord_flip()+
  theme(legend.position = "none",axis.title.y=element_blank(),
      axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x =element_text(size=10) )


cazysize_plot


#read the 20 overlap overrep_pcwde file
p<-read.delim("20overrep_pcwde.tsv")
head(p)
#based on previous OG_ID phylobarplot - selected 27OG_OIDs which showed physalacriaceae specific gene copies
#subset those 25 OG_IDs to make a plot for
colnames(p)
#p<-p[-c(3,4,5,6,7,33,34)]
#head(p)
p1<-melt(p)
head(p1)
colnames(p1)<-c("Species","Lifestyle","OG_ID","CAZy")
head(p1)
#intersect(p1$Species,splist$tree.tip.label)
p3<-left_join(splist,p1,by="Species")
head(p3)
p3$Species<-as.character(p3$Species)
p3$Species<-factor(p3$Species,levels=unique(p3$Species))
head(p3)

p2<-ggplot(p3,aes(x=OG_ID,y=Species,size=ifelse(CAZy==0,NA,CAZy),fill=Lifestyle))+
  geom_point(aes(alpha=0.4,color=Lifestyle))+
  theme_classic()+scale_fill_brewer(palette="Dark2",guide="none")+
  scale_color_brewer(palette = "Dark2",guide="none")+ggtitle("")+
  scale_x_discrete(guide=guide_axis(angle=90))
p2

p2<-ggplot(p3,aes(x=OG_ID,y=Species,size=ifelse(CAZy==0,NA,CAZy),fill=Lifestyle))+
  geom_point(aes(alpha=0.4,color=Lifestyle))+
  theme_classic()+scale_fill_brewer(palette="Dark2",guide="none")+
  scale_color_brewer(palette = "Dark2",guide="none")+ggtitle("")+
  scale_x_discrete(guide=guide_axis(angle=90))+
  theme(legend.position = "none",axis.title.y=element_blank(),
      axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x =element_text(size=10) )
p2

p2<-ggplot(p3,aes(x=Species,y=CAZy,fill=Lifestyle))+
  geom_bar(stat="identity",alpha=0.4,position="dodge")+
  theme_classic()+scale_fill_brewer(palette="Dark2")+
  ggtitle("")+coord_flip()+
  facet_grid(~OG_ID, scales = "free")+xlab("")+ylab("")+
  theme(strip.background=element_blank(),strip.text = element_text(angle=45,size=8),
        legend.position = "bottom",axis.title.y=element_blank(),
        axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x =element_text(size=5) )
p2

p4<-ggtree(tree)
#p4<-ggtree(tree,branch.length = "none",size=0.5)
p4<-p4 %<+% p3 +geom_tiplab(aes(color=Lifestyle),size=3,align=T)+
  scale_color_brewer(palette="Dark2")+theme(legend.position = "none")
p4

p2
p4
p4+p2+plot_layout(widths=c(3,7))
p4+p2+cazysize_plot+plot_layout(widths=c(3,9,0.2))



p<-read.delim("overrep_pcwde.tsv",row.names = 1)
head(p)
#based on previous OG_ID phylobarplot - selected 27OG_OIDs which showed physalacriaceae speific gene copies
#subset those 25 OG_IDs to make a plot for
colnames(p)
p<-p[-c(1,2,3,4,5,6,32,33)]
head(p)
colnames(p)

p<-as.matrix(p)[,2]
head(p)
#View(p)
fit<-phytools::fastAnc(tree,p,vars=T,CI=T)
head(fit)
tf<-data.frame(node=nodeid(tree, names(p)),
               trait=p)
nd<-data.frame(node=names(fit$ace),trait=fit$ace)
d<-rbind(tf,nd)
d$node<-as.numeric(d$node)
tree1<-full_join(tree,d,by="node")

pheno<-ggtree(tree1, aes(color=trait), continuous = "color", yscale = "trait")+ 
  scale_color_viridis_d()+theme_minimal()+geom_tiplab(aes(color=trait))
pheno







list.files()
allanot<-read.delim("All_annotations.tsv", quote="")
head(allanot)
colnames(allanot)
allanot<-unique(allanot[c(2,3,4,6)])
head(allanot)
colnames(allanot)
anot<-aggregate(.~OrthoID, allanot, toString)
head(anot)
tail(anot)
write.table(anot,"OG_ID_anot.tsv",quote=F,row.names = F,sep="\t")



#read the CAZY_family file - this will be used to create bubble plots
tree<-read.nexus("132sps_final_rooted.nex")
plotTree(tree)
splist<-data.frame(tree$tip.label)
head(splist)
colnames(splist)<-c("Species") #order of species to be used

cazyfam<-read.delim("CAZy_bubl_mainfams.tsv")
head(cazyfam)

cazyfam1<-melt(cazyfam)
head(cazyfam1)
colnames(cazyfam1)<-c("Species","Lifestyle","Type","CAZy","Number")
head(cazyfam1)

p3<-left_join(splist,cazyfam1,by="Species")
head(p3)
p3$Species<-as.character(p3$Species)
p3$Species<-factor(p3$Species,levels=unique(p3$Species))
head(p3)


p2<-ggplot(p3,aes(x=CAZy,y=Species,size=ifelse(Number==0,NA,Number),
                  fill=Lifestyle))+
  geom_point(aes(alpha=0.4,color=Lifestyle))+
  theme_classic()+scale_fill_brewer(palette="Dark2")+
  scale_color_brewer(palette = "Dark2")+ggtitle("")+
  scale_x_discrete(guide=guide_axis(angle=90))+
  facet_wrap(~Type,scales="free")
p2

head(p3)
p2<-ggplot(p3,aes(x=Species,y=Number,fill=Lifestyle))+
  geom_bar(stat="identity",alpha=0.4,position="dodge")+
  theme_classic()+scale_fill_brewer(palette="Dark2")+
  ggtitle("")+coord_flip()+
  facet_grid(.~Type +CAZy, scales = "free")+xlab("")+ylab("")
  
p2



###cicrular barplot with 7 OG_IDs sig overrepresented in BaSE and SR (each) vs WR+LD
#read the tree first - so that we can sort the order of the bubble plot according to the order of tips in the tree
tree<-read.nexus("132sps_final_rooted.nex")
plotTree(tree)
splist<-data.frame(tree$tip.label)
head(splist)
colnames(splist)<-c("Species") #order of species to be used


#check first which OG_IDs to plot
p<-read.delim("SR_BaSR_common.tsv")
head(p)
p_sorted<-left_join(splist,p,by="Species")
head(p_sorted)
colnames(p_sorted)
first<-p_sorted[c(1,2,3)]
head(first)
sec<-p_sorted[c(1,2,4)]
thrd<-p_sorted[c(1,2,5)]
forth<-p_sorted[c(1,2,6)]
fifth<-p_sorted[c(1,2,7)]
sixth<-p_sorted[c(1,2,8)]
seventh<-p_sorted[c(1,2,9)]




vert<-ggtree(tree)+geom_tiplab(size=2,align=T)
vert
vert1<-vert+geom_fruit(data=first,geom=geom_bar,mapping=aes(y=Species,x=OG0000701,fill=Lifestyle,alpha=0.1),pwidth=1,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")+
  geom_fruit(data=sec,geom=geom_bar,mapping=aes(y=Species,x=OG0000781,fill=Lifestyle,alpha=0.1),pwidth=1,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")+
  geom_fruit(data=thrd,geom=geom_bar,mapping=aes(y=Species,x=OG0003599,fill=Lifestyle,alpha=0.1),pwidth=1,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")+
  geom_fruit(data=forth,geom=geom_bar,mapping=aes(y=Species,x=OG0006508,fill=Lifestyle,alpha=0.1),pwidth=1,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")+
 geom_fruit(data=fifth,geom=geom_bar,mapping=aes(y=Species,x=OG0006823,fill=Lifestyle,alpha=0.1),pwidth=1,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")+
 geom_fruit(data=sixth,geom=geom_bar,mapping=aes(y=Species,x=OG0007180,fill=Lifestyle,alpha=0.1),pwidth=1,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")+
 geom_fruit(data=seventh,geom=geom_bar,mapping=aes(y=Species,x=OG0007224,fill=Lifestyle,alpha=0.1),pwidth=1,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")

vert1


circ<-ggtree(tree,layout="circular")+geom_tiplab(size=1)
circ<-ggtree(tree,layout="circular") #without tip labels

circ
circ1<-circ+geom_fruit(data=first,geom=geom_bar,mapping=aes(y=Species,x=OG0000701,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")+
  geom_fruit(data=sec,geom=geom_bar,mapping=aes(y=Species,x=OG0000781,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")+
  geom_fruit(data=thrd,geom=geom_bar,mapping=aes(y=Species,x=OG0003599,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")+
  geom_fruit(data=forth,geom=geom_bar,mapping=aes(y=Species,x=OG0006508,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")+
  geom_fruit(data=fifth,geom=geom_bar,mapping=aes(y=Species,x=OG0006823,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")+
  geom_fruit(data=sixth,geom=geom_bar,mapping=aes(y=Species,x=OG0007180,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")+
  geom_fruit(data=seventh,geom=geom_bar,mapping=aes(y=Species,x=OG0007224,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")

circ1

circ1<-circ+
  geom_fruit(data=seventh,geom=geom_bar,mapping=aes(y=Species,x=OG0007224,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")+
  geom_fruit(data=sixth,geom=geom_bar,mapping=aes(y=Species,x=OG0007180,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")+
  geom_fruit(data=fifth,geom=geom_bar,mapping=aes(y=Species,x=OG0006823,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")+
  geom_fruit(data=forth,geom=geom_bar,mapping=aes(y=Species,x=OG0006508,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")+
  geom_fruit(data=thrd,geom=geom_bar,mapping=aes(y=Species,x=OG0003599,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")+
  geom_fruit(data=sec,geom=geom_bar,mapping=aes(y=Species,x=OG0000781,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")+
  geom_fruit(data=first,geom=geom_bar,mapping=aes(y=Species,x=OG0000701,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")
  

circ1


library(forcats)
library(gghighlight)

x<-read.delim("Enrich.tsv", header = T)
head(x)
p<-ggplot(x,aes(x=GeneRatio, y=fct_reorder(Name,GeneRatio)))+
  geom_point(aes(size=GeneRatio, color = FDR))+
  theme_bw(base_size=14)+
  scale_color_gradient2(limits=c(0,0.05), low="green", mid = "orange", high="yellow", midpoint = 0.025)+
  ylab(NULL) + ggtitle("CAZy_enrichment")
p4<-p+facet_grid(.~Type)+theme_minimal()
p4
#
x<-read.delim("Enrichfig.tsv", header = T)
head(x)

p<-ggplot(x,aes(x=Geneperc, y=fct_reorder(Name,Geneperc)))+
  geom_point(aes(size=Geneperc, color = FDR))+
  theme_bw(base_size=10)+
  scale_color_gradient2(limits=c(0,0.05), low="cadetblue3", mid = "yellowgreen", high="orange", midpoint = 0.025)+
  ylab(NULL) + ggtitle("CAZy_enrichment")+theme_minimal()+facet_grid(.~Type)+
  gghighlight(Text>10000,
              unhighlighted_params = list(colour=NULL,alpha=0.3),calculate_per_facet = T)

p

x<-read.delim("Enrichfig_BaSRvsLDWR.tsv", header = T)
head(x)

p<-ggplot(x,aes(x=Geneperc, y=fct_reorder(Name,Geneperc)))+
  geom_point(aes(size=Geneperc, color = FDR))+
  theme_bw(base_size=10)+
  scale_color_gradient2(limits=c(0,0.05), low="cadetblue3", mid = "yellowgreen", high="orange", midpoint = 0.025)+
  ylab(NULL) + ggtitle("CAZy_enrichment")+theme_minimal()+facet_grid(.~Type)+
  gghighlight(Text=="two",
              unhighlighted_params = list(colour=NULL,alpha=0.3),calculate_per_facet = T)

p


###enriched phylobar/phyloheatmap
tree<-read.nexus("132sps_final_rooted.nex")
plotTree(tree)
tree$tip.label
plot(tree)
splist<-data.frame(tree$tip.label)
head(splist)
colnames(splist)<-c("Species")
head(splist)

p<-read.delim("Co-shared_OGs_CN.tsv")
head(p)
p_sorted<-left_join(splist,p,by="Species")
head(p_sorted)
colnames(p_sorted)
first<-p_sorted[c(1,2,3)]
head(first)
sec<-p_sorted[c(1,2,4)]
thrd<-p_sorted[c(1,2,5)]
forth<-p_sorted[c(1,2,6)]
fifth<-p_sorted[c(1,2,7)]
sixth<-p_sorted[c(1,2,8)]
seven<-p_sorted[c(1,2,9)]
eight<-p_sorted[c(1,2,10)]
ninth<-p_sorted[c(1,2,11)]
tenth<-p_sorted[c(1,2,12)]
eleven<-p_sorted[c(1,2,13)]
twel<-p_sorted[c(1,2,14)]
thirt<-p_sorted[c(1,2,15)]
fort<-p_sorted[c(1,2,16)]
fift<-p_sorted[c(1,2,17)]
sixteen<-p_sorted[c(1,2,18)]


data<-ggtree(tree)
data<-data$data

lf<-read.delim("Lifestyle_asco.tsv")
head(lf)
colnames(lf)[1]<-"label"


circ<-ggtree(tree,layout="circular",branch.length = "none")+geom_tiplab(size=1)
circ<-ggtree(tree,layout="circular",color="darkgray") #without tip labels
circ$data
circ
dat<-left_join(lf,circ,by="label")
head(dat)
dat1<-dat[c(1,2,4)]
head(dat1)
tail(dat1)


circ_col<-circ+geom_hilight(data=dat1,mapping=aes(node=node,fill=LF),extendto=3,alpha=0.5)+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "bottom")
circ_col

circ1<-circ_col+new_scale_fill()+
  geom_fruit(data=sixteen,geom=geom_bar,mapping=aes(y=Species,x=OG0009063,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=fift,geom=geom_bar,mapping=aes(y=Species,x=OG0008824,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=fort,geom=geom_bar,mapping=aes(y=Species,x=OG0007224,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=thirt,geom=geom_bar,mapping=aes(y=Species,x=OG0007180,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=twel,geom=geom_bar,mapping=aes(y=Species,x=OG0006856,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=eleven,geom=geom_bar,mapping=aes(y=Species,x=OG0006823,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=tenth,geom=geom_bar,mapping=aes(y=Species,x=OG0006508,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=ninth,geom=geom_bar,mapping=aes(y=Species,x=OG0003599,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=eight,geom=geom_bar,mapping=aes(y=Species,x=OG0002730,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=seven,geom=geom_bar,mapping=aes(y=Species,x=OG0001684,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=sixth,geom=geom_bar,mapping=aes(y=Species,x=OG0001517,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=fifth,geom=geom_bar,mapping=aes(y=Species,x=OG0001290,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=forth,geom=geom_bar,mapping=aes(y=Species,x=OG0001157,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=thrd,geom=geom_bar,mapping=aes(y=Species,x=OG0000893,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=sec,geom=geom_bar,mapping=aes(y=Species,x=OG0000781,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=first,geom=geom_bar,mapping=aes(y=Species,x=OG0000701,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  scale_fill_brewer(palette="Set2")+theme(legend.position = "none")


circ1


vert<-ggtree(tree,color="darkgray")+geom_tiplab(size=4)
vert
vert<-ggtree(tree,color="darkgray") #without tip labels
data<-vert$data

vert

dat<-left_join(lf,data,by="label")
head(dat)
dat1<-dat[c(1,2,4)]
head(dat1)
tail(dat1)
vert_col<-vert+geom_hilight(data=dat1,mapping=aes(node=node,fill=LF),extendto=3,alpha=0.5)+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "bottom")
vert_col

vert1<-vert_col+new_scale_fill()+
  geom_fruit(data=sixteen,geom=geom_bar,mapping=aes(y=Species,x=OG0009063,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=fift,geom=geom_bar,mapping=aes(y=Species,x=OG0008824,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=fort,geom=geom_bar,mapping=aes(y=Species,x=OG0007224,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=thirt,geom=geom_bar,mapping=aes(y=Species,x=OG0007180,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=twel,geom=geom_bar,mapping=aes(y=Species,x=OG0006856,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=eleven,geom=geom_bar,mapping=aes(y=Species,x=OG0006823,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=tenth,geom=geom_bar,mapping=aes(y=Species,x=OG0006508,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=ninth,geom=geom_bar,mapping=aes(y=Species,x=OG0003599,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=eight,geom=geom_bar,mapping=aes(y=Species,x=OG0002730,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=seven,geom=geom_bar,mapping=aes(y=Species,x=OG0001684,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=sixth,geom=geom_bar,mapping=aes(y=Species,x=OG0001517,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=fifth,geom=geom_bar,mapping=aes(y=Species,x=OG0001290,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=forth,geom=geom_bar,mapping=aes(y=Species,x=OG0001157,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=thrd,geom=geom_bar,mapping=aes(y=Species,x=OG0000893,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=sec,geom=geom_bar,mapping=aes(y=Species,x=OG0000781,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=first,geom=geom_bar,mapping=aes(y=Species,x=OG0000701,fill=LF),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  scale_fill_brewer(palette="Dark2")+theme(axis.title.x = element_text(size=1))

vert1

list.files()
getwd()
head(p)
p1<-melt(p)
head(p1)
p_sort<-left_join(splist,p1,by="Species")
head(p_sort)
p_sort$Species<-as.character(p_sort$Species)
p_sort$Species<-factor(p_sort$Species,levels=unique(p_sort$Species))
ggplot(p_sort)+geom_bar(aes(x=value,y=Species,fill=LF),stat="identity")+
  scale_fill_brewer(palette="Dark2")+facet_grid(~variable,scale="free")
