#################################################################################################################################
#############################                                                                       #############################
#############################         Phylogenetic Principal Component Analysis for PCWDEs          #############################
#############################                                                                       #############################
#################################################################################################################################

library(phytools)
library(tibble)
library(ape)
library(ggplot2)
library(readr)
library(tidyr)
library(tibble)
library(tidyverse)
library(ggtree)
library(gridExtra)
library(dplyr)
library(reshape2)
library(ggrepel)
library(RColorBrewer)
library(plyr)


###STEP1 : Species tree fr Dataset2
##read the tree for 131 sps
tree<-read.nexus("131sps_final_rooted.nex")
plotTree(tree)
tree$tip.label


###STEP2 : Species list for Phylogenetic PCAs
##list of wood-decay Basidiomycota
##this is required for making the pruned tree- only WR, LD lifestyles and Physalacriaceae
species<-c("Aurde3_1","Exigl1","Ricmel1","Schpa1","Lopni1","Stehi1","Phlrad1","Bjead1_1","Phchr2","Pycsa1",
           "Trave1","Dicsqu464_1","Lenti7_1","Polar1","Fibsp1","Ptegra1","Pleery1","PleosPC15_2","Cligib1",
           "Hypma","Lepnud1","Amath1","Volvo1","Agabi_varbisH97_2","Macfu1","Crula1","Cyastr2","Copmic2",
           "Copci_AmutBmut1","Copmar1","Panpap1","Phoaln1","Hypsu1","Phocon1","Galma1","Gymjun1","Agrped1",
           "Psicub1_1","Psiser1","Myccro1","Mycepi1","Mycflo1","Marfi1","Gyman1","Gymlu1","Lened_B_1_1",
           "Rhobut1_1","Hymrad1","Oudmuc1","Cylto1","Flave","Flaro","Guyne1","Armtab1",
           "Armect1","Armbor1","Armost1","Armosto1","Armnab1","Armcep1","Armga1","Armmel1",
           "Armme1_1","Armfus","Armnov1","Armfum1","Armlut1")

#make a pruned tree with above species
ptree<-drop.tip(tree,tree$tip.label[-match(species, tree$tip.label)])
ptree$tip.label
plot(ptree)


###STEP3 : Phylogenetic PCA using phytools
##This example is for Cellulases
##The file shows copy numbers of Cellulase related PCWDEs in species from Dataset2
##The gene copy numbers were normalized to proteome sizes
#for other substrates, file name can be changed from "cellul.tsv" to other substrates (eg. pectin.tsv)
Cellu<-read.delim("cellu.tsv", header = T, row.names = 1)
head(Cellu)
pCellu<-Cellu%>%filter(row.names(Cellu) %in% species)
pCellu
pCellu<-pCellu[, colSums(pCellu != 0) > 0]
colnames(Cellu)
colnames(pCellu)

x2=pCellu[,colnames(pCellu)!="Lifestyle"]
head(x2)
X<-phyl.pca(ptree, x2, method="BM", mode="cov")
##to plot the phylopca, phytool uses biplot
biplot(X, var.axes=T,main="Cellu")


####for plotting phylo-PCA in ggplots
S<-as_tibble(X$S)
S$name<-rownames(X$S)
pCellu$name<-rownames(pCellu)
x2<-as_tibble(pCellu)
S2<-bind_cols(S,x2[match(S$name, x2$name),"Lifestyle"])

##calculate the percentage variance along x and y axis
diag(X$Eval)/sum(X$Eval)*100

##plotting all species names
Cellu_pca<-ggplot(S2,aes(PC1,PC2,color=Lifestyle,label=name))+geom_text(size=4)+
  theme_bw()+scale_color_manual(values = c("#1B9E77","#E7298A","#666666"))+
  scale_fill_manual(values = c("#1B9E77","#E7298A","#666666"))+
  stat_ellipse(geom="polygon", aes(fill=Lifestyle), alpha=0.1, show.legend = T,level=0.6)+
  ggtitle("Cellulases")+xlab("PC1 (49.82%)")+ylab("PC2 (18.01%)")+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20,color="steelblue"),
        legend.position = "none",axis.text=element_text(color="steelblue",size=12),
        axis.title = element_text(color="steelblue",size=15))
Cellu_pca

##plotting only Physalacriaceae species names
Cellu_pca1<-ggplot(S2,aes(PC1,PC2,color=Lifestyle,label=name))+geom_point(size=2)+
  theme_bw()+scale_color_manual(values = c("#1B9E77", "#E7298A","#666666"))+
  scale_fill_manual(values = c("#1B9E77", "#E7298A","#666666"))+
  stat_ellipse(geom="polygon", aes(fill=Lifestyle), alpha=0.1, show.legend = T,level=0.6)+
  ggtitle("Cellulases")+xlab("PC1 (49.82%)")+ylab("PC2 (18.01%)")+
  geom_text_repel(aes(label=ifelse(Lifestyle=="Ba_SR",name,""),size=4))+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20,color="steelblue"),
        legend.position = "none",axis.text=element_text(color="steelblue",size=12),
        axis.title = element_text(color="steelblue",size=15))
Cellu_pca1


##pca loading factors
load<-X$L
write.table(load,"Cellu_loadings.tsv",sep="\t",quote=F,row.names = T)
L<-as_tibble(X$L)
L$name<-rownames(X$L)
pCellu$name<-rownames(pCellu)
x3<-as_tibble(pCellu)
L2<-bind_cols(L,x3[match(L$name, x3$name),"Lifestyle"])

##list the range of x and y axis respectively from the phylo-pca
ggplot_build(Cellu_pca)$layout$panel_scales_x[[1]]$range$range
ggplot_build(Cellu_pca)$layout$panel_scales_y[[1]]$range$range


cellu_loadings<-ggplot(L2,aes(PC1,PC2,color=Lifestyle, label=name))+
  geom_segment(aes(x=0,y=0,xend=PC1,yend=PC2), arrow=arrow(length = unit(0.1,"in")),color="blue",alpha=0.5)+
  geom_text(size=5)+theme_bw()+scale_color_manual(values = c("#1B9E77","#E7298A","#666666"))

##overlay region of pca on loadings
cellu_loadings1<-cellu_loadings+ggtitle("Cellulases Loadings")+xlab("PC1")+ylab("PC2")+
  geom_rect(data=L2, aes(xmin=-0.09734199  , xmax=0.35982449, ymin=-0.09416557, ymax=0.09969399),
            color="red",fill=NA,alpha=0.9,inherit.aes = FALSE)+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20,color="steelblue"),
        legend.position = "none",axis.text=element_text(color="steelblue",size=12),
        axis.title = element_text(color="steelblue",size=15))
  
final_plot<-Cellu_pca+Cellu_loadings1

##same can be done for other substrate files


