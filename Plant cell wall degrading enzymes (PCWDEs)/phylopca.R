rm(list=ls())
setwd("C:/Users/Dell/Documents/TopGO/Armi_compgen_updated/CAZy_PCWDEs/")
list.files()

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
library(gghighlight)
library(ggnewscale)
library(ggrepel)
library(rayshader)
library(RColorBrewer)
library(plyr)
library(ggConvexHull)
library(ggpubr)
########################################################################
#######################_______CAZY figure________#######################
########################################################################
#read the tree for 132 sps
tree<-read.nexus("132sps_final_rooted.nex")
plotTree(tree)
tree$tip.label

tree<-drop.tip(tree,"Cysmur1")
tree$tip.label
write.tree(tree,"dropped_Cysmur1.tree")
#list of basidiomycetes - this is required for making the pruned tree
#species<-c("Calco1","Calvi1","Rhiso1","Botbo1","Tulca1","Aurde3_1",
           "Exigl1","Ricmel1","Schpa1","Lopni1","Hetan2","Stehi1",
           "Jaaar1","Glotr1_1","Helsul1","Neole1","Phlrad1","Bjead1_1",
           "Phchr2","Pycsa1","Trave1","Dicsqu464_1","Lenti7_1","Polar1",
           "Daequ1","Fompi3","Laesu1","Wolco1","Fibra1","PosplRSB12_1",
           "Fibsp1","Pilcr1","SerlaS7_3_2","Conol1","Conpu1","Rhives1",
           "Rhivi1","Suibr2","Suilu4","Hydpi2","Paxin1","Paxru2","Sclci1",
           "Pismi1","Pisti1","Ptegra1","Pleery1","PleosPC15_2","FlaPMI526_1",
           "Fishe1","Auramp1","Schco3","Cligib1","Hypma","Lepnud1","Amamu1",
           "Amath1","Plucer1","Volvo1","Agabi_varbisH97_2","Macfu1","Crula1",
           "Cyastr2","Lacam2","Lacbi2","Copmic2","Copci_AmutBmut1","Copmar1",
           "Panpap1","Crevar1","Hebcy2","Phoaln1","Hypsu1","Phocon1","Galma1",
           "Gymjun1","Agrped1","Psicub1_1","Psiser1","Myccro1","Mycepi1","Mycflo1",
           "Chopu","Cysmur1","Marfi1","Monpe1_1","Denbi1","Gyman1","Gymlu1","Lened_B_1_1",
           "Rhobut1_1","Hymrad1","Oudmuc1","Cylto1","Flave","Flaro","Guyne1","Armtab1",
           "Armect1","Armbor1","Armost1","Armosto1","Armnab1","Armcep1","Armga1","Armmel1",
           "Armme1_1","Armfus","Armnov1","Armfum1","Armlut1")
#############################################################################################################
#list ofwood-decay basidiomycetes - this is required for making the pruned tree- only WR,BR,LD,BaSR
#species<-c("Calco1","Calvi1","Aurde3_1","Exigl1","Ricmel1","Schpa1","Lopni1","Stehi1",
           "Glotr1_1","Helsul1","Neole1","Phlrad1","Bjead1_1",
           "Phchr2","Pycsa1","Trave1","Dicsqu464_1","Lenti7_1","Polar1",
           "Daequ1","Fompi3","Laesu1","Wolco1","Fibra1","PosplRSB12_1",
           "Fibsp1","SerlaS7_3_2","Conol1","Conpu1","Hydpi2","Ptegra1","Pleery1","PleosPC15_2",
           "Fishe1","Cligib1","Hypma","Lepnud1",
           "Amath1","Volvo1","Agabi_varbisH97_2","Macfu1","Crula1",
           "Cyastr2","Copmic2","Copci_AmutBmut1","Copmar1",
           "Panpap1","Phoaln1","Hypsu1","Phocon1","Galma1",
           "Gymjun1","Agrped1","Psicub1_1","Psiser1","Myccro1","Mycepi1","Mycflo1",
           "Cysmur1","Marfi1","Gyman1","Gymlu1","Lened_B_1_1",
           "Rhobut1_1","Hymrad1","Oudmuc1","Cylto1","Flave","Flaro","Guyne1","Armtab1",
           "Armect1","Armbor1","Armost1","Armosto1","Armnab1","Armcep1","Armga1","Armmel1",
           "Armme1_1","Armfus","Armnov1","Armfum1","Armlut1")

####___USE THIS #############
####___USE THIS #############
####___USE THIS #############
#list ofwood-decay basidiomycetes - this is required for making the pruned tree- only WR,LD,BaSR - without BR
species<-c("Aurde3_1","Exigl1","Ricmel1","Schpa1","Lopni1","Stehi1",
           "Phlrad1","Bjead1_1",
           "Phchr2","Pycsa1","Trave1","Dicsqu464_1","Lenti7_1","Polar1",
           "Fibsp1","Ptegra1","Pleery1","PleosPC15_2",
           "Cligib1","Hypma","Lepnud1",
           "Amath1","Volvo1","Agabi_varbisH97_2","Macfu1","Crula1",
           "Cyastr2","Copmic2","Copci_AmutBmut1","Copmar1",
           "Panpap1","Phoaln1","Hypsu1","Phocon1","Galma1",
           "Gymjun1","Agrped1","Psicub1_1","Psiser1","Myccro1","Mycepi1","Mycflo1",
           "Marfi1","Gyman1","Gymlu1","Lened_B_1_1",
           "Rhobut1_1","Hymrad1","Oudmuc1","Cylto1","Flave","Flaro","Guyne1","Armtab1",
           "Armect1","Armbor1","Armost1","Armosto1","Armnab1","Armcep1","Armga1","Armmel1",
           "Armme1_1","Armfus","Armnov1","Armfum1","Armlut1")
#"Cysmur1" REMOVE Cysmur1
#make a pruned tree for only basidiomycetes

ptree<-drop.tip(tree,tree$tip.label[-match(species, tree$tip.label)])
ptree$tip.label
plot(ptree)
#physac<-c("Hymrad1","Oudmuc1","Cylto1","Flave","Flaro","Guyne1","Armtab1",
#"Armect1","Armbor1","Armost1","Armosto1","Armnab1","Armcep1","Armga1","Armmel1",
#"Armme1_1","Armfus","Armnov1","Armfum1","Armlut1")
#ptree<-drop.tip(tree,tree$tip.label[-match(physac, tree$tip.label)])
#plotTree(ptree)
#write.tree(ptree, file = "physac_tree.tree", append = FALSE)

#################################################################################################################
#check and customize color for species lifestyle
#display.brewer.pal(n = 8, name = 'Dark2') #shows the color palette
#brewer.pal(n = 8, name = "Dark2") #shows the name of the colors
#1B9E77 BaSR # for phylopca
#D95F02 BR # for phylopca
#7570B3 ECM
#E7298A LD # for phylopca
#66A61E PATH
#E6AB02 SR
#A6761D UR
#666666 WR # for phylopca

#phylopca for pruned species tree i.e. only basidios without BR

#Cellu_expn
Cellu_expn<-read.delim("Cellu_expn.tsv", header = T, row.names = 1)
head(Cellu_expn)
pCellu_expn<-Cellu_expn%>%filter(row.names(Cellu_expn) %in% species)
pCellu_expn
pCellu_expn<-pCellu_expn[, colSums(pCellu_expn != 0) > 0]
colnames(Cellu_expn)
colnames(pCellu_expn)

x2=pCellu_expn[,colnames(pCellu_expn)!="Lifestyle"]
head(x2)
X<-phyl.pca(ptree, x2, method="BM", mode="cov")
biplot(X, var.axes=T,main="Cellu_expn_w/o_BR")
biplot(X, var.axes=T,main="Cellu_expn w/o BR",ylim=c(-0.1,0.15),col=c("grey","blue"))

####for species PCA
S<-as_tibble(X$S)
S$name<-rownames(X$S)
pCellu_expn$name<-rownames(pCellu_expn)
x2<-as_tibble(pCellu_expn)
S2<-bind_cols(S,x2[match(S$name, x2$name),"Lifestyle"])
#with_clouds
diag(X$Eval)/sum(X$Eval)*100
Cellu_expn_pruned<-ggplot(S2,aes(PC1,PC2,color=Lifestyle,label=name))+geom_text(size=4)+
  theme_bw()+scale_color_manual(values = c("#1B9E77","#E7298A","#666666"))+
  scale_fill_manual(values = c("#1B9E77","#E7298A","#666666"))+
  stat_ellipse(geom="polygon", aes(fill=Lifestyle), alpha=0.1, show.legend = T,level=0.6)+
  ggtitle("Cellulases")+xlab("PC1 (49.82%)")+ylab("PC2 (18.01%)")+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20,color="steelblue"),
        legend.position = "none",axis.text=element_text(color="steelblue",size=12),
        axis.title = element_text(color="steelblue",size=15))
Cellu_expn_pruned

Cellu_expn_pruned1<-ggplot(S2,aes(PC1,PC2,color=Lifestyle,label=name))+geom_point(size=2)+
  theme_bw()+scale_color_manual(values = c("#1B9E77", "#E7298A","#666666"))+
  scale_fill_manual(values = c("#1B9E77", "#E7298A","#666666"))+
  stat_ellipse(geom="polygon", aes(fill=Lifestyle), alpha=0.1, show.legend = T,level=0.6)+
  ggtitle("Cellulases")+xlab("PC1 (49.82%)")+ylab("PC2 (18.01%)")+
  geom_text_repel(aes(label=ifelse(Lifestyle=="Ba_SR",name,""),size=4))+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20,color="steelblue"),
        legend.position = "none",axis.text=element_text(color="steelblue",size=12),
        axis.title = element_text(color="steelblue",size=15))
Cellu_expn_pruned1


######for loadings
load<-X$L
write.table(load,"Cellu_expn_wo_BR_loadings.tsv",sep="\t",quote=F,row.names = T)
L<-as_tibble(X$L)
L$name<-rownames(X$L)
pCellu_expn$name<-rownames(pCellu_expn)
x3<-as_tibble(pCellu_expn)
L2<-bind_cols(L,x3[match(L$name, x3$name),"Lifestyle"])
p<-ggplot(L2,aes(PC1,PC2,color=Lifestyle, label=name))+
  geom_segment(aes(x=0,y=0,xend=PC1,yend=PC2), arrow=arrow(length = unit(0.1,"in")),color="blue",alpha=0.5)+
  geom_text(size=5)+theme_bw()+scale_color_manual(values = c("#1B9E77","#E7298A","#666666"))
Cellu_expn_loadings<-p+ggtitle("Cellu_expn")+xlab("PC1")+ylab("PC2")
Cellu_expn_loadings


##list the range of x and y axis respectively
ggplot_build(Cellu_expn_pruned)$layout$panel_scales_x[[1]]$range$range
ggplot_build(Cellu_expn_pruned)$layout$panel_scales_y[[1]]$range$range

Cellu_expn_loadings<-p+ggtitle("Cellulases Loadings")+xlab("PC1")+ylab("PC2")+
  geom_rect(data=L2, aes(xmin=-0.09734199  , xmax=0.35982449, ymin=-0.09416557, ymax=0.09969399),
            color="red",fill=NA,alpha=0.9,inherit.aes = FALSE)+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20,color="steelblue"),
        legend.position = "none",axis.text=element_text(color="steelblue",size=12),
        axis.title = element_text(color="steelblue",size=15))
  
c1<-Cellu_expn_pruned+Cellu_expn_loadings
c1


#Hemicellulases
Hemicellu<-read.delim("Hemicellu.tsv", header = T, row.names = 1)
head(Hemicellu)
pHemicellu<-Hemicellu%>%filter(row.names(Hemicellu) %in% species)
pHemicellu
pHemicellu<-pHemicellu[, colSums(pHemicellu != 0) > 0]
colnames(Hemicellu)
colnames(pHemicellu)

x2=pHemicellu[,colnames(pHemicellu)!="Lifestyle"]
head(x2)
X<-phyl.pca(ptree, x2, method="BM", mode="cov")
biplot(X, var.axes=T,main="Hemicellulases_w/o_BR")
biplot(X, var.axes=T,main="Hemicellulases w/o BR",col=c("grey","blue"))
X$L #loadings
X$S

####for species PCA
S<-as_tibble(X$S)
S$name<-rownames(X$S)
pHemicellu$name<-rownames(pHemicellu)
x2<-as_tibble(pHemicellu)
S2<-bind_cols(S,x2[match(S$name, x2$name),"Lifestyle"])
p<-ggplot(S2,aes(PC1,PC2,color=Lifestyle, label=name))+
  geom_text(size=5)+theme_bw()+scale_color_manual(values = c("#1B9E77","#E7298A","#666666"))
p1<-p+ggtitle("Hemicellulases")+xlab("PC1")+ylab("PC2")
p1

#with_clouds
diag(X$Eval)/sum(X$Eval)*100
Hemicellu_pruned<-ggplot(S2,aes(PC1,PC2,color=Lifestyle,label=name))+geom_text(size=4)+
  theme_bw()+scale_color_manual(values = c("#1B9E77","#E7298A","#666666"))+
  scale_fill_manual(values = c("#1B9E77","#E7298A","#666666"))+
  stat_ellipse(geom="polygon", aes(fill=Lifestyle), alpha=0.1, show.legend = T,level=0.6)+
  ggtitle("Hemicellulases")+xlab("PC1 (43.46%)")+ylab("PC2 (20.46%)")+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20,color="steelblue"),
        legend.position = "none",axis.text=element_text(color="steelblue",size=12),
        axis.title = element_text(color="steelblue",size=15))
Hemicellu_pruned

Hemicellu_pruned1<-ggplot(S2,aes(PC1,PC2,color=Lifestyle,label=name))+geom_point(size=2)+
  theme_bw()+scale_color_manual(values = c("#1B9E77", "#E7298A","#666666"))+
  scale_fill_manual(values = c("#1B9E77", "#E7298A","#666666"))+
  stat_ellipse(geom="polygon", aes(fill=Lifestyle), alpha=0.1, show.legend = T,level=0.6)+
  ggtitle("Hemicellulases")+xlab("PC1 (43.46%)")+ylab("PC2 (20.46%)")+ 
  geom_text_repel(aes(label=ifelse(Lifestyle=="Ba_SR",name,""),size=4))+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20,color="steelblue"),
        legend.position = "none",axis.text=element_text(color="steelblue",size=12),
        axis.title = element_text(color="steelblue",size=15))
Hemicellu_pruned1


######for loadings
load<-X$L
write.table(load,"Hemicellu_wo_BR_loadings.tsv",sep="\t",quote=F,row.names = T)
L<-as_tibble(X$L)
L$name<-rownames(X$L)
pHemicellu$name<-rownames(pHemicellu)
x3<-as_tibble(pHemicellu)
L2<-bind_cols(L,x3[match(L$name, x3$name),"Lifestyle"])
p<-ggplot(L2,aes(PC1,PC2,color=Lifestyle, label=name))+
  geom_segment(aes(x=0,y=0,xend=PC1,yend=PC2), arrow=arrow(length = unit(0.1,"in")),color="blue",alpha=0.5)+
  geom_text(size=5)+theme_bw()+scale_color_manual(values = c("#1B9E77","#E7298A","#666666"))
Hemicellu_loadings<-p+ggtitle("Hemicellulases")+xlab("PC1")+ylab("PC2")
Hemicellu_loadings


##list the range of x and y axis respectively
ggplot_build(Hemicellu_pruned)$layout$panel_scales_x[[1]]$range$range
ggplot_build(Hemicellu_pruned)$layout$panel_scales_y[[1]]$range$range

Hemicellu_loadings<-p+ggtitle("Hemicellulases Loadings")+xlab("PC1")+ylab("PC2")+
  geom_rect(data=L2, aes(xmin=-0.09018479, xmax=0.06844143, ymin=-0.09375855, ymax=0.04590103),
            color="red",fill=NA,alpha=0.9,inherit.aes = FALSE)+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20,color="steelblue"),
        legend.position = "none",axis.text=element_text(color="steelblue",size=12),
        axis.title = element_text(color="steelblue",size=15))
c2<-Hemicellu_pruned+Hemicellu_loadings
c2


#Pectinases
Pectin<-read.delim("Pectin.tsv", header = T, row.names = 1)
head(Pectin)
pPectin<-Pectin%>%filter(row.names(Pectin) %in% species)
pPectin
pPectin<-pPectin[, colSums(pPectin != 0) > 0]
colnames(Pectin)
colnames(pPectin)

x2=pPectin[,colnames(pPectin)!="Lifestyle"]
head(x2)
X<-phyl.pca(ptree, x2, method="BM", mode="cov")
biplot(X, var.axes=T,main="Pectinases_w/o_BR")
biplot(X, var.axes=T,main="Pectinases w/o BR",col=c("grey","blue"))
X$L #loadings
X$S

####for species PCA
S<-as_tibble(X$S)
S$name<-rownames(X$S)
pPectin$name<-rownames(pPectin)
x2<-as_tibble(pPectin)
S2<-bind_cols(S,x2[match(S$name, x2$name),"Lifestyle"])
p<-ggplot(S2,aes(PC1,PC2,color=Lifestyle, label=name))+
  geom_text(size=5)+theme_bw()+scale_color_manual(values = c("#1B9E77","#E7298A","#666666"))
p1<-p+ggtitle("Pectinases")+xlab("PC1")+ylab("PC2")
p1

#with_clouds
diag(X$Eval)/sum(X$Eval)*100

Pectin_pruned<-ggplot(S2,aes(PC1,PC2,color=Lifestyle,label=name))+geom_text(size=4)+
  theme_bw()+scale_color_manual(values = c("#1B9E77","#E7298A","#666666"))+
  scale_fill_manual(values = c("#1B9E77","#E7298A","#666666"))+
  stat_ellipse(geom="polygon", aes(fill=Lifestyle), alpha=0.1, show.legend = T,level=0.6)+
  ggtitle("Pectinases")+xlab("PC1 (44.02%)")+ylab("PC2 (28.64%)")+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20,color="steelblue"),
        legend.position = "none",axis.text=element_text(color="steelblue",size=12),
        axis.title = element_text(color="steelblue",size=15))
  
Pectin_pruned

Pectin_pruned1<-ggplot(S2,aes(PC1,PC2,color=Lifestyle,label=name))+geom_point(size=2)+
  theme_bw()+scale_color_manual(values = c("#1B9E77", "#E7298A","#666666"))+
  scale_fill_manual(values = c("#1B9E77", "#E7298A","#666666"))+
  stat_ellipse(geom="polygon", aes(fill=Lifestyle), alpha=0.1, show.legend = T,level=0.6)+
  ggtitle("Pectinases")+xlab("PC1 (44.02%)")+ylab("PC2 (28.64%)")+ 
  geom_text_repel(aes(label=ifelse(Lifestyle=="Ba_SR",name,""),size=4))+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20,color="steelblue"),legend.position = "none",axis.text=element_text(color="steelblue",size=12),         axis.title = element_text(color="steelblue",size=15))
Pectin_pruned1


######for loadings
load<-X$L
write.table(load,"Pectin_wo_BR_loadings.tsv",sep="\t",quote=F,row.names = T)
L<-as_tibble(X$L)
L$name<-rownames(X$L)
pPectin$name<-rownames(pPectin)
x3<-as_tibble(pPectin)
L2<-bind_cols(L,x3[match(L$name, x3$name),"Lifestyle"])
p<-ggplot(L2,aes(PC1,PC2,color=Lifestyle, label=name))+
  geom_segment(aes(x=0,y=0,xend=PC1,yend=PC2), arrow=arrow(length = unit(0.1,"in")),color="blue",alpha=0.5)+
  geom_text(size=5)+theme_bw()+scale_color_manual(values = c("#1B9E77","#E7298A","#666666"))
Pectin_loadings<-p+ggtitle("Pectinases")+xlab("PC1")+ylab("PC2")
Pectin_loadings


##list the range of x and y axis respectively
ggplot_build(Pectin_pruned)$layout$panel_scales_x[[1]]$range$range
ggplot_build(Pectin_pruned)$layout$panel_scales_y[[1]]$range$range

Pectin_loadings<-p+ggtitle("Pectinases Loadings")+xlab("PC1")+ylab("PC2")+
  geom_rect(data=L2, aes(xmin=-0.08733975, xmax=0.08478027, ymin=-0.03962775, ymax=0.12566630),
            color="red",fill=NA,alpha=0.9,inherit.aes = FALSE)+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20,color="steelblue"),         legend.position = "none",axis.text=element_text(color="steelblue",size=12),         axis.title = element_text(color="steelblue",size=15))

c3<-Pectin_pruned+Pectin_loadings
c3

#Ligninases
Lignin<-read.delim("Lignin.tsv", header = T, row.names = 1)
head(Lignin)
pLignin<-Lignin%>%filter(row.names(Lignin) %in% species)
pLignin
pLignin<-pLignin[, colSums(pLignin != 0) > 0]
colnames(Lignin)
colnames(pLignin)

x2=pLignin[,colnames(pLignin)!="Lifestyle"]
head(x2)
X<-phyl.pca(ptree, x2, method="BM", mode="cov")
biplot(X, var.axes=T,main="Ligninases_w/o_BR")
biplot(X, var.axes=T,main="Ligninases w/o BR",col=c("grey","blue"))

####for species PCA
S<-as_tibble(X$S)
S$name<-rownames(X$S)
pLignin$name<-rownames(pLignin)
x2<-as_tibble(pLignin)
S2<-bind_cols(S,x2[match(S$name, x2$name),"Lifestyle"])
p<-ggplot(S2,aes(PC1,PC2,color=Lifestyle, label=name))+
  geom_text(size=5)+theme_bw()+scale_color_manual(values = c("#1B9E77","#E7298A","#666666"))
p1<-p+ggtitle("Ligninases")+xlab("PC1")+ylab("PC2")
p1

#with_clouds
diag(X$Eval)/sum(X$Eval)*100
Lignin_pruned<-ggplot(S2,aes(PC1,PC2,color=Lifestyle,label=name))+geom_text(size=4)+
  theme_bw()+scale_color_manual(values = c("#1B9E77","#E7298A","#666666"))+
  scale_fill_manual(values = c("#1B9E77","#E7298A","#666666"))+
  stat_ellipse(geom="polygon", aes(fill=Lifestyle), alpha=0.1, show.legend = T,level=0.6)+
  ggtitle("Ligninases")+xlab("PC1 (59.98%)")+ylab("PC2 (34.47%)")+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20,color="steelblue"),         legend.position = "none",axis.text=element_text(color="steelblue",size=12),         axis.title = element_text(color="steelblue",size=15))
Lignin_pruned

Lignin_pruned1<-ggplot(S2,aes(PC1,PC2,color=Lifestyle,label=name))+geom_point(size=2)+
  theme_bw()+scale_color_manual(values = c("#1B9E77", "#E7298A","#666666"))+
  scale_fill_manual(values = c("#1B9E77", "#E7298A","#666666"))+
  stat_ellipse(geom="polygon", aes(fill=Lifestyle), alpha=0.1, show.legend = T,level=0.6)+
  ggtitle("Ligninases")+xlab("PC1 (59.98%)")+ylab("PC2 (34.47%)")+ 
  geom_text_repel(aes(label=ifelse(Lifestyle=="Ba_SR",name,""),size=4))+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20,color="steelblue"),         legend.position = "none",axis.text=element_text(color="steelblue",size=12),         axis.title = element_text(color="steelblue",size=15))
Lignin_pruned1


######for loadings
load<-X$L
write.table(load,"Lignin_wo_BR_loadings.tsv",sep="\t",quote=F,row.names = T)
L<-as_tibble(X$L)
L$name<-rownames(X$L)
pLignin$name<-rownames(pLignin)
x3<-as_tibble(pLignin)
L2<-bind_cols(L,x3[match(L$name, x3$name),"Lifestyle"])
p<-ggplot(L2,aes(PC1,PC2,color=Lifestyle, label=name))+
  geom_segment(aes(x=0,y=0,xend=PC1,yend=PC2), arrow=arrow(length = unit(0.1,"in")),color="blue",alpha=0.5)+
  geom_text(size=5)+theme_bw()+scale_color_manual(values = c("#1B9E77","#E7298A","#666666"))
Lignin_loadings<-p+ggtitle("Ligninases")+xlab("PC1")+ylab("PC2")
Lignin_loadings


##list the range of x and y axis respectively
ggplot_build(Lignin_pruned)$layout$panel_scales_x[[1]]$range$range
ggplot_build(Lignin_pruned)$layout$panel_scales_y[[1]]$range$range

Lignin_loadings<-p+ggtitle("Ligninases_Loadings")+xlab("PC1")+ylab("PC2")+
  geom_rect(data=L2, aes(xmin=-0.05282943, xmax=0.14357546, ymin=-0.10411443, ymax=0.09144819),
            color="red",fill=NA,alpha=0.9,inherit.aes = FALSE)+theme(plot.title = element_text(hjust = 0.5,face="bold",size=20,color="steelblue"),         legend.position = "none",axis.text=element_text(color="steelblue",size=12),         axis.title = element_text(color="steelblue",size=15))

c4<-Lignin_pruned+Lignin_loadings
c4

#Putative_ligninases
Putative_lignin<-read.delim("Putative_lignin.tsv", header = T, row.names = 1)
head(Putative_lignin)
pPutative_lignin<-Putative_lignin%>%filter(row.names(Putative_lignin) %in% species)
pPutative_lignin
pPutative_lignin<-pPutative_lignin[, colSums(pPutative_lignin != 0) > 0]
colnames(Putative_lignin)
colnames(pPutative_lignin)

x2=pPutative_lignin[,colnames(pPutative_lignin)!="Lifestyle"]
head(x2)
X<-phyl.pca(ptree, x2, method="BM", mode="cov")
biplot(X, var.axes=T,main="Putative_ligninases_w/o_BR")
biplot(X, var.axes=T,main="Putative_ligninases w/o BR",ylim=c(-0.08,0.04),col=c("grey","blue"))

####for species PCA
S<-as_tibble(X$S)
S$name<-rownames(X$S)
pPutative_lignin$name<-rownames(pPutative_lignin)
x2<-as_tibble(pPutative_lignin)
S2<-bind_cols(S,x2[match(S$name, x2$name),"Lifestyle"])

#with_clouds
diag(X$Eval)/sum(X$Eval)*100

Putative_lignin_pruned<-ggplot(S2,aes(PC1,PC2,color=Lifestyle,label=name))+geom_text(size=4)+
  theme_bw()+scale_color_manual(values = c("#1B9E77","#E7298A","#666666"))+
  scale_fill_manual(values = c("#1B9E77","#E7298A","#666666"))+
  stat_ellipse(geom="polygon", aes(fill=Lifestyle), alpha=0.1, show.legend = T,level=0.6)+
  ggtitle("Putative_ligninases")+xlab("PC1 (82.94%)")+ylab("PC2 (6.46%)")+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20,color="steelblue"),
        legend.position = "none",axis.text=element_text(color="steelblue",size=12),
        axis.title = element_text(color="steelblue",size=15))
Putative_lignin_pruned
Putative_lignin_pruned1<-ggplot(S2,aes(PC1,PC2,color=Lifestyle,label=name))+geom_point(size=2)+
  theme_bw()+scale_color_manual(values = c("#1B9E77", "#E7298A","#666666"))+
  scale_fill_manual(values = c("#1B9E77", "#E7298A","#666666"))+
  stat_ellipse(geom="polygon", aes(fill=Lifestyle), alpha=0.1, show.legend = T,level=0.6)+
  ggtitle("Putative_ligninases")+xlab("PC1 (82.94%)")+ylab("PC2 (6.45%)")+ 
  geom_text(aes(label=ifelse(Lifestyle=="Ba_SR",name,""),size=4))+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20,color="steelblue"),
        legend.position = "none",axis.text=element_text(color="steelblue",size=12),
        axis.title = element_text(color="steelblue",size=15))
Putative_lignin_pruned1


######for loadings
load<-X$L
write.table(load,"Putative_lignin_wo_BR_loadings.tsv",sep="\t",quote=F,row.names = T)
L<-as_tibble(X$L)
L$name<-rownames(X$L)
pPutative_lignin$name<-rownames(pPutative_lignin)
x3<-as_tibble(pPutative_lignin)
L2<-bind_cols(L,x3[match(L$name, x3$name),"Lifestyle"])
p<-ggplot(L2,aes(PC1,PC2,color=Lifestyle, label=name))+
  geom_segment(aes(x=0,y=0,xend=PC1,yend=PC2), arrow=arrow(length = unit(0.1,"in")),color="blue",alpha=0.5)+
  geom_text(size=5)+theme_bw()+scale_color_manual(values = c("#1B9E77","#E7298A","#666666"))
Putative_lignin_loadings<-p+ggtitle("Putative_ligninases")+xlab("PC1")+ylab("PC2")
Putative_lignin_loadings


##list the range of x and y axis respectively
ggplot_build(Putative_lignin_pruned)$layout$panel_scales_x[[1]]$range$range
ggplot_build(Putative_lignin_pruned)$layout$panel_scales_y[[1]]$range$range

Putative_lignin_loadings<-p+ggtitle("Putative_ligninases_Loadings")+xlab("PC1")+ylab("PC2")+
  geom_rect(data=L2, aes(xmin=-0.2773609, xmax=0.1036482, ymin=-0.05397133, ymax=0.04380062),
            color="red",fill=NA,alpha=0.9,inherit.aes = FALSE)+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20,color="steelblue"),         legend.position = "none",axis.text=element_text(color="steelblue",size=12),         axis.title = element_text(color="steelblue",size=15))

c5<-Putative_lignin_pruned+Putative_lignin_loadings
c5

ggarrange(c1,c3,c2,c4,c5,ncol=2,nrow=3) #saved at 20*15
ggarrange(Cellu_expn_pruned1,Pectin_pruned1,Hemicellu_pruned1,Lignin_pruned1) #saved at 16*16


####co-enriched PCWDEs
co<-read.delim("Co-shared_OGs_CN.tsv")
head(co)
co1<-melt(co)
head(co1)
co1<-co1[which(co1$Species!="Cysmur1"),]

unique(co1$LF)
levs<-c("BaSR","WR","LD","SR_Asco","BR","UR","ECM","PATH")
ggplot(co1,aes(x=factor(LF,levs),y=value,fill=LF))+geom_boxplot(alpha=0.5)+facet_wrap(variable ~.,scales="free")+
  scale_fill_brewer(palette = "Dark2")+theme_bw()

levs1<-c("BaSR","WR","LD","SR_Asco")
co2<-filter(co1,LF %in% levs1)
head(co2)
ggplot(co2,aes(x=factor(LF,levs1),y=value,fill=LF))+geom_boxplot(alpha=0.5)+facet_wrap(variable ~.,scales="free")+
  scale_fill_manual(values = c("#1B9E77","#E7298A","#E6AB02","#666666"))+theme_bw()

box<-ggplot(co2,aes(x=factor(LF,levs1),y=value,fill=LF))+geom_boxplot(alpha=0.5)+facet_grid(~variable)+
  scale_fill_manual(values = c("#1B9E77","#E7298A","#E6AB02","#666666"),labels=c("Physalacriaceae","LD","ASCO","WR"))+theme_bw()+
  theme(axis.text.x = element_blank(),#panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.text.y=element_text(color="steelblue",size=12),axis.ticks.x=element_blank(),
        strip.text = element_text(color="steelblue",size=12,face="bold"),
        axis.title = element_text(color="steelblue",size=15),legend.position = "bottom")+
  xlab("")+ylab("Number of gene copies")+
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12,14),limits=c(0,14))
box

pca<-ggarrange(Cellu_expn_pruned1,Pectin_pruned1,Hemicellu_pruned1,Lignin_pruned1,nrow=1) #saved at 16*16
pca
ggarrange(pca,box,nrow=2,heights=c(1,0.8), labels = c("A","B")) #saved at 18.02*8.36

#Suberin
Suberin<-read.delim("Suberin.tsv", header = T, row.names = 1)
head(Suberin)
pSuberin<-Suberin%>%filter(row.names(Suberin) %in% species)
pSuberin
pSuberin<-pSuberin[, colSums(pSuberin != 0) > 0]
colnames(Suberin)
colnames(pSuberin)

x2=pSuberin[,colnames(pSuberin)!="Lifestyle"]
head(x2)
X<-phyl.pca(ptree, x2, method="BM", mode="cov")
biplot(X, var.axes=T,main="Suberin_w/o_BR")
biplot(X, var.axes=T,main="Suberin w/o BR",ylim=c(-0.05,0.01),col=c("grey","blue"))

####for species PCA
S<-as_tibble(X$S)
S$name<-rownames(X$S)
pSuberin$name<-rownames(pSuberin)
x2<-as_tibble(pSuberin)
S2<-bind_cols(S,x2[match(S$name, x2$name),"Lifestyle"])

#with_clouds
Suberin_pruned<-ggplot(S2,aes(PC1,PC2,color=Lifestyle,label=name))+geom_text(size=4)+
  theme_bw()+scale_color_manual(values = c("#1B9E77","#E7298A","#666666"))+
  scale_fill_manual(values = c("#1B9E77","#E7298A","#666666"))+
  stat_ellipse(geom="polygon", aes(fill=Lifestyle), alpha=0.1, show.legend = T,level=0.6)+
  ggtitle("Suberin")+xlab("PC1")+ylab("PC2")
Suberin_pruned
Suberin_pruned1<-ggplot(S2,aes(PC1,PC2,color=Lifestyle,label=name))+geom_point(size=2)+
  theme_bw()+scale_color_manual(values = c("#1B9E77", "#E7298A","#666666"))+
  scale_fill_manual(values = c("#1B9E77", "#E7298A","#666666"))+
  stat_ellipse(geom="polygon", aes(fill=Lifestyle), alpha=0.1, show.legend = T,level=0.6)+
  ggtitle("Suberin")+xlab("PC1")+ylab("PC2")+ 
  geom_text(aes(label=ifelse(Lifestyle=="Ba_SR",name,""),size=4))+
  theme(legend.position = "none")
Suberin_pruned1


######for loadings
load<-X$L
write.table(load,"Suberin_wo_BR_loadings.tsv",sep="\t",quote=F,row.names = T)
L<-as_tibble(X$L)
L$name<-rownames(X$L)
pSuberin$name<-rownames(pSuberin)
x3<-as_tibble(pSuberin)
L2<-bind_cols(L,x3[match(L$name, x3$name),"Lifestyle"])
p<-ggplot(L2,aes(PC1,PC2,color=Lifestyle, label=name))+
  geom_segment(aes(x=0,y=0,xend=PC1,yend=PC2), arrow=arrow(length = unit(0.1,"in")),color="blue",alpha=0.5)+
  geom_text(size=5)+theme_bw()+scale_color_manual(values = c("#1B9E77","#E7298A","#666666"))
Suberin_loadings<-p+ggtitle("Suberin")+xlab("PC1")+ylab("PC2")
Suberin_loadings


##list the range of x and y axis respectively
ggplot_build(Suberin_pruned)$layout$panel_scales_x[[1]]$range$range
ggplot_build(Suberin_pruned)$layout$panel_scales_y[[1]]$range$range

Suberin_loadings<-p+ggtitle("Suberin_Loadings")+xlab("PC1")+ylab("PC2")+
  geom_rect(data=L2, aes(xmin=-0.01026302, xmax=0.07977887, ymin=-0.04771960, ymax=0.01281996),
            color="red",fill=NA,alpha=0.9,inherit.aes = FALSE)

c7<-Suberin_pruned1+Suberin_loadings
c7
##saved pdf at 16*8

############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
###co-expanded IPRs
ipr<-read.delim("ipr_all_phylobar_uniq.tsv", header = T, row.names = 1)
head(ipr)
rownames(ipr)
colnames(ipr)

x2=ipr[,colnames(ipr)!="Lifestyle"]
head(x2)
X<-phyl.pca(tree, x2, method="BM", mode="cov")
biplot(X, var.axes=T,main="ipr")
#biplot(X, var.axes=T,main="ipr",ylim=c(-0.08,0.09))
#biplot(X, var.axes=T,main="ipr",ylim=c(-0.08,0.09),xlim=c(-0.01,0.12),col=c("grey","blue"))
X$L #loadings
X$S

####for species PCA
S<-as_tibble(X$S)
S$name<-rownames(X$S)
ipr$name<-rownames(ipr)
x2<-as_tibble(ipr)
S2<-bind_cols(S,x2[match(S$name, x2$name),"Lifestyle"])
p<-ggplot(S2,aes(PC1,PC2,color=Lifestyle, label=name))+
  geom_text(size=5)+theme_bw()+scale_color_brewer(palette = "Dark2")
p1<-p+ggtitle("IPR")+xlab("PC1")+ylab("PC2")
p1
#with_clouds
ipr1<-ggplot(S2,aes(PC1,PC2,color=Lifestyle,label=name))+geom_text(size=4)+
  theme_bw()+scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  stat_ellipse(geom="polygon", aes(fill=Lifestyle), alpha=0.1, show.legend = T,level=0.7)+
  ggtitle("IPR")+xlab("PC1")+ylab("PC2")
ipr1
ipr2<-ggplot(S2,aes(PC1,PC2,color=Lifestyle,label=name))+geom_point(size=2)+
  theme_bw()+scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  stat_ellipse(geom="polygon", aes(fill=Lifestyle), alpha=0.1, show.legend = T,level=0.7)+
  ggtitle("Co-expanded IPRs")+xlab("PC1")+ylab("PC2")+ 
  geom_text(aes(label=ifelse(Lifestyle=="Ba_SR",name,""),size=4))+
  theme(legend.position = "none")
ipr2


######for loadings
X$L
L<-as_tibble(X$L)
L$name<-rownames(X$L)
ipr$name<-rownames(ipr)
x3<-as_tibble(ipr)
L2<-bind_cols(L,x3[match(L$name, x3$name),"Lifestyle"])
p<-ggplot(L2,aes(PC1,PC2,color=Lifestyle, label=name))+
  geom_text(size=2)+theme_bw()+scale_color_manual(values = c("#1B9E77", "#D95F02","#E7298A","#666666"))
ipr_loadings<-p+ggtitle("IPR_loadings")+xlab("PC1")+ylab("PC2")
ipr_loadings

ipr2+ipr_loadings
##list the range of x and y axis respectively
ggplot_build(ipr2)$layout$panel_scales_x[[1]]$range$range
ggplot_build(ipr2)$layout$panel_scales_y[[1]]$range$range

ipr_load_spec<-p+ggtitle("IPR_Loadings")+xlab("PC1")+ylab("PC2")+
  xlim(-868.9077,163.0451)+ylim(-353.2816,968.6956)
ipr_load_spec

ipr2_loadings<-p+ggtitle("ligninases_Loadings")+xlab("PC1")+ylab("PC2")+
  geom_rect(data=L2, aes(xmin=-868.9077, xmax=163.0451, ymin=-353.2816, ymax=968.6956),
            color="red",fill=NA,alpha=0.9,inherit.aes = FALSE)

ipr2+ipr_loadings
ipr2+ipr2_loadings

######cazy enrichment#########
x<-read.delim("Enrichfig_BaSRvsLDWR.tsv", header = T)
head(x)
tail(x)
plot<-ggplot(x,aes(x=Geneperc, y=fct_reorder(Name,Geneperc)))+
  geom_point(aes(size=Geneperc, color = FDR))+
  theme_bw(base_size=10)+
  scale_color_gradient2(limits=c(0,0.05), low="cadetblue3", mid = "yellowgreen", high="orange", midpoint = 0.025)+
  ylab(NULL) + ggtitle("CAZy_enrichment")+theme_minimal()+facet_grid(.~Type)+
  gghighlight(Text=="two",
              unhighlighted_params = list(colour=NULL,alpha=0.3),calculate_per_facet = T)
plot
head(x)
plot<-ggplot(x,aes(x=Geneperc, y=fct_reorder2(Name,Type=="BaSR_LDWR",Geneperc,.desc = F)))+
  geom_point(aes(size=Geneperc, color = FDR))+
  geom_segment(aes(xend=0, yend=Name,color=FDR,size=2))+
  theme_minimal(base_size = 10)+
  scale_color_gradient2(limits=c(0,0.05), low="cadetblue3", mid = "yellowgreen", 
                        high="orange", midpoint = 0.025)+
  ylab(NULL) + ggtitle("CAZy enrichment")+facet_grid(.~Type)+
  gghighlight(Text=="two",
              unhighlighted_params = list(colour=NULL,alpha=0.3),calculate_per_facet = T)

plot


####SR and BaSR common OGIDs##########

plotTree(tree)
splist<-data.frame(tree$tip.label)
head(splist)
colnames(splist)<-c("Species") #order of species to be used
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

data<-ggtree(tree)
data$data
vert<-ggtree(tree)+geom_tiplab(size=2,align=T)
vert


vert1<-vert+geom_fruit(data=first,geom=geom_bar,mapping=aes(y=Species,x=OG0000701,fill=Lifestyle,alpha=0.1),pwidth=1,orientation="y",stat="identity")+scale_fill_brewer(palette="Dark2")+
  geom_fruit(data=sec,geom=geom_bar,mapping=aes(y=Species,x=OG0000781,fill=Lifestyle,alpha=0.1),pwidth=1,orientation="y",stat="identity")+
  geom_fruit(data=thrd,geom=geom_bar,mapping=aes(y=Species,x=OG0003599,fill=Lifestyle,alpha=0.1),pwidth=1,orientation="y",stat="identity")+
  geom_fruit(data=forth,geom=geom_bar,mapping=aes(y=Species,x=OG0006508,fill=Lifestyle,alpha=0.1),pwidth=1,orientation="y",stat="identity")+
  geom_fruit(data=fifth,geom=geom_bar,mapping=aes(y=Species,x=OG0006823,fill=Lifestyle,alpha=0.1),pwidth=1,orientation="y",stat="identity")+
  geom_fruit(data=sixth,geom=geom_bar,mapping=aes(y=Species,x=OG0007180,fill=Lifestyle,alpha=0.1),pwidth=1,orientation="y",stat="identity")+
  geom_fruit(data=seventh,geom=geom_bar,mapping=aes(y=Species,x=OG0007224,fill=Lifestyle,alpha=0.1),pwidth=1,orientation="y",stat="identity")

vert1<-vert+geom_fruit(data=first[which(first$OG0000701>0),],geom=geom_point,mapping=aes(y=Species,x=1,size=OG0000701,color=Lifestyle),alpha=0.5)+
  geom_fruit(data=sec[which(sec$OG0000781>0),],geom=geom_point,mapping=aes(y=Species,x=2,size=OG0000781,color=Lifestyle),alpha=0.5)+
  geom_fruit(data=thrd[which(thrd$OG0003599>0),],geom=geom_point,mapping=aes(y=Species,x=3,size=OG0003599,color=Lifestyle),alpha=0.5)+
  geom_fruit(data=forth[which(forth$OG0006508>0),],geom=geom_point,mapping=aes(y=Species,x=4,size=OG0006508,color=Lifestyle),alpha=0.5)+
  geom_fruit(data=fifth[which(fifth$OG0006823>0),],geom=geom_point,mapping=aes(y=Species,x=5,size=OG0006823,color=Lifestyle),alpha=0.5)+
  geom_fruit(data=sixth[which(sixth$OG0007180>0),],geom=geom_point,mapping=aes(y=Species,x=6,size=OG0007180,color=Lifestyle),alpha=0.5)+
  geom_fruit(data=seventh[which(seventh$OG0007224>0),],geom=geom_point,mapping=aes(y=Species,x=7,size=OG0007224,color=Lifestyle),alpha=0.5)+
  scale_color_brewer(palette="Dark2")+scale_size(range=c(0,15),breaks = c(0.04, 0.16))

vert1


circ<-ggtree(tree,layout="circular")+geom_tiplab(size=1)
circ<-ggtree(tree,layout="circular") #without tip labels
circ$data
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
  geom_fruit(data=seventh,geom=geom_bar,mapping=aes(y=Species,x=OG0007224,fill=Lifestyle),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=sixth,geom=geom_bar,mapping=aes(y=Species,x=OG0007180,fill=Lifestyle),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=fifth,geom=geom_bar,mapping=aes(y=Species,x=OG0006823,fill=Lifestyle),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=forth,geom=geom_bar,mapping=aes(y=Species,x=OG0006508,fill=Lifestyle),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=thrd,geom=geom_bar,mapping=aes(y=Species,x=OG0003599,fill=Lifestyle),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=sec,geom=geom_bar,mapping=aes(y=Species,x=OG0000781,fill=Lifestyle),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  geom_fruit(data=first,geom=geom_bar,mapping=aes(y=Species,x=OG0000701,fill=Lifestyle),alpha=0.8,pwidth=0.3,orientation="y",stat="identity")+
  scale_fill_brewer(palette="Dark2")


circ1

lf<-read.delim("Lignin.tsv")
head(lf)
lf<-lf[c(1,2)]
colnames(lf)<-c("label","Lifestyle1")
head(lf)
#data<-vert$data
data<-circ1$data
head(data)
tail(data)
dat<-left_join(lf,data,by="label")
head(dat)
dat1<-dat[c(1,2,4)]
head(dat1)
tree
#t<-ggtree(tree,layout="roundrect",color="darkgray")
t<-ggtree(tree,layout="circular",color="darkgray")
t
p<-t+geom_hilight(data=dat1,mapping=aes(node=node,fill=Lifestyle1,color=Lifestyle1),extendto=3,alpha=0.5)+
  scale_fill_brewer(palette = "Set2")+
  scale_color_brewer(palette = "Set2")+theme(legend.position = "bottom")
#+geom_tiplab(color="darkgray", size=3)
p
#p+geom_tiplab(size=3,offset=0.3,align=F,color="darkgray",linewidth=0.0001)
p1<-p+new_scale_fill()+geom_fruit(data=first,geom=geom_bar,mapping=aes(y=Species,x=OG0000701,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+
  geom_fruit(data=sec,geom=geom_bar,mapping=aes(y=Species,x=OG0000781,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+
  geom_fruit(data=thrd,geom=geom_bar,mapping=aes(y=Species,x=OG0003599,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+
  geom_fruit(data=forth,geom=geom_bar,mapping=aes(y=Species,x=OG0006508,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+
  geom_fruit(data=fifth,geom=geom_bar,mapping=aes(y=Species,x=OG0006823,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+
  geom_fruit(data=sixth,geom=geom_bar,mapping=aes(y=Species,x=OG0007180,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+
  geom_fruit(data=seventh,geom=geom_bar,mapping=aes(y=Species,x=OG0007224,fill=Lifestyle,alpha=0.1),pwidth=0.2,orientation="y",stat="identity")+
  scale_fill_brewer(palette="Dark2")
p1
p1+theme(legend.title = element_blank())


p1


#circ1<-circ+
geom_fruit(data=seventh,geom=geom_bar,mapping=aes(y=Species,x=OG0007224,fill=Lifestyle),alpha=0.9,pwidth=0.3,orientation="y",stat="identity")+
  scale_fill_manual(values=c("#56b4e9","#a36e46","#a777ba","#f578e4","#459900","#e69f00","#999999","#e86958"))+
  geom_fruit(data=sixth,geom=geom_bar,mapping=aes(y=Species,x=OG0007180,fill=Lifestyle),alpha=0.9,pwidth=0.3,orientation="y",stat="identity")+scale_fill_manual(values=c("#56b4e9","#a36e46","#a777ba","#f578e4","#459900","#e69f00","#999999","#e86958"))+
  geom_fruit(data=fifth,geom=geom_bar,mapping=aes(y=Species,x=OG0006823,fill=Lifestyle),alpha=0.9,pwidth=0.3,orientation="y",stat="identity")+scale_fill_manual(values=c("#56b4e9","#a36e46","#a777ba","#f578e4","#459900","#e69f00","#999999","#e86958"))+
  geom_fruit(data=forth,geom=geom_bar,mapping=aes(y=Species,x=OG0006508,fill=Lifestyle),alpha=0.9,pwidth=0.3,orientation="y",stat="identity")+scale_fill_manual(values=c("#56b4e9","#a36e46","#a777ba","#f578e4","#459900","#e69f00","#999999","#e86958"))+
  geom_fruit(data=thrd,geom=geom_bar,mapping=aes(y=Species,x=OG0003599,fill=Lifestyle),alpha=0.9,pwidth=0.3,orientation="y",stat="identity")+scale_fill_manual(values=c("#56b4e9","#a36e46","#a777ba","#f578e4","#459900","#e69f00","#999999","#e86958"))+
  geom_fruit(data=sec,geom=geom_bar,mapping=aes(y=Species,x=OG0000781,fill=Lifestyle),alpha=0.9,pwidth=0.3,orientation="y",stat="identity")+scale_fill_manual(values=c("#56b4e9","#a36e46","#a777ba","#f578e4","#459900","#e69f00","#999999","#e86958"))+
  geom_fruit(data=first,geom=geom_bar,mapping=aes(y=Species,x=OG0000701,fill=Lifestyle),alpha=0.9,pwidth=0.3,orientation="y",stat="identity")+scale_fill_manual(values=c("#56b4e9","#a36e46","#a777ba","#f578e4","#459900","#e69f00","#999999","#e86958"))
circ1


library(ggpubr)

dev.off()
figure<-ggarrange(cellu_pruned,Hemicellu_pruned,pectin_pruned,lignin_pruned,nrow=1,ncol=4)
figure
p1
#fig2<-ggarrange(plot,circ1,nrow=1,ncol=2)
fig2<-ggarrange(plot,p1,nrow=1,ncol=2)
fig2


fig<-ggarrange(figure,fig2,nrow=2,ncol=1)
fig #size 25*14
##

plotTree(tree)

library(reshape2)
library(data.table)
detach("package:data.table",unload = T)
library(ggtreeExtra)
library(reshape)
library(reshape2)
library(ggforce)
library(tidyr)
plotTree(tree)
gg<-ggtree(tree)+ggtitle("Shared CAZymes from Ascomycota")+geom_tiplab(size=3, align = T)
gg

splist<-data.frame(gg$data)
head(splist)

setnames(splist, "label", "Species") #order of species to be used

p<-read.delim("SR_BaSR_common.tsv")
head(p)
p_sorted<-left_join(splist,p,by="Species")
head(p_sorted)
colnames(p_sorted)
head(p_sorted)
p1<-p_sorted[c(4,6,7,10:17)]
head(p1)
p2<-melt(p1,id=c("Species","x","y","Lifestyle"))
head(p2)
tail(p2)
p3<-na.omit(p2)
head(p3)
tail(p3)

  
p4<-ggplot(p3[which(p3$value>0),],aes(x=1,y=y,fill=Lifestyle))+
  geom_point(aes(x=1,y=y,fill=Lifestyle,size=value, color=Lifestyle), alpha=0.5)+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~variable,scales = "free", nrow =7)+
  scale_color_brewer(palette = "Dark2")+scale_size(range=c(0,10),breaks = c(0.04, 0.16))+
  theme_void()+theme(legend.position = "bottom")+coord_flip()
p4
p3

gg+p4
ggarrange(gg,p4,nrow=2)

head(p_sorted)
p_sort<-melt(p_sorted,id=c("parent","node","branch.length","Species","isTip","x","y","branch","angle","Lifestyle"))
head(p_sort)
p5<-ggtree(tree)%<+%p_sorted[which(p_sorted$OG0000701>0),]
p5<-ggtree(tree)%<+%p_sort[which(p_sort$value>0),]
p5
p6<-p5+geom_point(aes(x=x,y=y,fill=Lifestyle,size=value, color=Lifestyle), alpha=0.5)+
  scale_fill_brewer(palette = "Dark2")+facet_wrap(~variable, scales = "free")+
  scale_color_brewer(palette = "Dark2")+scale_size(range=c(0,7),breaks = c(0.04, 0.16))+theme_void()

p6



pp<-read.delim("fig1_stats.tsv")
head(pp)
colnames(pp)
p<-pp[c(2,3,4,5,6,7,8,10)] #stats
head(p)

p1<-melt(p)
head(p1)
colnames(p1)<-c("Species","Lifestyle","Desc","Count")
head(p1)

p3<-left_join(splist,p1,by="Species")
head(p3)
p3$Species<-as.character(p3$Species)
p3$Species<-factor(p3$Species,levels=unique(p3$Species))
head(p3)

p2<-ggplot(p3,aes(x=Species,y=Count,fill=Lifestyle))+
  geom_bar(stat="identity",alpha=0.4,position="dodge")+
  theme_bw()+scale_fill_brewer(palette="Dark2")+
  ggtitle("Figure 1")+coord_flip()+
  facet_grid(~Desc, scales = "free")+theme(strip.text.x = element_text(angle=0,hjust = 0), strip.background = element_blank())
p2


ipr<-read.delim("3col.tsv",header=F)
head(ipr)
colnames(ipr)<-c("ProtID","IPRID","Desc")
ipr<-unique(ipr)
head(ipr)
ipr$Species<-str_extract(ipr$ProtID,"[^_]+")
head(ipr)
iprcn<-table(ipr$Species,ipr$IPRID)
write.table(iprcn,"ipr_cn_table.tsv",sep="\t",quote=F)


###for removing overlapping/duplicated domains
ipr<-read.delim("001new.tsv",header=F)
head(ipr)
colnames(ipr)<-c("ProtID","start","end","IPRID","Desc")
ipr<-unique(ipr)
head(ipr)
tail(ipr)
View(ipr)
ipr$Species<-str_extract(ipr$ProtID,"[^_]+")
head(ipr)
iprcn<-table(ipr$Species,ipr$IPRID)
write.table(iprcn,"ipr_cn_table_dup.tsv",sep="\t",quote=F)
ipr_desc<-unique(ipr[c(4,5)])
head(ipr_desc)
write.table(ipr_desc,"ipr_anot.tsv",sep="\t",quote=F, row.names = F)

##quinin degradation and pectin degradation

tree<-read.nexus("132sps_final_rooted.nex")
plotTree(tree)
splist<-data.frame(tree$tip.label)
head(splist)
colnames(splist)<-c("Species") #order of species to be used

p<-read.delim("quin_pect.tsv")
head(p)
colnames(p)
p_quin<-p[c(1,2,3,4,9,11)]
head(p_quin)
p_pect<-p[c(1,2,5,6,7,8,10)]
head(p_pect)

p1<-melt(p_quin)
p1<-melt(p_pect)
head(p1)
colnames(p1)<-c("Species","Lifestyle","Desc","Count")
head(p1)


##all 159 IPRs

tree<-read.nexus("132sps_final_rooted.nex")
plotTree(tree)
splist<-data.frame(tree$tip.label)
head(splist)
colnames(splist)<-c("Species") #order of species to be used

pp<-read.delim("ipr_59common.tsv")
head(pp)
colnames(pp)
p<-pp[c(1,2,3:26)] #As
p<-pp[c(1,2,27:40)] #Bs_Cs
p<-pp[c(1,2,41:59)] #Ds_Es_Fs_Gs
p<-pp[c(1,2,60:76)] #HIJKLs
p<-pp[c(1,2,77:89)] #MNs
p<-pp[c(1,2,90:115)] #Ps
p<-pp[c(1,2,116:137)] #RSs
p<-pp[c(1,2,138:161)] #TUVXYZs
head(p)

p1<-melt(p)
head(p1)
colnames(p1)<-c("Species","Lifestyle","Desc","Count")
head(p1)

p3<-left_join(splist,p1,by="Species")
head(p3)
p3$Species<-as.character(p3$Species)
p3$Species<-factor(p3$Species,levels=unique(p3$Species))
head(p3)

p2<-ggplot(p3,aes(x=Species,y=Count,fill=Lifestyle))+
  geom_bar(stat="identity",alpha=0.4,position="dodge")+
  theme_bw()+scale_fill_brewer(palette="Dark2")+
  ggtitle("159_IPR_IDs_common_BaSR and SR")+coord_flip()+
  facet_grid(~Desc, scales = "free")+theme(strip.text.x = element_text(angle=90,hjust = 0), strip.background = element_blank())
p2


####for IPR enrichement
###SRvsLD+WR
list.files()
main<-read.delim("Table_main_IPR.tsv")
head(main)

main$A<-main$SR
main$B<-main$LD+main$WR
head(main)
main$C<-sum(main$SR)-main$A
main$D<-sum(main$SR,main$LD,main$WR)-main$B-main$A-main$C
head(main)
main$E<-sum(main$A,main$B,main$C,main$D)
main1<-main[c(1,10,11,12,13)]
head(main1)
EM=main1 #this is th enrichment matrix
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
write.table(EM, file='SRvsLDWR_IPR1.tsv', quote=FALSE, sep='\t', row.names = F)


###
physac<-read.tree("physac_tree.tree")
plotTree(physac)
physac_tree<-ggtree(physac,color="MediumAquamarine",layout="roundrect",
                    size=0.5,branch.length = "none",linetype="dotted")+
  geom_tiplab(size=4,offset=-1,color="MediumAquamarine")+
  ggtitle("Physalacriaceae")

physac_tree<-ggtree(physac,color="steelblue",layout="roundrect",
                    size=0.5,branch.length = "none",linetype="dotted")+
  geom_tiplab(size=4,offset=-1,color="steelblue")+
  ggtitle("Physalacriaceae")

physac_tree
data<-physac_tree$data
write.table(data,"physac_nodedata.tsv",row.names = F, quote=F, sep="\t")
getwd()
splist<-data.frame(physac$tip.label)
head(splist)
colnames(splist)<-c("Species")

pp<-read.delim("fig1_physac.tsv")
head(pp)
pp$Annotated_Sect<-pp$Secretome-pp$Secretome_no.annot
colnames(pp)
pp$Annotated_SSP<-pp$SSPs-pp$SSPs_no.annot
colnames(pp)

p<-pp
head(p)
gg$data

p1<-melt(p)
head(p1)
colnames(p1)<-c("Names","Species","Lifestyle","Desc","Count")
head(p1)

p3<-left_join(splist,p1,by="Species")
head(p3)
p3$Species<-as.character(p3$Species)
p3$Species<-factor(p3$Species,levels=unique(p3$Species))
head(p3)

busco<-p3[which(p3$Desc=="BUSCO"),]
head(busco)
busco1<-ggplot(busco,aes(x=Lifestyle,y=Species,size=Count),
               fill=Lifestyle)+
  geom_point(aes(x=Lifestyle,y=Species,size=Count),alpha=0.4,color="steelblue")+
  geom_text(aes(label=Count),alpha=0.4,color="steelblue",size=5)+
  theme_void()+
  ggtitle("BUSCO")+theme(legend.position = "none")+
  scale_size(range=c(0,15),breaks = c(84,100))
busco1


prot_size<-p3[which(p3$Desc=="Proteome_size"|p3$Desc=="InterPro"),]
head(prot_size)
prot_size1<-ggplot(prot_size,aes(x=Species,y=Count))+
  geom_bar(stat="identity",alpha=0.4,position="dodge",fill="yellowgreen")+coord_flip()+
  theme_bw()+ggtitle("Proteome sizes")+
  theme(strip.text.x = element_text(angle=90,hjust = 0), strip.background = element_blank())
prot_size1

p
head(p3)
sect<-p3[which(p3$Desc=="Secretome"|p3$Desc=="Annotated_Sect"),]
head(sect)
sect1<-ggplot(sect,aes(x=Species,y=Count))+
  geom_bar(stat="identity",alpha=0.4,position="dodge",fill="orange")+
  theme_bw()+
  ggtitle("Secretome")+coord_flip()+
  theme(strip.text.x = element_text(angle=90,hjust = 0), strip.background = element_blank())
sect1

head(p3)
ssp<-p3[which(p3$Desc=="SSPs"|p3$Desc=="Annotated_SSP"),]
head(ssp)
ssp1<-ggplot(ssp,aes(x=Species,y=Count))+
  geom_bar(stat="identity",alpha=0.4,position="dodge",fill="indianred")+
  theme_bw()+
  ggtitle("SSPs")+coord_flip()+
  theme(strip.text.x = element_text(angle=90,hjust = 0), strip.background = element_blank())
ssp1

head(p)
cazy<-p3[which(p3$Desc=="CAZy"|p3$Desc=="PCWDE"),]
head(cazy)
cazy1<-ggplot(cazy,aes(x=Species,y=Count))+
  geom_bar(stat="identity",alpha=0.4,position="dodge",fill="tan3")+
  theme_bw()+
  ggtitle("CAZy")+coord_flip()+
  theme(strip.text.x = element_text(angle=90,hjust = 0), strip.background = element_blank())
cazy1

one<-physac_tree
one<-p1+theme(legend.position = "none")
two<-busco1+prot_size1+cazy1
three<-sect1+ssp1
one+plot_spacer()+two+three+plot_layout(nrow=1)

two+three
#saved as 25*15 pdf
library(patchwork)
