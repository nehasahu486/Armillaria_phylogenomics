
rm(list=ls())
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ape)
library(phytools)
library(ggtree)
library(reshape2)
library(forcats)
library(ggpubr)
library(dplyr)
library(stringr)
library(UpSetR)
library(phytools)
library(ggtree)
library(ape)

#######################################################################################
#######################################################################################
############_____________Species specific and core genes in Armillaria________#########
#######################################################################################
#######################################################################################
order<-read.delim("Order_of_sps.tsv") #list of species to arrange their order in plots
head(order)  
ord<-(order$Species)


up<-read.delim("Armi_for_upset.tsv") #orthofinder counts for Physalacriaceae species
head(up)
#select Orthogroups with atleast 70% of the Armillaria species present
up$Total<-rowSums(up[c(2:16)])
head(up)
up1<-up[c(1:16)][which(up$Guyne1==1 & up$Total>=12),]
head(up1)
length(rownames(up1))
a<-colnames(up1[c(2:16)])
b<-colnames(up1[c(2:15)])
c<-colnames(up1[c(2,3,5:16)])

upset(up1,nsets=16,main.bar.color="black",mainbar.y.label = "Number of Orthogroups",
      shade.color = "grey",
      mb.ratio=c(0.55,0.45),keep.order = T,order.by = "degree",decreasing=T,
      query.legend = "bottom", nintersects =NA,sets=ord,text.scale = c(2, 2, 2, 2, 2, 2))
#saved at 10*6

upset(up1,nsets=16,main.bar.color="black",mainbar.y.label = "Number of Orthogroups",
      shade.color = "grey",
      mb.ratio=c(0.55,0.45),keep.order = T,order.by = "freq",decreasing=T,
      query.legend = "bottom", nintersects =5,sets=ord,
      queries = list(list(query = intersects,params =a, color = "steelblue",active=T),
                     list(query = intersects,params =b, color = "mediumseagreen",active=T),
                     list(query = intersects,params =c, color = "orchid",active=T)))
#saved at 10*6


##annotations
anot<-read.delim("67sps_OGID_ProtID_IPRs.tsv")
head(anot)
anot1<-filter(anot,Ortho67 %in% up1$Var1)
head(anot1)
unique(anot1$Species)
write.table(anot1,"Check_DEGs.tsv",sep="\t",row.names = F,quote=F)

exp<-read.delim("All_DEGs.tsv")
head(exp)
exp1<-left_join(anot1,exp,by="ProtID")
head(exp1)
exp2<-exp1[which(exp1$Setup!="NA"),]
head(exp2)
write.table(exp2,"Exp_DEGs.tsv",sep="\t",row.names = F,quote=F)
ggplot(exp2)+geom_bar(aes(y=Experiment,fill=Status))+facet_wrap(~Setup,scales="free")+
   scale_fill_manual(values=c("steelblue1","lightgreen"))+theme_classic()+ggtitle("Expression of core-Armillaria genes")


head(up1)
anot<-read.delim("OG_IDs_67sps and % of IPR_IDs.tsv")
head(anot)
colnames(anot)[1]<-c("Var1")

anot1<-left_join(up1,anot,by="Var1")
head(anot1)
unique(anot1$Species)
write.table(anot1,"core_OG_IDs.tsv",sep="\t",row.names = F,quote=F)

core<-read.delim("core_OG_IDs.tsv")
head(core)
colnames(core)
core1<-unique(core[c(1,18,19)])
head(core1)
core2<-aggregate(.~Ortho67,core1,paste,collapse="|")
head(core2)
colnames(core2)[1]<-c("Var1")
head(up1)
length(rownames(up1))
core3<-left_join(up1,core2,by="Var1")
head(core3)
write.table(core3,"core_OG_IDsanot.tsv",sep="\t",row.names = F,quote=F)
