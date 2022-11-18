#################################################################################################################################
#############################                                                                       #############################
#############################   Visualization of COMPARE output files for Dataset 1 (66 species)    #############################
#############################                                                                       #############################
#################################################################################################################################

library(phytools)
library(dplyr)
library(ape)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(reshape2)
library(ggforce)

###STEP1 : 
#Read the 66 species tree first. 
tree<-read.tree("species.tree.nwk")
plotTree(tree)
nodelabels() #check node labels
##Following is the LadderizeTree function from COMPARE scripts
##To ensure matching of the node numbers from COMPARE output and species tree, Use this function to ladderize the species tree 
LadderizeTree <- function(tree, temp_file = "temp", orientation = "left"){
  if(file.exists(paste0("./", temp_file))){
    stop("The chosen temporary file exists! Please choose an other temp_file name")
  }
  if(orientation == "left"){
    right <- FALSE
  }else{
    right <- TRUE
  }
  tree_temp <- ladderize(tree, right = right)
  write.tree(tree_temp, file = paste0("./", temp_file, ".tre"))
  tree_lad <- read.tree(paste0("./", temp_file, ".tre"))
  file.remove(paste0("./", temp_file, ".tre"))
  return(tree_lad)
}

##this is the Ladderzied tree with node numbers matching with COMPARE output
tree<-LadderizeTree(tree)
plotTree(tree)

##Check how the tree looks in ggtree
tree_dim<-ggtree(tree)
tree_dim+geom_tiplab(offset = 0.01,size=7)+geom_text(aes(label=node),color="red",hjust=-0.2, size=7) #Check how the tree looks
##save the above tree -> pdf 20*30 portrait



###STEP2 : Exract tree data from the species tree containing the node number and their location in x and y axis
tree_axis<-tree_dim$data 
tree_axis



###STEP3 : Read the compare output file - we need only allDLT_data_allCLus.tsv for SUMMARY plot
##the allDLT_data_allCLus.tsv is a summary file for sum of net gains, duplications, losses and copy numbers at each node of the species tree
data<-read.delim("allDLT_data_allClus.tsv")
head(data)



###STEP4 : Plotting attributes like Net Gains, Copy numbers, Duplications and Losses
##NET GAINS
ng<-data[c(1,4)] #we need only 1st and 4th column
head(ng)# check if the columns are correct, i.e. treenode and sumNet_gains
colnames(ng)
colnames(ng)<-c("node","Net_Gains_color")
ng$Net_Gains <-abs(ng$Net_Gains_color) #changes all negative values to positive
head(ng)

#join the ng data with tree_axis data
data_ng<-left_join(ng,tree_axis,by="node")
head(data_ng)
a<-min(data_ng$Net_Gains)
b<-median(data_ng$Net_Gains)
c<-max(data_ng$Net_Gains)
a
b
c
head(data_ng)
ng_tree<-ggtree(tree, color="grey69",layout = "roundrect")%<+%data_ng
ng_tree

ngtree<-ng_tree+geom_point(aes(x=x,y=y,size=Net_Gains,color=Net_Gains_color>0,fill=Net_Gains_color>0,stroke=0),alpha=0.5)+
  geom_text(aes(label=Net_Gains),size=4,color="black")+
  scale_size(range=c(0,25), breaks = c(a,b,c))+
  scale_color_manual(values=c("orange","cadetblue3"), name=NULL,labels=c("Contraction", "Expansion"))+
  ggtitle("Summary Net Gains")+theme(legend.title = element_blank())+
  theme(legend.position = "none")
ngtree

##For Physalacriaceae clade
ngtree_physac<-viewClade(ngtree,MRCA(ngtree,"Armga1","Hymrad1"))
ngtree_physac
##save pdf at 10*15 portrait

###COPY NUMBERS
##PLOT copy numbers _ for main fig
head(data)
cn<-data[c(1,5)] #we need only 1st and 5th column
head(cn)# check if the columns are correct, i.e. treenode and sumNet_gains
colnames(cn)
colnames(cn)<-c("node","Copynum")
head(cn)
max(cn$Copynum)
min(cn$Copynum)
#join the cn data with tree_axis data
data_cn<-left_join(cn,tree_axis,by="node")
head(data_cn)
max(data_cn$Copynum)
min(data_cn$Copynum)

a<-min(data_cn$Copynum)
b<-median(data_cn$Copynum)
c<-max(data_cn$Copynum)
a
b
c
head(data_cn)
cn_tree<-ggtree(tree, color="lightsteelblue2",layout = "roundrect")%<+%data_cn
cn_tree

cntree<-cn_tree+geom_point(aes(x=x,y=y,size=Copynum,stroke=0),color="mediumorchid",alpha=0.5)+
  geom_text(aes(label=Copynum),size=3, color="black", alpha=0.85)+
  scale_size(range=c(0,15), breaks = c(a,b,c))+
  geom_tiplab(offset=0.02,size=3)+
  ggtitle("Ancestral sizes")+theme(legend.title = element_blank())+
  theme(legend.position = "bottom")
cntree

##For Physalacriaceae clade
cntree_physac<-viewClade(cntree,MRCA(ngtree,"Armga1","Hymrad1"))
cntree_physac
##save pdf at 10*15 portrait
##similarly data from other columns such as "Duplications" and "Losses" can be plotted onto the tree



#STEP5 : For plotting data for any specific cluster(s), use the "catDLT.tsv" file
data<-read.delim("catDLT.tsv")
head(data)
#select the interesting clusters (eg. intradiol dioxygenases)
data1<-data[(data$Cluster=="Cluster0000303"|data$Cluster=="Cluster0006804"|data$Cluster=="Cluster0007729"),]

head(data1)
tail(data1)
sapply(data1,class)
data1[,1]
data1[c(2:7)]
data1[,1]<-sapply(data1[,1],as.character)
data1[c(2:7)]<-sapply(data1[c(2:7)],as.numeric)
sapply(data1,class)
head(data1)
summary(data1)
data1<-data1[c(1:7)]
head(data1)
sapply(data1,class)
data2<-melt(data1)
head(data2)
data3<-aggregate(.~treenode, data1, sum)
head(data3)
data3[,1]<-sapply(data3[,1],as.numeric)
head(data3)
data4<-data3[order(data3$treenode),]
head(data4) #this is the final sum of interesting cluster for visualization

##PLOT Copynums
head(data4)
cn<-data4[c(1,7)]
head(cn)# check if the columns are correct, i.e. treenode and sumNet_gains
colnames(cn)
colnames(cn)<-c("node","copyNum")
head(cn)

#join the cn data with tree_axis data
data_cn<-left_join(cn,tree_axis,by="node")
head(data_cn)
a<-min(data_cn$copyNum)
b<-median(data_cn$copyNum)
c<-max(data_cn$copyNum)
a
b
c
sapply(data_cn,class)
head(data_cn)
cn_tree<-ggtree(tree,color="steelblue",layout="roundrect")%<+%data_cn
cn_tree+geom_tiplab(offset=0.05,size=3)+geom_text(aes(label=ifelse(copyNum>0,as.numeric(copyNum),"")),size=3, color="black", alpha=0.85)+
  geom_point(aes(x=x,y=y,size=copyNum),color="steelblue",alpha=0.3)+
  scale_size(range=c(0,15), breaks = c(1,15,30))+ggtitle("Intradiol ring cleavage dioxygenases")
##save pdf at 10*15 portrait
