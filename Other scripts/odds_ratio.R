library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ape)
library(reshape2)
library(phytools)
library(ggtree)
library(reshape2)
library(plyr)
library(forcats)
library(ggpubr)
library(ggrepel)
library(gplots)

#all DEG files
counts<-read.delim("DEG_counts.tsv")
head(counts)
unique(counts$Setup)
#counts<-counts[which(counts$Setup!="Developmental"),]
#head(counts)
#unique(counts$Setup)
unique(counts$Species)
#counts$Downreg<-counts$Downreg*(-1)
#counts1<-melt(counts[c(1,5,4,2,3)])
#head(counts1)
#tail(counts1)
#ggplot(counts1)+geom_bar(aes(y=Experiment,x=value,fill=variable),stat="identity",alpha=0.8)+
#  facet_grid(Setup~Species,scales="free")+theme_classic()+scale_fill_manual(values=c("seagreen","mediumorchid1"))

deg<-read.delim("All_DEGs.tsv")
head(deg)
unique(deg$Setup)
#deg<-deg[which(deg$Setup!="Developmental"),]
#head(deg)
#unique(deg$Setup)
colnames(deg)<-c("ProtID","logFC","Experiment","Setup","Status","Sample")
head(deg)

func<-read.delim("Functions.tsv")
head(func)

unique(func$Function)
func1<-unique(func %>% 
               mutate(Function = strsplit(as.character(Function), ",")) %>% 
               unnest(Function))
head(func1)
write.table(func1,"Func_segg.tsv",quote=F,sep="\t",row.names = F)

func_all<-unique(aggregate(.~ProtID,func1,toString))
head(func_all)
write.table(func_all,"Func_all.tsv",quote=F,sep="\t",row.names = F)

##subset
##read_selected function files
func2<-read.delim("Func_counts1.tsv") #made after manual checking
head(func2)
tail(func2)
unique(func2$Plot_sep)
fun<-unique(func2[c(1,4,5)])
head(fun)
fun1<-unique(aggregate(.~ProtID,fun,toString))
head(fun1)
write.table(fun1,"ProtIDs_agg_func.tsv",quote=F,row.names = F,sep="\t")

func3<-func2[!is.na(func2$Plot),]
head(func3)
tail(func3)
unique(func3$Plot_sep)
func4<-left_join(func3,deg,by="ProtID")
func4$Species<-sub('_[^_]*$','',func4$ProtID)
unique(func4$Species)
head(func4)
func4_1<-unique(func4[c(1,12,5)])
head(func4_1)
func5<-data.frame(table(func4_1$Species,func4_1$Plot_sep))
head(func5)
head(func4)
func4_2<-unique(func4[c(1,12,5,8)])
head(func4_2)
func6<-data.frame(table(func4_2$Experiment,func4_2$Plot_sep))
head(func6)
unique(func6$Var2)
func7<-left_join(func5,func6,by="Var2")
head(func7)
colnames(func7)<-c("Spec","Function","Total","Experiment","DEGs")
head(func7)

head(counts)
func8<-left_join(func7,counts,by="Experiment")
head(func8)
func9<-func8[(func8$Spec==func8$Species),]
head(func9)

write.table(func9,"All_func_degs1.tsv",quote=F,sep="\t",row.names = F)

head(func9)
main<-func9[c(2,5,3,6,9,8,4,7)]
head(main)
unique(main$Function)
write.table(main,"Func_big1.tsv",quote=F,sep="\t",row.names = F)



##########_____enrichment____################
main<-read.delim("Func_big1.tsv")
head(main)
main$A<-main$DEGs
main$B<-main$Total#-main$DEGs
main$C<-main$Tot_DEG-main$DEGs
main$D<-main$Proteome-main$A-main$B-main$C
head(main)

Fstat<-sapply(1:nrow(main[c(1,9:12)]), function(i2) fisher.test(matrix(unlist(main[i2,9:12]),nrow = 2)))
main$Pvalue<-unlist(t(Fstat)[,1])
main$Odds_ratio<-unlist(t(Fstat)[,3])
main$FDR <- p.adjust(main$Pvalue, method= "BH")
PT<-as_tibble(main$Pvalue)
FDR<-as_tibble(main$FDR)
head(main)
szig<-FDR
szig[szig[1]<0.0001,]<-4
szig[szig[1]<0.001,]<-3
szig[szig[1]<0.05,]<-2
szig[szig[1]<=1,]<-0
szigm<-szig
szigm[szigm[1]>1,]<-2
szigm[szigm[1]<1,]<-1
main$sig<-szig$value
main$sigm<-szigm$value
main$Dir<-main$A/main$B > as.numeric(main$C)/as.numeric(main$D)
head(main)
main$Status<-ifelse(grepl(pattern = "Upreg", main$Experiment), "Upregs", "Downregs")
write.table(main, file='Enriched_degs1.tsv', quote=FALSE, sep='\t', row.names = F)


###plot p-values
main<-read.delim("Enriched_degs1.tsv") #select comparisons from Enriched_degs file
unique(main$Function)

lev<-c("Cellulose","Hemicellulose","Pectin","Lignin", "Suberin",
       "Chitin","Glucan","CPPs","CAP domain/superfam","CDA","Catalase","LysM","SOD",
       "Glutathione","Lactone","Ergothione","Peroxiredoxin","Redoxin", "Thioredoxin","HT","Anot_SSP","Unanot_SSPs")
head(main)
ggplot(main[which(main$Status=="Upregs"),],aes(y=Total,x=factor(Function,lev),fill=Setup))+
  geom_bar(stat="identity",position="dodge")+theme_classic()+scale_fill_brewer(palette = "Spectral")+
  theme(axis.text.x=element_text(angle=90))+ylab("")+xlab("")+
  ggtitle("Functions enriched in Upregulated genes")

upreg<-ggplot(main[which(main$Status=="Upregs"),],aes(x=Experiment,y=factor(Function,lev),alpha=Dir,color=FDR))+
  geom_point(size=3)+theme_classic()+scale_color_gradient2(low="red",mid="red",high="blue")+
  facet_wrap(~Setup,scales="free",nrow=1)+theme(axis.text.x=element_text(angle=90))+ylab("")+xlab("")+
  ggtitle("Functions enriched in Upregulated genes")
upreg

down<-ggplot(main[which(main$Status=="Downregs"),],aes(x=Experiment,y=factor(Function,lev),alpha=Dir,color=FDR))+
  geom_point(size=3)+theme_classic()+scale_color_gradient2(low="red",mid="red",high="blue")+
  facet_wrap(~Setup,scales="free",nrow=1)+theme(axis.text.x=element_text(angle=90))+ylab("")+xlab("")+
  ggtitle("Functions enriched in Downregulated genes")

ggarrange(upreg,down,nrow=2)


####Plot odds ratio
plot<-read.delim("Enriched_degs1.tsv") #select comparisons from Enriched_degs file
head(plot)
unique(plot$Status)

x<-dcast(plot,Experiment~Function,value.var = "Odds_ratio")
head(x)
write.table(x,"plot_hm.tsv",quote=F,sep="\t",row.names=F)


x1<-read.delim("plot_hm1.tsv",row.names = 1)  #edited before to order columns
head(x1)
rownames(x1)
mycols = colorRampPalette(c("blue","white","red"))(100) #this is for the color.. blue being lowest and red being the highest

x2<-as.matrix(x1[c(1:10,13),]) #upregs 3exp
x3<-as.matrix(x1[c(27:36,39),]) #downregs 3exp
z <- heatmap.2(x2, dendrogram = "none", trace = "none",  
               col = mycols, keysize = 1, key.title = "", Colv = FALSE,Rowv = F,
               scale = c("col"), cexRow = 0.75, cexCol = 1,  
               adjRow = c(0,NA), sepcolor="white" ,margins = c(10,8), 
               lwid = c(2,12),
               lhei= c(5,15), main ="enriched_upregs",
               rowsep = c(5,7,9,10),sepwidth = c(0.1,0.1)) 

z1 <- heatmap.2(x3, dendrogram = "none", trace = "none",  
                col = mycols, keysize = 1, key.title = "", Colv = FALSE,Rowv = F,
                scale = c("col"), cexRow = 0.75, cexCol = 1,  
                adjRow = c(0,NA), sepcolor="white" ,margins = c(10,8), 
                lwid = c(2,12),
                lhei= c(5,15), main ="enriched_downregs",
                rowsep = c(5,7,9,10),sepwidth = c(0.1,0.1)) 

rownames(x1)

x2<-as.matrix(x1[c(16:26),]) #upregs 3exp
x3<-as.matrix(x1[c(42:52),]) #downregs
head(x2)
head(x3)

dev.off()

