library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(purrr)


######################################################################################################################################
######################################################################################################################################
####_____________________________________________Max taxon occupancy for each HT event_____________________________________________###
######################################################################################################################################
######################################################################################################################################

th2<-read.delim("full_lineage_hits100_percid.tsv")
head(th2)
unique(th2$order)
sapply(th2,class)
th2$Category<-if_else(th2$family %in% c("Physalacriaceae"),"Physalacriaceae",th2$phylum)
selec<-c("Physalacriaceae","Basidiomycota","Ascomycota")
head(th2)
#th2$Category<-if_else(th2$Category %in% selec,th2$Category,"Other_taxa")
#head(th2)

th4<-data.frame(table(th2$query,th2$Category))
head(th4)
#write.table(th4,"lineage_count.tsv",quote=F,sep="\t",row.names = F)

unique(th4$Var2)
sapply(th4,class)
th4$Var2<-as.character(th4$Var2)
unique(th4$Var2)
head(th4)
th5<-th4[which(th4$Freq>0),]
head(th5)
th6<-arrange(th5,Var1,-Freq)
head(th6)
#write.table(th6,"lineage_count_table.tsv",quote=F,sep="\t")


th7<-th6 %>% 
  group_by(Var1) %>%
  arrange(-Freq) %>% #arranges based on decreasing Frequency of hits in each taxonomic category
  slice_head(n=1) %>%
  ungroup()
head(th7)
length(unique(th7$ProtID))
colnames(th7)<-c("ProtIDs","Category","Top_hits")
write.table(th7,"lineage_MAXcount.tsv",quote=F,sep="\t")


trans<-read.delim("Check_transfers.tsv")
head(trans)
colnames(trans)
trans1<-left_join(trans,th7,by="ProtIDs")
head(trans1)
colnames(trans1)
write.table(trans1,"Trans_tophits.tsv",quote=F,row.names = F,sep="\t")

