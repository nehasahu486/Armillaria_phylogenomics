library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(purrr)

setwd("C:/Users/Dell/Documents/TopGO/Armi_compgen_updated/HGT/plant_bact/all_fun_plant_bac/generax_crisscross_Tree/uniref/old/")

list.files()
rm(list=ls())
######################################################################################################################################
######################################################################################################################################
#####_____________________________________________Following steps were done in cluster_____________________________________________###
######################################################################################################################################
######################################################################################################################################

#STEP1: Protein IDs for 675 OGs and their sequences for mmseqs

###files needed
##675_OG.tsv <- list of all 675 OG IDs
##OGID_ProtID_330.tsv <- 3 column file with OG_ID, ProtID and Species abbreviations
##01Physac.fas <-fasta file with sequences of all Physalacriaceae proteins
og<-read.delim("OGID_ProtID_330.tsv")
head(og)
og675<-read.delim("675_OGs.tsv")
head(og675)
og1<-left_join(og675,og,by="mcl")
write.table(og1,"675_ProtIDs.tsv",quote=F,sep="\t",row.names=F)
og2<-og1[2]
head(og2)
write.table(og2,"ProtID_list.tsv",quote=F,sep="\t",row.names=F)

##to extract fasta sequences of proteins from these 675 clusters; type the following in terminal
#seqtk subseq 01Physac.fas ProtID_list.tsv > HT_out.fas

##use the file HT_out.fas as query and 6ada78c1_uniref100tax_v2022_01.fasta as subject and run mmseqs easy-search
## the subject file is uniref sequences from Balazs, with taxid(s) in the headers
#mmseqs easy-search HT_out.fas 6ada78c1_uniref100tax_v2022_01.fasta HT_out.tsv --threads 60 tmp

#STEP2: Complete lineage for taxid
##first check if you have ncbitax2lin in terminal by typing "ncbitax2link -h"
##if nothing pops up, install ncbitax2lin by typing "pip install -U ncbitax2lin"
## if its already installed, then to match the taxid with the complete lineage; do following steps in terminal
##retrive the latest taxonomy dump file from ncbi by the following lines
##I used this, 13_Apr_2022 dated file
#wget -N ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz

#mkdir -p new_taxdump && tar zxf new_taxdump.tar.gz -C ./new_taxdump
##Now the folder new_taxdump, will have the dump files, from which we need only 2 files - nodes.dmp and names.dmp
##Create a dictionary for tax_ids and their full lineage info by typing the following:
#ncbitax2lin --nodes-file new_taxdump/nodes.dmp --names-file new_taxdump/names.dmp
#the output will be saved in a compressed file "ncbi_lineages_[date].csv.gz"
##decompress this file with gzip -d "ncbi_lineages_[date_of_utcnow].csv.gz"
##Finally "ncbi_lineages_[date_of_utcnow].csv" is the file with tax_ids and their complete lineages..this can be used by itself later as well, without doing the previous three steps!
##read this file
ncbi<-read.csv("ncbi_lineages_2022-10-10.csv")
head(ncbi)
ncbi$tax_id<-as.character(ncbi$tax_id) #set tax_ids as characters
sapply(ncbi,class) 
colnames(ncbi)
ncbi1<-ncbi[c(1:9)]
head(ncbi1)
ncbi2<-ncbi1[which(ncbi1$superkingdom==""),]

#ncbi$merge<-pmap(ncbi[2:69],paste,sep="|")

##Format the mmseqs output to keep top 100 hits for each query, arranged by evalues
trns<-read.delim("Transfers.tsv")
head(trns)
colnames(trns)
events<-trns[c(1,5)]
head(events)
trns1<-separate_rows(events,Reciever_clade,sep=", ",convert = T)
head(trns1)
colnames(trns1)[2]<-"ProtIDs"

big<-read.delim("HT_out.tsv",header=F)
head(big)
colnames(big)<-c("query","subject","perc_id","align_len","mismatches","gap_opening","query_start","query_end","sub_start","sub_end","evalue","bitscore")
head(big)

max(big$perc_id)
min(big$perc_id)

#big<-big[which(big$perc_id>=0.4),]
big1_1<-big[,c(1,2,3,4,11)]
head(big1_1)

big1<-filter(big1_1,query %in% trns1$ProtIDs)
length(unique(big1$query))


big2<-big1 %>% 
  group_by(query) %>%
  arrange(evalue) %>%
  slice_head(n=100) %>%
  ungroup()
head(big2)

table(big2$query)
test<-big2[which(big2$query=="Armbor1_1005786"),]
head(test)
min(test$evalue)
tail(test)
max(test$evalue)


head(big2)
big2[c('subject','tax_id','Kingdom')]<-str_split_fixed(big2$subject,':',3)
head(big2)
big2$subject<-gsub("UniRef100_","",big2$subject)
big2$tax_id<-gsub("t","",big2$tax_id)
head(big2)
unique(big2$tax_id)
write.table(big2,"filtered_uniref100_perc_ids.tsv",quote=F,sep="\t",row.names = F)
#write.table(big2,"filtered_uniref500.tsv",quote=F,sep="\t",row.names = F)

#write.table(big2,"filtered_uniref250.tsv",quote=F,sep="\t",row.names = F)

#write.table(big2,"filtered_uniref.tsv",quote=F,sep="\t",row.names = F)

th<-big2
head(th)
length(unique(th$query))


th1<-left_join(th,ncbi,by="tax_id")
head(th1)
unique(th1$kingdom)
colnames(th1)
th2<-th1[c(1:14,44,59,60,64:66)]
head(th2)
colnames(th2)

write.table(th2,"full_lineage_hits100_percid.tsv",quote=F,sep="\t",row.names = F)

#write.table(th1,"full_lineage_ALL_hits.tsv",quote=F,sep="\t",row.names = F)
#write.table(th2,"full_lineage_hits250.tsv",quote=F,sep="\t",row.names = F)
#write.table(th2,"full_lineage_hits500.tsv",quote=F,sep="\t",row.names = F)

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
  arrange(-Freq) %>%
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


unique(th6$Var2)
head(th6)
max(th6$Freq)
sapply(th6,class)
th6$Freq<-as.character(th6$Freq)
tabl<-dcast(th6,Var1~Var2)
head(tabl)
tabl[is.na(tabl)] <- 0
head(tabl)
tabl<-tabl[c(1,5,3,2,4)]
head(tabl)
colnames(tabl)[1]<-"ProtIDs"
head(evnt)
evnt2<-left_join(tabl,evnt,by="ProtIDs")
head(evnt2)
evnt3<-evnt2[c(6,2,3,4,5)]
head(evnt3)
evnt3<-na.omit(evnt3)
evnt3<-arrange(evnt3,event)
head(evnt3)
sapply(evnt3,class)
cols.num <- c("Physalacriaceae","Basidiomycota","Ascomycota","Other_taxa")
evnt3[cols.num] <- sapply(evnt3[cols.num],as.integer)
evnt4<-aggregate(.~event,evnt3,mean)
head(evnt4)
head(evnt3)
evnt4$Total<-evnt4$Physalacriaceae+evnt4$Basidiomycota+evnt4$Ascomycota+evnt4$Other_taxa
head(evnt4)
max(evnt4$Total)
evnt4$Physac<-(evnt4$Physalacriaceae/evnt4$Total)*100
evnt4$Basid<-(evnt4$Basidiomycota/evnt4$Total)*100
evnt4$Asco<-(evnt4$Ascomycota/evnt4$Total)*100
evnt4$OT<-(evnt4$Other_taxa/evnt4$Total)*100
head(evnt4)
write.table(evnt4,"Tophits%.tsv",quote=F,sep="\t",row.names = F)
##th7
th7<-read.delim("lineage_MAXcount.tsv")
head(th7)
unique(th7$Category)
th8<-data.frame(table(th7$Category))
th9<-th8[which(th8$Freq<=10000 & th8$Freq >30),]
head(th9)
length(unique(th9$Var1))
ggplot(th9)+geom_bar(aes(x=fct_reorder(Var1,Freq,.desc=T),y=Freq,fill=Var1),stat="identity")+
  scale_fill_brewer(palette = "Set1")+theme_bw()+theme(axis.text = element_text(color="black",size=14,angle=90),legend.position = "none")+
  xlab("")+ylab("Number of Top hits in each category")


##Transfer with 70% bootstrap support
##reciever clade has 70% Basido
##donor clade has 70% Asco
##sister clade of transfer is 70% Asco

trns<-read.delim("Transfers.tsv")
head(trns)
trns$event<-paste0("event",seq.int(nrow(trns)))
head(trns)
colnames(trns)
events<-trns[c(6,4)]
head(events)

trns1<-separate_rows(trns,Reciever_clade,sep=", ",convert = T)
head(trns1)
trns1$Species<-sub('_[^_]*$','',trns1$Reciever_clade)
head(trns1)
colnames(trns1)[4]<-"ProtIDs"

trans2<-left_join(trns1,events,by="event")
head(trans2)
phys<-c("Oudmuc1","Hymrad1","Cylto1","FUN","Flave","Guyne1","Armtab1","Armect1","Armnab1","Armnov1","Armmel1","Armme1_1","Armosto1","Armost1","Armcep1","Armfus","Armfum1","Armga1","Armbor1","Armlut1")
evnt<-trans2[c(4,6)]
head(evnt)

trans3<-filter(trans2,trans2$Species %in% phys)
unique(trans3$Species)
write.table(trans3,"Check_transfers.tsv",sep="\t",row.names = F,quote=F)
