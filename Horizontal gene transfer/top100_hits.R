library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(purrr)

setwd("C:/Users/Dell/Documents/TopGO/Armi_compgen_updated/HGT/plant_bact/all_fun_plant_bac/generax_crisscross_Tree/uniref/old/")

list.files()
rm(list=ls())



######################################################################################################################################
#STEP1: Protein IDs for 675 OGs and their sequences for mmseqs
######################################################################################################################################
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

######################################################################################################################################
#STEP2: Complete lineage for taxid
######################################################################################################################################
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

######################################################################################################################################
#STEP3: Parsing mmseqs output file and retaining top 100 hits for each query protein
######################################################################################################################################
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
  slice_head(n=100) %>% #this selects the top 100 hits based on evalues.. change this number if more/less hits are needed
  ungroup()
head(big2)

table(big2$query)
test<-big2[which(big2$query=="Armbor1_1005786"),] #cross check for a few query proteins
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
