library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ape)
library(reshape2)
library(phytools)
library(ggtree)
library(reshape2)
library(forcats)
library(ggpubr)
library(ggforce)
library(data.table)


##############################################################################################################################
##############################################################################################################################

##########################_________For parsing mmseqs output and alien-index calculation_________#############################

##############################################################################################################################
##############################################################################################################################



#######################################################

############_________For e-value_________##############

#######################################################


##do the following steps on HPC in case of larger file sizes
#######Read the evalue file##########
setwd("/work/nsahu/132sps_Armi/HGT/all_Asco/")
rm(list=ls())
list.files()
library(dplyr)
library(data.table)
x<-fread("Physac_big_search",header=F)
|--------------------------------------------------|
  |==================================================|
colnames(x)<-c("query","subject","perc_id","align_len","mismatches","gap_opening","query_start","query_end","sub_start","sub_end","evalue","bitscore")
qlen<-len #len is file with protID and seq length
colnames(qlen)<-c("query","qlen")
x1<-left_join(x,qlen,by="query")
slen<-len
colnames(slen)<-c("subject","slen")
x2<-left_join(x1,slen,by="subject")
x3<-x2[which(x2$query!=x2$subject),]
head(x3)
x3$query_cov<-((x3$query_end-x3$query_start)/x3$qlen)*100
x3$sub_cov<-((x3$sub_end-x3$sub_start)/x3$slen)*100
write.table(x3,"Phyac_big_search_query_sub_cov.tsv",sep="\t",row.names=F,quote=F)
write.table(x3,"Phyac_big_search_query_sub_cov.tsv",sep="\t",row.names=F,quote=F)
x4<-x3[which(x3$query_cov>=30 & x3$sub_cov>=30 & x3$evalue<1e-5),]
write.table(x4,"Phyac_big_search_cov30_eval1e_5.tsv",sep="\t",row.names=F,quote=F)
min(x4$sub_cov)
30
min(x4$query_cov)
30
max(x4$evalue)
9.994e-06


####################
######################

big_file<-fread("Phyac_big_search_cov30_eval1e_5.tsv")
big<-big_file[,c(1,2,11)]
head(big)
colnames(big)<-c("query","subject","evalue")
big$query_sps<-sub('_[^_]*$','',big$query)
big$Species<-sub('_[^_]*$','',big$subject)
unique(big$query_sps)
unique(big$Species)

###remove hits if query and subject are from the same species
big1<-big[(big$query_sps!=big$Species),]
head(big1)
####read the lifestyle file for orders and categories
lf<-read.delim("Lifestyle.txt")
head(lf)
big1_1<-left_join(big1,lf,by="Species")
head(big1_1)
unique(big1_1$Class)
unique(big1_1$Category)
write.table(big1_1,"evalue_hits.tsv",quote=F,sep="\t",row.names=F)
head(big1_1)


####START FROM HERE (can be done without HPC)
list.files()
big1<-fread("evalue_hits_FILT.tsv")
head(big1)
big2<-big1 %>% 
  group_by(query,Category) %>%
  arrange(evalue) %>%
  slice(1) %>%
  ungroup()
head(big2)
write.table(big2,"evalue_hits1.tsv",quote=F,sep="\t",row.names = F)


#################################################################################################

##############________FOR ALIEN INDEX CALCULATION START FROM HERE________##############

#################################################################################################


####START FROM HERE
rm(list=ls())
filt<-read.delim("evalue_hits1.tsv")
head(filt)
filt1<-filt[c(1,7,3)]
head(filt1)
filt2<-unique(filt1)
head(filt2)
tail(filt2)
min(filt2$evalue)
detach(package:data.table)
filt3<-dcast(filt2,query~Category)
head(filt3)
##replace NA by 1
filt3[is.na(filt3)] <- 1
head(filt3)
write.table(filt3,"check.tsv",quote=F,sep="\t",row.names=F)
###AI calculation
filt3$AI<-log(filt3$Basid+1e-200)-log(filt3$Donor+1e-200)
filt3$AI_log10<-log10(filt3$Basid+1e-200)-log10(filt3$Donor+1e-200)
head(filt3)
min(filt3$AI_log10)
max(filt3$AI_log10)

#add anot and OG_IDs from 67 sps dataset
anot<-read.delim("67sps_OGID_ProtID_IPRs.tsv")
head(anot)
unique(anot$Species)
colnames(anot)<-c("OG_ID67","query","Recipient","IPR","Desc")
filt4<-left_join(filt3,anot,by="query")
head(filt4)
filt4$Type<-ifelse(filt4$AI_log10>1,"HT","VT")
head(filt4)
table(filt4$Type)
write.table(filt4,"AI_evalue_calc.tsv",quote=F,sep="\t",row.names = F)

