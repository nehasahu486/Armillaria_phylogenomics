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


rm(list=ls())
setwd("C:/Users/Dell/Documents/TopGO/Armi_compgen_updated/HGT/plant_bact/all_fun_plant_bac/cov_50_eval1E_5/")
list.files()



#######################################################

############_________For e-value_________##############

#######################################################


##do the following steps on cluster
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


####START FROM HERE _in your computer
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

##best hits for basidio and asco
head(filt4)
filt6<-filt4[c(1,4,5,6:10)]
head(filt6)
head(filt)
filt6_1<-left_join(filt6,filt,by="query")
head(filt6_1)
unique(filt6_1$Class)

asco<-filt6_1[which(filt6_1$Class=="Ascomycota"),]
head(asco)
write.table(asco[which(asco$Type=="HT"),],"Asco_best.tsv",quote=F,sep="\t",row.names = F)

zoo<-filt6_1[which(filt6_1$Class=="Zoopagomycotina"),]
head(zoo)
write.table(zoo[which(zoo$Type=="HT"),],"Zoopago_best.tsv",quote=F,sep="\t",row.names = F)

Mucor<-filt6_1[which(filt6_1$Class=="Mucoromycota"),]
head(Mucor)
write.table(Mucor[which(Mucor$Type=="HT"),],"Mucor_best.tsv",quote=F,sep="\t",row.names = F)

Chyt<-filt6_1[which(filt6_1$Class=="Chytidiomycota"),]
head(Chyt)
write.table(Chyt[which(Chyt$Type=="HT"),],"Chyt_best.tsv",quote=F,sep="\t",row.names = F)

EDF<-filt6_1[which(filt6_1$Class=="EDF"),]
head(EDF)
write.table(EDF[which(EDF$Type=="HT"),],"EDF_best.tsv",quote=F,sep="\t",row.names = F)

Plant<-filt6_1[which(filt6_1$Class=="Plant"),]
head(Plant)
write.table(Plant[which(Plant$Type=="HT"),],"Plant_best.tsv",quote=F,sep="\t",row.names = F)

Bacteria<-filt6_1[which(filt6_1$Class=="Bacteria"),]
head(Bacteria)
write.table(Bacteria[which(Bacteria$Type=="HT"),],"Bacteria_best.tsv",quote=F,sep="\t",row.names = F)

basid<-filt6_1[which(filt6_1$Class=="Basidiomycota"),]
head(basid)
write.table(basid[which(basid$Type=="VT"),],"Basido_best.tsv",quote=F,sep="\t",row.names = F)


#check the number of Horizontal transfers (HT) and Vertical transfers (VT) in each species
#read the merged best hits file
filt4<-read.delim("Merged_best_hits.tsv")
head(filt4)
a<-ggplot(filt4)+geom_bar(aes(y=Recipient,fill=Type))+
  facet_wrap(Division~Type, scales="free")+
  scale_fill_manual(values=c("steelblue","tan"))+theme_classic()+xlab("# of genes")
a




lf<-read.delim("Lifestyle.tsv")
head(lf)
filt5<-unique(left_join(filt4,lf[c(1,4)],by="Species"))
head(filt5)
filt7<-filt5[which(filt5$Type=="HT"),]
filt7<-filter(filt5,Type=="HT",Recipient!="Agabi_varbisH97_2",Recipient!="Copci_AmutBmut1",Recipient!="Schco3")
head(filt7)
write.table(filt7,"Merged_bestHT_lf.tsv",quote=F,sep="\t",row.names = F)

unique(filt7$Recipient)
unique(filt7$Division)
unique(filt7$Class)


recip<-data.frame(table(filt7$Recipient))
div<-data.frame(table(filt7$Division))
class<-data.frame(table(filt7$Class))
Species<-data.frame(table(filt7$Species))

merg1<-dplyr::bind_rows(list(A_Recipient=recip,B_Division=div,Class=class,D_Species=Species),.id="Facet")
head(merg1)
colnames(merg1)<-c("Facet","y_list","Count")
write.table(merg1,"overall_plot.tsv",quote=F,sep="\t",row.names = F)

merg1<-read.delim("overall_plot.tsv")
unique(merg1$Facet)


physac<-ggplot(merg1[which(merg1$Facet=="A_Recipient"),])+geom_bar(aes(x=Count,y=y_list,fill=Color),alpha=0.5,stat="identity")+
  theme_minimal()+scale_size(range=c(0,1))+scale_fill_manual(values="#1B9E77")+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_discrete(expand=c(0,1))
physac


aa<-ggplot(merg1[which(merg1$AI>200),])+geom_point(aes(x="",y=fct_reorder(y_list,Count),size=Count,color=Color,alpha=0.5))+
  theme_minimal()+scale_size(range=c(0,30))+scale_color_brewer(palette = "Spectral")+
  facet_row(vars(Facet),scales="free",space="free")+
  geom_text(aes(x="",y=y_list,label=Count,size=40))+
  theme(legend.position = "none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),axis.text.y =element_text(size=12,color="black"),strip.text.x = element_blank())+
  scale_y_discrete(expand=c(0,2))+xlab("Overall HGT donors")+ylab("")
aa
physac+aa
#geom_text(aes(x="",y=y_list,label=ifelse(Count>=50,Count,''),size=10))
#check the distribution of Horizontal transfers (HT) and Vertical transfers (VT) in each species
head(filt4)

ggplot(filt4,aes(x=AIlog10_basido_asco,color=Typelog10,fill=Typelog10))+
  geom_histogram(aes(y=log10(..count..)),alpha=0.5,bins=150)+facet_grid(Recipient~.,scales = "free")+
  theme_minimal()+geom_vline(xintercept = 0,linetype="dashed")+
  scale_fill_manual(values=c("steelblue","tan"))+scale_color_manual(values=c("steelblue","tan"))+
  theme(strip.text.y = element_text(angle = 0))+
  xlab("Alien index Basdio vs Asco")+
  ylab("log10_transformed_counts")

ggplot(filt4,aes(x=log10AI_basido_asco,color=Typelog10,fill=Type))+
  geom_histogram(aes(y=log10(..count..)),alpha=0.5,bins=150)+facet_grid(Recipient~.,scales = "free")+
  theme_minimal()+geom_vline(xintercept = 0,linetype="dashed")+
  scale_fill_manual(values=c("steelblue","tan"))+scale_color_manual(values=c("steelblue","tan"))+
  theme(strip.text.y = element_text(angle = 0))+
  xlab("Alien index Basdio vs Asco")+
  ylab("log10_transformed_counts")
##
head(filt4)

dupl<-read.delim("Armi_dup_clusters.tsv")
head(dupl)
colnames(dupl)<-c("OG_ID67","Armi")
filt5<-left_join(filt4,dupl,by="OG_ID67")
head(filt5)
filt5$Armi<-filt5$Armi %>% replace_na('Not_dupl')

ggplot(filt5)+geom_point(aes(x=AI_basido_asco,y=Recipient,color=Armi))+
  scale_color_manual(values=c("darkgreen","grey"))+theme_minimal()+
  geom_vline(xintercept = 0,linetype="dashed")

#######################################################

############_________For bitscores_________##############

#######################################################

#######Read the bitscore file##########
big<-fread("bitscore_all.tsv")
head(big)
big$query_sps<-sub('_[^_]*$','',big$query)
###remove hits if query and subject are from the same species
big1<-big[(big$query_sps!=big$Species),]
head(big1)
big2<-big1 %>% 
  group_by(query,Category) %>%
  arrange(desc(bitscore)) %>%
  slice(1) %>%
  ungroup()
head(big2)
write.table(big2,"bitscore_hits.tsv",quote=F,sep="\t",row.names = F)


#################################################################################################

##############________FOR HGT INDEX CALCULATION START FROM HERE________##############

#################################################################################################

filt<-read.delim("bitscore_hits.tsv")
head(filt)
filt1<-filt[c(1,5,3)]
head(filt1)
filt2<-unique(filt1)
head(filt2)
tail(filt2)
min(filt2$bitscore)
filt3<-dcast(filt2,query~Category)
head(filt3)
##replace NA by 1
filt3[is.na(filt3)] <- 1
head(filt3)
###AI calculation
filt3$HI_basido_asco<-filt3$Ascomycota-filt3$Basidiomycota
head(filt3)


#add anot and OG_IDs from 67 sps dataset
anot<-read.delim("67sps_OGID_ProtID_IPRs.tsv")
head(anot)
unique(anot$Species)
colnames(anot)<-c("OG_ID67","query","Species","IPR","Desc")
filt4<-left_join(filt3,anot,by="query")
head(filt4)
unique(filt4$Species)
filt4$Type<-ifelse(filt4$HI_basido_asco>0,"HT","VT")
head(filt4)
filt4$Recipient<-sub('_[^_]*$','',filt4$query)
write.table(filt4,"HGTindex_bitscore_calc.tsv",quote=F,sep="\t",row.names = F)

#check the number of Horizontal transfers (HT) and Vertical transfers (VT) in each species
head(filt4)
c<-ggplot(filt4)+geom_bar(aes(y=Recipient,fill=Type))+
  facet_wrap(~Type, scales="free")+
  scale_fill_manual(values=c("steelblue","tan"))+theme_classic()+xlab("# of genes_BITSCORE")
c

a+c

#check the distribution of Horizontal transfers (HT) and Vertical transfers (VT) in each species
head(filt4)

ggplot(filt4,aes(x=HI_basido_asco,color=Type,fill=Type))+
  geom_histogram(aes(y=log10(..count..)),alpha=0.5,bins=200)+facet_grid(Recipient~.,scales = "free")+
  theme_minimal()+geom_vline(xintercept = 0,linetype="dashed")+
  scale_fill_manual(values=c("steelblue","tan"))+scale_color_manual(values=c("steelblue","tan"))+
  theme(strip.text.y = element_text(angle = 0))+
  xlab("HGT index Basdio vs Asco")+
  ylab("log10_transformed_counts")


###CHECK how many genes from the Armi clade expansions are HGT
exp<-read.delim("Armi_dup_clusters.tsv")
head(exp)
colnames(exp)<-c("OG_ID67","Status")
length(exp$OG_ID67)

ht<-read.delim("AI_evalue_calc.tsv")
head(ht)
unique(ht$Recipient)
physac<-c("Armbor1","Armga1","Armme1_1","Armnab1","Armost1","Armosto1",
"Armtab1","Flave","Hymrad1","Oudmuc1","Armect1",
"Cylto1","Flaro","Guyne1","Armfus","Armnov1",
"Armlut1","Armmel1","Armcep1","Armfum1")
ht1<-filter(ht,ht$Recipient %in% physac)

htt<-left_join(ht1,exp,by="OG_ID67")
head(htt)
htt$Status<-gsub('ARMI', 'Expansion in Armillaria clade', htt$Status)
htt$Status[is.na(htt$Status)] <- c("Not_expanded")


ggplot(htt,aes(y=Recipient,fill=Status))+geom_bar()+facet_wrap(~Typelog10,scales="free")+
  theme_classic()+scale_fill_manual(values=c("tan","steelblue"))+xlab("Gene Transfer and Expansions in Armillaria MRCA")

ggplot(ht1,aes(y=Recipient,fill=Status))+geom_bar()+facet_wrap(~Typelog10,scales="free")+
  theme_classic()+scale_fill_manual(values=c("tan","steelblue"))+xlab("Gene Transfer and Expansions in Armillaria MRCA")


ht1<-ht[c(1,2,6)]
head(ht1)
colnames(ht1)<-c("ProtID","AI","Type")

ht2<-left_join(dup,ht1,by="ProtID")
head(ht2)
length(unique(ht2$Ortho67))
sapply(ht2,class)

bar<-data.frame(table(ht2$Species))
head(bar)
bar1<-data.frame(table(ht2$Species,ht2$Type))
head(bar1)
bar2<-left_join(bar1,bar,by="Var1")
head(bar2)
bar3<-bar2[c(1,4,3)]
colnames(bar3)<-c("Species","Total","HT")
head(bar3)
bar3$Other<-bar3$Total-bar3$HT
head(bar3)
bar3_1<-filter(bar3,bar3$Species %in% ht$query_sps)
head(bar3_1)

ggplot(bar3_1,aes(y=Species,x=Total))+
  geom_tile()

head(bar3_1)
bar4<-melt(bar3_1[c(1,2,3)])
head(bar4)
tail(bar4)
bar4$log10<-log10(bar4$value)
ggplot(bar4,aes(y=Species,x=log10,fill=variable))+
  geom_point()
  
