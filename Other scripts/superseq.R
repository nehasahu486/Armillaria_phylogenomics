#library(devtools)
#install_github("StoreyLab/superSeq", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
#3
library(superSeq)
library(subSeq)
library(Biobase)
rm(list=ls())

 
#the counts file has been already formated and filtered to the following parameters 
#(the same parameters used for DEG calculation):
#remove rowsum 0
#remove lowly expressed genes, keep transcripts with counts >=5 in atleast 4 samples
#remove JW18 sample - outlier in MDS plot
#counts <- data.matrix(counts[rowSums(counts) >= 5, ]) #since we have prefiltered, we do not need this step

counts<-read.delim("Armlut_rawreads",row.names = 1) 
head(counts)
tail(counts)
class(counts)
length(rownames(counts)) #check how many genes it has, for Armlut it should have 9318 genes

design<-read.delim("Armlut_design",row.names = 1)
head(design)

options("max.print"=10000000)

#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
####for 24hrs-Control samples
counts1<-counts[,(1:8)] #take the columns of the control and 24 hrs samples 
head(counts1) #check


samp<-c("Control","24hrs") #list the pair which you need to compare

design1<-filter(design,Sample%in%samp)
design1
tail(design1)

#proportions: what you would like to check, eg if you have 60M library depth then this is proportion 1 
#and in this simulation you are inquiring what happens if you would have just 30M read (0.5 proportion)
#or 6M (0.1 proportion)
#proportions <- 10 ^ seq(-2,1,0.01)
proportions <- seq(0.01,1,0.01)
proportions

ss = subsample(counts = counts1,treatment=design1$Sample,proportions = proportions,method=c("voomLimma"),seed=12345,replications=3)


ss_sum <- summary(ss)
head(ss_sum)
write.table(ss_sum,"24hrs_control_subsample.tsv",sep="\t",quote=F,row.names = F)
# apply superSeq model
ss_obj <- superSeq(ss_sum)
pred<-ss_obj$predictions
head(pred)
pred$Sample<-c("24hrs_control")
pred
write.table(pred,"24hrs_controlpred.tsv",sep="\t",quote=F,row.names = F)
#capture.output(ss_obj,file="24hrs_control_sum.txt")
# plot results
plot(ss_obj,main="24hrs_control")

# summarise results
summary(ss_obj)
capture.output(summary(ss_obj),file="24hrs_control_summary.txt")
#save(list=ls(),file="zsnal.saved")



#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
####for 24hrs_memb-Control samples
head(counts)
counts1<-counts[,c(1:4,9:12)] #take the columns of the control and 24 hrs samples 
head(counts1) #check


samp<-c("Control","24hrs_memb") #list the pair which you need to compare

design1<-filter(design,Sample%in%samp)
design1
tail(design1)

proportions <- seq(0.01,1,0.01)
proportions

ss = subsample(counts = counts1,treatment=design1$Sample,proportions = proportions,method=c("voomLimma"),seed=12345,replications=3)

ss_sum <- summary(ss)
head(ss_sum)
write.table(ss_sum,"24hrs_memb_control_subsample.tsv",sep="\t",quote=F,row.names = F)
# apply superSeq model
ss_obj <- superSeq(ss_sum)
pred<-ss_obj$predictions
head(pred)
pred$Sample<-c("24hrs_memb_control")
pred
write.table(pred,"24hrs_memb_controlpred.tsv",sep="\t",quote=F,row.names = F)
# plot results
plot(ss_obj,main="24hrs_memb_control")

# summarise results
summary(ss_obj)
capture.output(summary(ss_obj),file="24hrs_memb_control_summary.txt")


#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
####for 48hrs-Control samples
head(counts)
counts1<-counts[,c(1:4,13:16)] #take the columns of the control and 24 hrs samples 
head(counts1) #check


samp<-c("Control","48hrs") #list the pair which you need to compare

design1<-filter(design,Sample%in%samp)
design1
tail(design1)

proportions <- seq(0.01,1,0.01)
proportions

ss = subsample(counts = counts1,treatment=design1$Sample,proportions = proportions,method=c("voomLimma"),seed=12345,replications=3)

ss_sum <- summary(ss)
head(ss_sum)
write.table(ss_sum,"48hrs_control_subsample.tsv",sep="\t",quote=F,row.names = F)
# apply superSeq model
ss_obj <- superSeq(ss_sum)
pred<-ss_obj$predictions
head(pred)
pred$Sample<-c("48hrs_control")
pred
write.table(pred,"48hrs_controlpred.tsv",sep="\t",quote=F,row.names = F)
# plot results
plot(ss_obj,main="48hrs_control")

# summarise results
summary(ss_obj)
capture.output(summary(ss_obj),file="48hrs_control_summary.txt")


#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
####for 1week-Control samples
head(counts)
counts1<-counts[,c(1:4,17:19)] #take the columns of the control and 24 hrs samples 
head(counts1) #check


samp<-c("Control","1week") #list the pair which you need to compare

design1<-filter(design,Sample%in%samp)
design1
tail(design1)

proportions <- seq(0.01,1,0.01)
proportions

ss = subsample(counts = counts1,treatment=design1$Sample,proportions = proportions,method=c("voomLimma"),seed=12345,replications=3)

ss_sum <- summary(ss)
head(ss_sum)
write.table(ss_sum,"1week_control_subsample.tsv",sep="\t",quote=F,row.names = F)
# apply superSeq model
ss_obj <- superSeq(ss_sum)
pred<-ss_obj$predictions
head(pred)
pred$Sample<-c("1week_control")
pred
write.table(pred,"1week_controlpred.tsv",sep="\t",quote=F,row.names = F)
# plot results
plot(ss_obj,main="1week_control")

# summarise results
summary(ss_obj)
capture.output(summary(ss_obj),file="1week_control_summary.txt")


#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
####for 2week-Control samples
head(counts)
counts1<-counts[,c(1:4,20:23)] #take the columns of the control and 24 hrs samples 
head(counts1) #check


samp<-c("Control","2week") #list the pair which you need to compare

design1<-filter(design,Sample%in%samp)
design1
tail(design1)

proportions <- seq(0.01,1,0.01)
proportions

ss = subsample(counts = counts1,treatment=design1$Sample,proportions = proportions,method=c("voomLimma"),seed=12345,replications=3)

ss_sum <- summary(ss)
head(ss_sum)
write.table(ss_sum,"2week_control_subsample.tsv",sep="\t",quote=F,row.names = F)
# apply superSeq model
ss_obj <- superSeq(ss_sum)
pred<-ss_obj$predictions
head(pred)
pred$Sample<-c("2week_control")
pred
write.table(pred,"2week_controlpred.tsv",sep="\t",quote=F,row.names = F)
# plot results
plot(ss_obj,main="2week_control")


# summarise results
summary(ss_obj)
capture.output(summary(ss_obj),file="2week_control_summary.txt")

####################################################################################
####################################################################################
####################################################################################
####################################################################################
##rbind prediction files for all samples

predfiles = purrr::map(Sys.glob("*pred.tsv"), read.delim) 
pf1<-data.table::rbindlist(predfiles)
head(pf1)
tail(pf1)
unique(pf1$Sample)
pf1$power<-ceiling(pf1$estimated_power*100)
head(pf1)

unique(pf1$Sample)

#saved pdf at 9*8
head(pf1)
p<-ggplot(pf1[which(pf1$method=="voomLimma"),],aes(x=proportion, y=power))+
  geom_vline(xintercept = 1, linetype="dashed",color="grey",size=0.8)+
  geom_line(aes(x=proportion, y=power,color=Sample),size=1,alpha=0.7)+
  geom_point(aes(x=proportion, y=power,color=Sample),alpha=0.7)+
  theme_bw()+
  xlab("Proportion of read depth")+ylab("Percentage of predicted DEGs")+
  ggtitle("Transcriptomics")+scale_color_brewer(palette = "Set1")
p

p1<-ggplot(pf1[which(pf1$method=="voomLimma"),],aes(x=proportion, y=predicted))+
  geom_vline(xintercept = 1, linetype="dashed",color="grey",size=0.8)+
  geom_line(aes(x=proportion, y=predicted,color=Sample),alpha=0.7)+
  geom_point(aes(x=proportion, y=predicted, color=Sample),alpha=0.7)+
  theme_bw()+
  xlab("Proportion of read depth")+ylab("Number of predicted DEGs")+
  ggtitle("Transcriptomics")+scale_color_brewer(palette = "Set1")
p1
