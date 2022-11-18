############################################################################################################
#¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯#
#________Prediction of Secretome and Small Secreted Proteins (SSPs) in 131 species (Dataset 2)_____________#
#__________________________________________________________________________________________________________#
############################################################################################################
¯

###STEP1 : Concetenate the proteome files of 131 species into one fasta file (linux)
#cat *.fas > 131sps_all.fas (here 131sps_all.fas is the output file) -this file needs checking before further processing
#grep -c ">" 131sps_all.fas (checks the number of lines containing ">" character)
#grep -c "^>" 131sps_all.fas (checks the number of lines starting with ">", this should be same as the previous step - if not - CHECK why!) 
#once it is confirmed that fasta headers are in correct number proceed with the following steps


###STEP2 : Read the concatenated proteome fasta for 131 species. For this use "seqinr" package in R
setwd("my/directory")
library(seqinr)
allprot<-read.fasta("131sps_all.fas")
#next step will convert the big fasta file into 3 columns - ProtID, Sequence, and Length
protlen<-data.frame(ProtID=names(allprot), Sequence=unlist(getSequence(allprot, as.string=T)), Length=(getLength(allprot)))
write.table(protlen,"ProtIDs_len_all.tsv", quote=F, sep="\t", row.names = F)

#filter the ProtIDs which have Length <=300 amino acids; these will be SSPs
len300<-data.frame(protlen[protlen$Length<=300,])
write.table(len300,"ProtIDs_len300.tsv", quote=F, sep="\t", row.names = F)

###STEP3 : Make fasta file for all sequences with <=300aa using the ProtID_len300 file - this will be used in SignalP
len300<-read.delim("ProtIDs_len300.tsv")
head(len300)
protfasta300<-allprot[names(allprot) %in% len300$ProtID]
head(protfasta300)
write.fasta(sequences=protfasta300,names=names(protfasta300), file.out="fname.fasta", open="w",as.string=F)

##STEP4 : Run SignalP v.4.1 with the "fname.fasta" file (for SSPs) OR "131sps_all.fas" (for complete Secretome)
##select the -short output style - for one line per protein type output
###___command for SignalP__### /my/directory/signalp-4.1/signalp -t euk -f short fname.fasta > fname_short_out
##"fname_short_out" is the SignalP output file
#The same can be done for -summary and long type output as well - these would be the detailed output files

##STEP5: Use the output from SignalP to create a fasta file of sequences WITH a signal peptide
#this fasta file will be used in TMHMM to predict the presence of transmembrane domains
#In the columns provided in the SignalP ouput there will be a column "?" , which has the predictions Y and N - yes and no for signal peptides
#I changed the "?" to "Pred" beforehand and selected the ones which were marked "Y"
sigP_seq<-read.fasta("fname.fasta")
sigP<-data.frame(read.delim("fname_short_out"))
colnames(sigP)
levels(sigP$Pred)
sigP_Y<-data.frame(sigP[sigP$Pred=="Y",])
forTM<-sigP_seq[names(sigP_seq)%in%sigP_Y$name]
head(forTM)
write.fasta(sequences=forTM,names=names(forTM),file.out="forTM.fasta",open='w',as.string=F)

##STEP6: Use the "forTM.fasta" file to run TMHMM v2.0
###using -short output style - for one line per protein type output
###___command for TMHMM___### /my/directory/tmhmm-2.0c/bin/tmhmm forTM.fasta -short > TM_output.tsv
##"TM_output.tsv" is the output file from TMHMM

##STEP7: Use the output from TMHMM to create a fasta file of sequences with ZERO trasmembrane helices
#In the columns provided in the TMHMM ouput there will be a column "PredHel", with number of Transmembrane helices predicted in each sequence
#selected sequences with "0" transmembrane helices
#this fasta file will be used in WolfPsort v0.2 to predict the subcellular localization
TM_seq<-read.fasta("forTM.fasta")
head(TM_seq)
TM<-data.frame(read.delim("TM_output.tsv"))
head(TM)
colnames(TM)
levels(TM$PredHel)
TM_0<-data.frame(TM[TM$PredHel==0,])
forWP<-TM_seq[names(TM_seq)%in%TM_0$ProtID]
head(forWP)
write.fasta(sequences=forWP,names=names(forWP), file.out="forWP.fasta",open='w',as.string=F)

##STEP8: Use the "forWP.fasta" file to run WolfPsort v0.2
###___command for WolfPsort___### /my/directory/wolfpsort/bin/runWolfPsortSummary fungi < forWP.fasta >  WP_out.tsv
#Use the output from WolfPsort i.e. "WP_out.tsv" to select sequences with extracellular localization ("extr")
#This will be the final list of SSPs (or Secreted proteins if "131sps_all.fas" was used at STEP4)

