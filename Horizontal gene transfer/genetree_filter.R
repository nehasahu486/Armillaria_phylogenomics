###########################################################################################################################################################
#############################   Filter potential HT clades based on support values and taxonomic composition in gene trees    #############################
###########################################################################################################################################################

##This script was orginally written by Zsolt Merenyi (https://github.com/zsmerenyi)

rm(list=ls())

library(DECIPHER)
library(ape)
library(readr)
library(parallel)
library(stringr)
library(phytools)
library(phangorn)
library(dplyr)

##Create the function with parameters to select clades from genetree

HGTDUP3<-function(fa,boot=0.7,ABG=0.6,AT=0.6,PT=0.7,abra=F){
  options(stringsAsFactors=FALSE)
  evt<-read.tree(fa)
  evt<-LadderizeTree(evt)
  cat(fa)
  anevek<-sapply(str_split(rev(str_split(fa,"\\/")[[1]])[1] ,"//."), "[", 1)
  anevek<-sapply(str_split(anevek ,"\\."), "[", 1)
  anevek<-gsub("_events","",anevek)
  #TFA<-TransFile[TransFile$File==anevek,]
  
  ##construct the shared possible subtree
  stv<-Sys.time()
  st<-subtrees(evt) 
  Sys.time()-stv
  
  names(st)<-as.numeric(length(evt$tip.label)+1:length(st))
  
  ET<-tibble(evt$node.label)
  colnames(ET)<-"support"
  ET$support<-as.numeric(ET$support)/100
  ET$node<-1:nrow(ET)+length(evt$tip.label)
  ETT<-ET 
  if(nrow(ETT)==0){ETTR<-NA}else{
    edt<-data.frame(evt$edge)
    edl<-lapply(ETT$node,function(h) edt[which(edt$X1==h),"X2"])
    edt0<-lfuzo(edl,1)
    ETT<-cbind(ETT,edt0)
    
    #x<-which(ETT$node==653) x=1
    gennevek<-function(x){
      #    gy<-NULL
      #for(x in 1:nrow(ETT)){
      mite<-ETT[x,]
      cl1<-st[as.character(mite$`1`)][[1]]
      cl2<-st[as.character(mite$`2`)][[1]]
      if(is.null(cl1)){
        c1sp<-tibble(protID=evt$tip.label[mite$`1`])
      }else{
        c1sp<-tibble(protID=cl1$tip.label)
      }
      if(is.null(cl2)){
        c2sp<-tibble(protID=evt$tip.label[mite$`2`])
      }else{
        c2sp<-tibble(protID=cl2$tip.label)
      }
      
      c1sp$sp<-spresz(c1sp$protID)
      c2sp$sp<-spresz(c2sp$protID)
      c1sp<-bind_cols(c1sp,TAXN[match(c1sp$sp,TAXN$nodelabel),"Category"])
      c2sp<-bind_cols(c2sp,TAXN[match(c2sp$sp,TAXN$nodelabel),"Category"])
      d1<-sum((cl1$node.label=="D")*1) # cladeon belul a D-k Duplikaciok szama
      d2<-sum((cl2$node.label=="D")*1)
      e1<-length(cl1$node.label)       # osszesen ennyi edge van a cladeon belul
      e2<-length(cl2$node.label)
      T1<-data.frame(table(c1sp$Category)) # kategoriak szama
      T2<-data.frame(table(c2sp$Category))
      #no of genes
      Asco1<-length(which(c1sp$Category %in% "Ascomycota"))
      Asco2<-length(which(c2sp$Category %in% "Ascomycota"))
      MZ1<-length(which(c1sp$Category %in% c("Mucoromycota","Zoopagomycota")))
      MZ2<-length(which(c2sp$Category %in% c("Mucoromycota","Zoopagomycota")))
      Phy1<-length(which(c1sp$Category %in% "Physac"))  
      Phy2<-length(which(c2sp$Category %in% "Physac"))
      G1<-sz(c1sp$protID)
      G2<-sz(c2sp$protID)
      sp1<-sz(c1sp$sp)
      sp2<-sz(c2sp$sp)
      #PGN1<-paste0(c1sp[which(c1sp$Category=="Physac"),"protID"][[1]],collapse="|")
      #PGN2<-paste0(c2sp[which(c2sp$Category=="Physac"),"protID"][[1]],collapse="|")
      
      PGN1<-paste0(c1sp$protID,collapse="|")
      PGN2<-paste0(c2sp$protID,collapse="|")
      
      dupP1<-length(which(duplicated(c1sp[which(c1sp$Category=="Physac"),"sp"])))
      dupP2<-length(which(duplicated(c2sp[which(c2sp$Category=="Physac"),"sp"])))
      
      gner<-c(mite$node,d1, d2, e1, e2, Asco1, Asco2, MZ1, MZ2, Phy1, Phy2, G1, G2, sp1, sp2, dupP1, dupP2,PGN1,PGN2)
      return(gner)
    }
    
    spstat0<-lapply(1:nrow(ETT),function(x) gennevek(x))
    spstat<-lfuzo(spstat0,1)
    colnames(spstat)<-c("node","d1", "d2", "e1", "e2", "Asco1", "Asco2", "MZ1", "MZ2", "Phy1", "Phy2", "G1", "G2", "sp1", "sp2", "dupP1", "dupP2","PGN1","PGN2")
    spstat<-chrtonum(spstat,2:17)  
    spstat$AA1<-(spstat$Asco1+spstat$MZ1)/spstat$G1
    spstat$AA2<-(spstat$Asco2+spstat$MZ2)/spstat$G2
    spstat$PHA1<-spstat$Phy1/spstat$G1
    spstat$PHA2<-spstat$Phy2/spstat$G2
    spstat
    #spstat[which(spstat$AA1>=0.5 & spstat$PHA2>=0.7),]
    #spstat[which(spstat$AA2>=0.5 & spstat$PHA1>=0.7),]
    
    x<-ETT$node[2]
    MAG<-lapply(ETT$node,function(x){ #for(x in ETT$node){
      if((length(evt$tip.label)+1)==x){hatter<-0}else{
        prevnode<-edt[which(edt$X2==x),1]
        masikag<-setdiff(edt[which(edt$X1==prevnode),"X2"],x)
        if(masikag<=length(evt$tip.label)){
          hspec<-spresz(evt$tip.label[masikag])
          H2<-TAXN[which(TAXN$nodelabel==hspec),"Category"]
          hatter<-sum((H2=="Ascomycota")*1)/nrow(H2)
        }else{
          hfa<-extract.clade(evt,masikag)
          hspec<-spresz(hfa$tip.label)
          H2<-TAXN[which(TAXN$nodelabel %in% hspec),"Category"]
          hatter<-sum((H2=="Ascomycota")*1)/nrow(H2)
        }
      }
      return(hatter)}
    )
    ETT$Ascobackgr<-unlist(MAG)
    
    ETT<-as_tibble(ETT)
    ETT$OG<-anevek
    colnames(ETT)[2]<-"TNODE"
    ETTR<-as_tibble(cbind(ETT,spstat))
    #save(ETTR,file=paste0(anevek,".saved"))
  
    write.table(ETTR, file=paste0('ETTR_',anevek,"_",SD,'.csv'), quote=FALSE, sep=',', row.names = F)
    
    #FILTER
    VETTR<-ETTR[which(ETTR$support>=boot & ETTR$Ascobackgr>=ABG),]
    VT1<-VETTR[which(VETTR$AA2>AT & VETTR$PHA1>PT),]
    VT2<-VETTR[which(VETTR$AA1>AT & VETTR$PHA2>PT),]
    EFI<-rbind(VT1,VT2)
  }
  
  if(abra){
    gtt<-as_tibble(data.frame(evt$tip.label))
    colnames(gtt)<-"protID"
    gtt$spec<-spresz(gtt$protID)
    gtt$col<-"#000000"
    ascok<-TAXN[which(TAXN$Category=="Ascomycota"),"nodelabel"][[1]]
    physek<-TAXN[which(TAXN$Category=="Physac"),"nodelabel"][[1]]
    gtt[which(gtt$spec %in% ascok), "col"]<-"#FF0000"	
    gtt[which(gtt$spec %in% physek), "col"]<-"#8A2BE2"	
    
    meret<-length(evt$tip.label)/10
    if(meret<=30){meret2<-1}else{meret2<-30/meret}
    
      File2<-paste0(anevek,SD,".pdf")
      pdf(file = File2, width=20,height=10+meret, pointsize = 18, family = "sans", bg = "white")
      plot.phylo(evt,tip.color=gtt$col,cex=0.6)
      nodelabels(evt$node.label,frame="none",col="blue",cex=1.0*meret2)
      nodelabels(frame="none",adj=c(1.5,0),cex=0.8*meret2)
      dev.off()
    
  }
  
  return(EFI)
}
###################
###################
###################
###################
###################
###################
###################
###################
HGTDUPPhys<-function(fa,boot=0.7,ABG=0.6,AT=0.6,PT=0.7,abra=F){
  options(stringsAsFactors=FALSE)
  evt<-read.tree(fa)
  evt<-LadderizeTree(evt)
  cat(fa)
  anevek<-sapply(str_split(rev(str_split(fa,"\\/")[[1]])[1] ,"//."), "[", 1)
  anevek<-sapply(str_split(anevek ,"\\."), "[", 1)
  anevek<-gsub("_events","",anevek)
  #TFA<-TransFile[TransFile$File==anevek,]
  
  # legyartjuk az osszes lehetseges subtreet
  stv<-Sys.time()
  st<-subtrees(evt) 
  Sys.time()-stv
  
  names(st)<-as.numeric(length(evt$tip.label)+1:length(st))
  
  ET<-tibble(evt$node.label)
  colnames(ET)<-"support"
  ET$support<-as.numeric(ET$support)/100
  ET$node<-1:nrow(ET)+length(evt$tip.label)
  ETT<-ET 
  if(nrow(ETT)==0){ETTR<-NA}else{
    edt<-data.frame(evt$edge)
    edl<-lapply(ETT$node,function(h) edt[which(edt$X1==h),"X2"])
    edt0<-lfuzo(edl,1)
    ETT<-cbind(ETT,edt0)
    
    #x<-which(ETT$node==653) x=1
    gennevek<-function(x){
      #    gy<-NULL
      #for(x in 1:nrow(ETT)){
      mite<-ETT[x,]
      cl1<-st[as.character(mite$`1`)][[1]]
      cl2<-st[as.character(mite$`2`)][[1]]
      if(is.null(cl1)){
        c1sp<-tibble(protID=evt$tip.label[mite$`1`])
      }else{
        c1sp<-tibble(protID=cl1$tip.label)
      }
      if(is.null(cl2)){
        c2sp<-tibble(protID=evt$tip.label[mite$`2`])
      }else{
        c2sp<-tibble(protID=cl2$tip.label)
      }
      
      c1sp$sp<-spresz(c1sp$protID)
      c2sp$sp<-spresz(c2sp$protID)
      c1sp<-bind_cols(c1sp,TAXN[match(c1sp$sp,TAXN$nodelabel),"Category"])
      c2sp<-bind_cols(c2sp,TAXN[match(c2sp$sp,TAXN$nodelabel),"Category"])
      d1<-sum((cl1$node.label=="D")*1) # cladeon belul a D-k Duplikaciok szama
      d2<-sum((cl2$node.label=="D")*1)
      e1<-length(cl1$node.label)       # osszesen ennyi edge van a cladeon belul
      e2<-length(cl2$node.label)
      T1<-data.frame(table(c1sp$Category)) # kategoriak szama
      T2<-data.frame(table(c2sp$Category))
      #no of genes
      Asco1<-length(which(c1sp$Category %in% "Ascomycota"))
      Asco2<-length(which(c2sp$Category %in% "Ascomycota"))
      MZ1<-length(which(c1sp$Category %in% c("Mucoromycota","Zoopagomycota")))
      MZ2<-length(which(c2sp$Category %in% c("Mucoromycota","Zoopagomycota")))
      Phy1<-length(which(c1sp$Category %in% "Physac"))  
      Phy2<-length(which(c2sp$Category %in% "Physac"))
      G1<-sz(c1sp$protID)
      G2<-sz(c2sp$protID)
      sp1<-sz(c1sp$sp)
      sp2<-sz(c2sp$sp)
      #PGN1<-paste0(c1sp[which(c1sp$Category=="Physac"),"protID"][[1]],collapse="|")
      #PGN2<-paste0(c2sp[which(c2sp$Category=="Physac"),"protID"][[1]],collapse="|")
      
      PGN1<-paste0(c1sp$protID,collapse="|")
      PGN2<-paste0(c2sp$protID,collapse="|")
      
      dupP1<-length(which(duplicated(c1sp[which(c1sp$Category=="Physac"),"sp"])))
      dupP2<-length(which(duplicated(c2sp[which(c2sp$Category=="Physac"),"sp"])))
      
      gner<-c(mite$node,d1, d2, e1, e2, Asco1, Asco2, MZ1, MZ2, Phy1, Phy2, G1, G2, sp1, sp2, dupP1, dupP2,PGN1,PGN2)
      return(gner)
    }
    
    spstat0<-lapply(1:nrow(ETT),function(x) gennevek(x))
    spstat<-lfuzo(spstat0,1)
    colnames(spstat)<-c("node","d1", "d2", "e1", "e2", "Asco1", "Asco2", "MZ1", "MZ2", "Phy1", "Phy2", "G1", "G2", "sp1", "sp2", "dupP1", "dupP2","PGN1","PGN2")
    spstat<-chrtonum(spstat,2:17)  
    spstat$AA1<-(spstat$Asco1+spstat$MZ1)/spstat$G1
    spstat$AA2<-(spstat$Asco2+spstat$MZ2)/spstat$G2
    spstat$PHA1<-spstat$Phy1/spstat$G1
    spstat$PHA2<-spstat$Phy2/spstat$G2
    spstat
    #spstat[which(spstat$AA1>=0.5 & spstat$PHA2>=0.7),]
    #spstat[which(spstat$AA2>=0.5 & spstat$PHA1>=0.7),]
    
    x<-ETT$node[2]
    MAG<-lapply(ETT$node,function(x){ #for(x in ETT$node){
      if((length(evt$tip.label)+1)==x){hatter<-0}else{
        prevnode<-edt[which(edt$X2==x),1]
        masikag<-setdiff(edt[which(edt$X1==prevnode),"X2"],x)
        if(masikag<=length(evt$tip.label)){
          hspec<-spresz(evt$tip.label[masikag])
          H2<-TAXN[which(TAXN$nodelabel==hspec),"Category"]
          hatter<-sum((H2=="Ascomycota")*1)/nrow(H2)
        }else{
          hfa<-extract.clade(evt,masikag)
          hspec<-spresz(hfa$tip.label)
          H2<-TAXN[which(TAXN$nodelabel %in% hspec),"Category"]
          hatter<-sum((H2=="Ascomycota")*1)/nrow(H2)
        }
      }
      return(hatter)}
    )
    ETT$Ascobackgr<-unlist(MAG)
    
    ETT<-as_tibble(ETT)
    ETT$OG<-anevek
    colnames(ETT)[2]<-"TNODE"
    ETTR<-as_tibble(cbind(ETT,spstat))
    #save(ETTR,file=paste0(anevek,".saved"))
    
    write.table(ETTR, file=paste0('ETTR_',anevek,"_",SD,'.csv'), quote=FALSE, sep=',', row.names = F)
    
    #FILTER
    VETTR<-ETTR[which(ETTR$support>=boot & ETTR$Ascobackgr>=ABG),]
    VT1<-VETTR[which(VETTR$AA2>AT & VETTR$PHA1>PT),]
    VT2<-VETTR[which(VETTR$AA1>AT & VETTR$PHA2>PT),]
    EFI<-rbind(VT1,VT2)
  }
  
  if(abra){
    gtt<-as_tibble(data.frame(evt$tip.label))
    colnames(gtt)<-"protID"
    gtt$spec<-spresz(gtt$protID)
    gtt$col<-"#000000"
    ascok<-TAXN[which(TAXN$Category=="Ascomycota"),"nodelabel"][[1]]
    physek<-TAXN[which(TAXN$Category=="Physac"),"nodelabel"][[1]]
    gtt[which(gtt$spec %in% ascok), "col"]<-"#FF0000"	
    gtt[which(gtt$spec %in% physek), "col"]<-"#8A2BE2"	
    
    meret<-length(evt$tip.label)/10
    if(meret<=30){meret2<-1}else{meret2<-30/meret}
    
    File2<-paste0(anevek,SD,".pdf")
    pdf(file = File2, width=20,height=10+meret, pointsize = 18, family = "sans", bg = "white")
    plot.phylo(evt,tip.color=gtt$col,cex=0.6)
    nodelabels(evt$node.label,frame="none",col="blue",cex=1.0*meret2)
    nodelabels(frame="none",adj=c(1.5,0),cex=0.8*meret2)
    dev.off()
    
  }
  
  return(EFI)
}
#----------------------------------------------------------------------------------------------------------
source("SFun220218.R")


list.files()

stree<-read.tree("inferred_species_tree.newick") 

GFlist<-list.files(path="genetrees/", pattern = "MT.tree",full.names = T)
GFlist
GFT<-as_tibble(GFlist)
GFT$OG<-gsub("MT.tree","",gsub("genetrees/","",GFT$value))
TAXN<-read_tsv("Nodes.tsv")

fa<-GFT[grep("OG0001085",GFT$OG),1][[1]] 
fa
fa<-GFT$value[[2]]
fa
#HGTDUP3(fa)



testf<-lapply(GFT$value,function(fa) HGTDUP3(fa,boot=0.7,ABG=0.7,AT=0.7,PT=0.7,abra=T))
testf
testt<-lfuzo(testf,1)
write.table(testt, file=paste0('TestT',SD,'.tsv'), quote=FALSE, sep='\t', row.names = F)





# for debug mode:
for(i in GFT$value){
  cat(i,"\n")
  HGTDUP3(fa,boot=80,ABG=0.6,AT=0.6,PT=0.7)
}
