library(ape)
library(RRphylo)
library(PDcalc)
library(phytools)
library(parallel)

rm(list=ls())
setwd("C:/Users/Dell/Documents/TopGO/Armi_compgen_updated/HGT/plant_bact/all_fun_plant_bac/generax_crisscross_Tree/polytomy_resolver/polytomic_trees/")
list.files()


#eg<-read.newick("OG0005935.newick")
#plot(eg)
#eg1<-bifurcatr(eg,runs=1)
#plot(eg1)

#
wl=list.files(pattern=".newick")
wl=list.files(pattern=".fas.treefile")

wl
names(wl)=gsub("\\.fas.treefile$","",wl)
wl

myfunc=function(treename)
{
  treefile=wl[[treename]]
  my_tree=read.newick(treefile)
  outfile=paste0(treename,".tree")
  bif_res=bifurcatr(my_tree, runs=1)
  write.tree(bif_res,outfile)
  return(outfile)
}
#debug(myfunc)
#myfunc("OG0005935")
#file.info("OG0005935_bf.tree")
res=mclapply(names(wl),FUN=myfunc)
res


