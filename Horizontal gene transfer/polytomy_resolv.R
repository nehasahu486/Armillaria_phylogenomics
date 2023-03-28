############################################################################################################################################
#############################  Resolve gene trees polytomies (nodes with more than two descendent branches)    #############################
############################################################################################################################################

##Used bifurcatr function from PDcalc package (https://rdrr.io/github/davidnipperess/PDcalc/man/bifurcatr.html) which resolves all polytomies randomly
library(ape)
library(RRphylo)
library(PDcalc)
library(phytools)
library(parallel)

##specify the directory that has polytomic gene trees

poly_tree=list.files(pattern=".newick")
poly_tree
names(poly_tree)=gsub("\\.fas.treefile$","",wl)
poly_tree


myfunc=function(treename)
{
  treefile=poly_tree[[treename]]
  my_tree=read.newick(treefile)
  outfile=paste0(treename,".tree")
  bif_res=bifurcatr(my_tree, runs=1)
  write.tree(bif_res,outfile)
  return(outfile)
}

res=mclapply(names(poly_tree),FUN=myfunc)
res

##the resolved gene trees will be save in the same original directory, with the suffix ".tree"
