SD<-format(Sys.time(),'_%m%d')

suppressMessages(library(readr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(tibble))
suppressMessages(library(RColorBrewer))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
set.seed(2);col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
col_vector <- unique(col_vector[c(unique(c(44,45,47,46,26,20,24,54,58,22,19,18,8,56,7,25,4,38,35,34,31,30,29,27,6,2,1,42)),c(1:60)[-which(c(1:60) %in% unique(c(47,46,26,20,24,54,58,22,19,18,8,56,7,25,4,38,35,34,31,30,29,27,6,2,1,42)))])] )


LadderizeTree <- function(tree, temp_file = "temp", orientation = "left"){
  if(file.exists(paste0("./", temp_file))){
    stop("The chosen temporary file exists! Please choose an other temp_file name")
  }
  if(orientation == "left"){
    right <- FALSE
  }else{
    right <- TRUE
  }
  tree_temp <- ladderize(tree, right = right)
  write.tree(tree_temp, file = paste0("./", temp_file, ".tre"))
  tree_lad <- read.tree(paste0("./", temp_file, ".tre"))
  file.remove(paste0("./", temp_file, ".tre"))
  return(tree_lad)
}

spresz<-function(s){
  s2<-str_split(s,"_")
  s3<-sapply(1:length(s2), function(x) paste0(s2[[x]][1:(length(s2[[x]])-1)],collapse = "_") )
  return(s3)
}

sz<-function(vmi){length(unique(vmi))}

chrtonum<-function(tabla,oszl=c(1:1)){
  orignames<-colnames(tabla)
  for(h in oszl){
    tabla[,h]<-as.numeric(tabla[,h][[1]])
  }
  colnames(tabla)<-orignames
  return(tabla)
}

# Wraps a single sentence
wrap_sentence <- function(string, width) {
  words <- unlist(strsplit(string, " "))
  fullsentence <- ""
  checklen <- ""
  for(i in 1:length(words)) {
    checklen <- paste(checklen, words[i])
    if(nchar(checklen)>(width+1)) {
      fullsentence <- paste0(fullsentence, "\n")
      checklen <- ""
    }
    fullsentence <- paste(fullsentence, words[i])
  }
  fullsentence <- sub("^\\s", "", fullsentence)
  fullsentence <- gsub("\n ", "\n", fullsentence)
  return(fullsentence)
}


t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

lfuzo<-function(x,ncpu){
  NL<-x #<-resout
  ll<-length(x)
  tag<-round(ll/ncpu,0)
  a=1
  LL<-list()
  for(i in 1:ncpu){
    b=a+(tag-1)
    LL[[i]]<-c(i,a,b)
    a=b+1
  }
  x<-LL[[1]]
  rag<-function(x){
    aktl<-NL[x[2]:x[3]]
    rag<-plyr::ldply(aktl,rbind)
    return(rag)
  }
  er<-mclapply(LL,FUN=rag,mc.cores = ncpu)
  errag<-plyr::ldply(er,rbind)
  errag<-tibble::as_tibble(errag)
  return(errag) 
}        

