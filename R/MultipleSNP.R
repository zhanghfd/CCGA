MultipleSNP <-
function(Gs,Y,Z,S,fs,par=NULL,link='logit',modified=TRUE,cl.cores=1){

  if(!is.matrix(Gs)){
    Gs = as.matrix(Gs);
  }
  G.list = unclass(data.frame(Gs));

  cl = makeCluster(cl.cores);

  cl.cores = as.integer(cl.cores);

  #SingleGeneticEffect(G,Y,Z,S,fs,par,link,modified)
  res = parSapply(cl, G.list, SingleGeneticEffect,Y,Z,S,fs,par,link,modified);

  stopCluster(cl);

  return(res);

}
