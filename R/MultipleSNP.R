MultipleSNP = function(Gs,Y,Z,S,fs,par=NULL,link='logit',modified=TRUE,cl.cores=1){

  SingleGeneticEffect = function(G){

    res = SingleSNP(Y,G,Z,S,fs,par,link,modified)$mpMLE;
    if(length(fs)==1){
      rr = as.numeric(res[3,]);
    }else{
      rr = as.numeric(res[4,]);
    }
    stat = (rr[1]/rr[2])^2;
    p.value = ifelse(stat<20,1-pchisq(stat,1),-pchisq(stat,1,log.p=TRUE));
    return(c(rr,p.value));
  }
  if(!is.matrix(Z)){
    Z = as.matrix(Z);
  }
  if(!is.data.frame(Gs)){
    Gs = as.data.frame(Gs);
  }
  Gs = as.list(Gs);

  res = as.data.frame(t(mclapply(Gs,FUN=SingleGeneticEffect,mc.cores=cl.cores)));
  names(res) = c('log-OR','SE','p-value');

  return(as.data.frame(t(res)));

}
