
SingleGeneticEffect <-
function(G,Y,Z,S,fs,par,link,modified){

  res = SingleSNP(Y,G,Z,S,fs,par,link,modified)$mpMLE;
  if(length(fs)==1){
    rr = as.numeric(res[3,]);
  }else{
    rr = as.numeric(res[4,]);
  }
  stat = (rr[1]/rr[2])^2;
  p.value = ifelse(stat>20,1-pchisq(stat,1),pchisq(stat,1,log.p=TRUE));
  return(c(rr,p.value));
}
