g.grr.probit <-
function(G,theta){
  res = rep(NA,length(G));
  res[G==0|G==2] = 2;
  res[G==1] = -4;
  return(res);
}
