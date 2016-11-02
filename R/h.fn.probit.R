h.fn.probit <-
function(par,p.G,X,S){
  
  n = nrow(X);
  h = rep(0,n);
  
  for(g in 0:2){
    h = h + pene.fn.probit(S,g,X,1,par)*p.G[g+1];
  }
  return(h);
}
