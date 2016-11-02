delta.fn.probit <-
function(dat,allpar,fs,X,S,ns){
  n = nrow(X);
  n.X = ncol(X);
  n.S = length(fs);
  n.p = length(allpar) - n.S;
  
  theta = allpar[1];
  par = allpar[2:n.p];
  lambda = allpar[-(1:n.p)];
  
  p.G = c((1-theta)^2,2*theta*(1-theta),theta^2);
  h = h.fn.probit(par,p.G,X,S);
  
  delta = rep(NA,n);
  for(s in 1:n.S){
    delta[S==s] = 1/ns[s]/(1+(h[S==s]-fs[s])*lambda[s]);
  }
  
  return(delta);
}
