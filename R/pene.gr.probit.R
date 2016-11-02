pene.gr.probit <-
function(S,G,X,Y,par){
  
  if(length(par)-ncol(X)==2){
    if(ncol(X)==1){
      tmp = par[1] + par[2]*G + as.numeric(X*par[-(1:2)]);
    }else{
      tmp = par[1] + par[2]*G + as.numeric(X%*%par[-(1:2)]);
    }
  }else{
    if(ncol(X)==1){
      tmp = par[1] + par[2]*S + par[3]*G + as.numeric(X*par[-(1:3)]);
    }else{
      tmp = par[1] + par[2]*S + par[3]*G + as.numeric(X%*%par[-(1:3)]);
    }
  }
  
  return((2*Y-1)*dnorm(tmp));
}
