h.gr.logit <-
function(par,p.G,p.G.1,X,S){
  
  n = nrow(X);
  n.X = ncol(X);
  n.p = length(par)+1;
  onestratum = (n.p-n.X==3);
  h.d = matrix(0,n,n.p);
  
  for(g in 0:2){
    tmp0 = pene.fn.logit(S,g,X,1,par);
    tmp1 = pene.gr.logit(S,g,X,1,par);
    if(onestratum){
      h.d[,1] = h.d[,1] + tmp0*p.G.1[g+1];
      h.d[,2] = h.d[,2] + tmp1*p.G[g+1];
      h.d[,3] = h.d[,3] + tmp1*p.G[g+1]*g;
      h.d[,-(1:3)] = h.d[,-(1:3)] + sweep(X,1,tmp1*p.G[g+1],'*');
    }else{
      h.d[,1] = h.d[,1] + tmp0*p.G.1[g+1];
      h.d[,2] = h.d[,2] + tmp1*p.G[g+1];
      h.d[,3] = h.d[,3] + tmp1*p.G[g+1]*S;
      h.d[,4] = h.d[,4] + tmp1*p.G[g+1]*g;
      h.d[,-(1:4)] = h.d[,-(1:4)] + sweep(X,1,tmp1*p.G[g+1],'*');
    }
  }
  return(h.d);
}
