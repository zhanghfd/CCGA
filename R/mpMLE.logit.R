mpMLE.logit <-
function(dat,para){ 
  
  n.S = para$n.S;
  ns = dat$ns;
  n = sum(ns);
  fs = para$fs;
  S = dat$S;
  X = dat$X;
  if(is.vector(X)){
    n.X = 1;
  }else{
    n.X = ncol(X);
  }
  G = dat$G;
  Y = dat$Y;
  if(n.S==1){
    par0 = c(para$theta,para$alpha,para$betaG,para$betaX);    
    n.p = 3+n.X;
    Z = cbind(1,G,X);
  }else{
    par0 = c(para$theta,para$alpha,para$betaS,para$betaG,para$betaX);
    n.p = 4+n.X;
    Z = cbind(1,S,G,X);
  }
  
  lambda = para$lambda;
  
  fn = function(allpar){
    theta = allpar[1];
    par = allpar[2:n.p];
    
    delta = delta.fn.logit(dat,c(allpar,lambda),fs,X,S,ns); 
    
    loglik = sum(log(pene.fn.logit(S,G,X,Y,par)*g.fn.logit(G,theta)))+sum(log(delta));
    
    return(-loglik);
  }
  
  gr = function(allpar){
    theta = allpar[1];
    theta = ifelse(theta < 1e-5,1e-5,ifelse(theta>1-1e-5,1-1e-5,theta));
    par = allpar[2:n.p];
    
    p.G = c((1-theta)^2, 2*theta*(1-theta), theta^2);
    p.G.1 = c(2*theta-2, 2-4*theta, 2*theta);
    
    h = h.fn.logit(par,p.G,X,S);
    h.d = h.gr.logit(par,p.G,p.G.1,X,S);
    
    tmp0 = tmp1 = ff = rep(NA,n);
    for(s in 1:n.S){
      tmp0[S==s] = 1+lambda[s]*(h[S==s]-fs[s]);
      tmp1[S==s] = lambda[s]/tmp0[S==s];
    }
    tmp2 = Y-pene.fn.logit(S,G,X,1,par);
    
    res = rep(NA,length(allpar));
    res[1] = sum(g.gr.logit(G,theta)/g.fn.logit(G,theta));
    res[2:n.p] = as.vector(tmp2 %*% Z);
    res[1:n.p] = res[1:n.p] - as.numeric(tmp1%*%h.d);
    
    return(-res);
  }
  
  score = function(allpar){
    theta = allpar[1];
    theta = ifelse(theta < 1e-5,1e-5,ifelse(theta>1-1e-5,1-1e-5,theta));
    par = allpar[2:n.p];
    
    p.G = c((1-theta)^2, 2*theta*(1-theta), theta^2);
    p.G.1 = c(2*theta-2, 2-4*theta, 2*theta);
    
    h = h.fn.logit(par,p.G,X,S);
    h.d = h.gr.logit(par,p.G,p.G.1,X,S);
    
    tmp0 = tmp1 = ff = rep(NA,n);
    for(s in 1:n.S){
      tmp0[S==s] = 1+lambda[s]*(h[S==s]-fs[s]);
      tmp1[S==s] = lambda[s]/tmp0[S==s];
    }
    tmp2 = Y-pene.fn.logit(S,G,X,1,par);
    
    res = matrix(NA,n,length(allpar));
    res[,1] = g.gr.logit(G,theta)/g.fn.logit(G,theta);
    res[,2:n.p] = sweep(Z,1,tmp2,'*');
    res[,1:n.p] = res[,1:n.p] - sweep(h.d,1,tmp1,'*');
    
    ress = 0;
    for(y in 0:1){
      for(s in 1:n.S){
        tmp = res[Y==y&S==s,];
        tmp = sweep(tmp,2,colMeans(tmp));
        ress = ress + t(tmp)%*%tmp;
      }
    }
    
    return(ress);
  }
  
  allpar = par0;
  
  fit = multiroot(f=gr, start = allpar);
  
  while(is.na(fit$estim.precis) | fit$estim.precis>1e-4){
    allpar[2] = par0[2]+runif(1,-3,0);
    fit = multiroot(f=gr, start = allpar);
  }
  par = fit$root;
  
  h = numericGradient(gr, par, eps=1e-6);
  infinv = solve(h);
  sigma = score(par);
  varcov = (infinv%*%sigma)%*%t(infinv);
  se = sqrt(diag(varcov));
  
  return(list(par=par,se=se));
}
