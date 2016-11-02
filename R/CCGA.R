
CCGA = function(Y,G,X,S,fs,par=NULL,link='logit',modified=TRUE){

  ## fs[i] = prevalence for stratum i

  n.S = length(fs);
  if(is.matrix(X)){
    n.X = ncol(X);
  }else{
    n.X = 1;
    X = as.matrix(X);
  }
  Scount = as.vector(table(S));

  S.sort = sort(unique(S));
  if(length(S.sort)!=n.S){
    stop('Numbers of penetrences and strata do not match!');
  }
  if(any(S.sort!=1:n.S)){
    if(n.S==1){
      stop('Stratum variable should be coded as 1');
    }else{
      stop(paste('Stratum variable should be coded as', 1,'...',n.S));
    }
  }
  uniqueY = sort(unique(Y));
  if(length(uniqueY)!=2){
    stop('There should be at least one case and one control');
  }
  if(uniqueY[1]!=0 | uniqueY[2]!=1){
    stop('Cases should be coded as 1 and controls should be coded as 0');
  }

  nn0 = nn1 = rep(NA,n.S);
  for(s in 1:n.S){
    nn1[s] = sum(S==s&Y==1);
    nn0[s] = sum(S==s&Y==0);
  }

  ns = nn1 + nn0;

  lambda = -1/(1-fs)+1/fs/(1-fs)*nn1/(nn0+nn1);
  para = list(n.S=n.S,n.X=n.X,fs=fs,lambda=lambda);

  data = as.data.frame(cbind(Y,G,X));
  names(data) = c('Y','G',paste('X',1:n.X,sep=''));

  f0 = as.numeric(Scount%*%fs)/sum(Scount);
  offset = log((1-f0)/f0)+log(sum(nn1)/sum(nn0));
  fit.log = glm(Y~.,family=binomial(link='logit'),data=data,offset=rep(offset,length(Y)));
  est.log = summary(fit.log)$coef[,1];
  se.log = summary(fit.log)$coef[,2];

  para$theta = mean(G)/2;
  para$fs = fs;
  para$lambda = -1/(1-fs)+1/fs/(1-fs)*nn1/(nn0+nn1);
  if(is.null(par)){
    para$alpha = est.log[1];
    para$betaG = est.log[2];
    para$betaX = est.log[-(1:2)];
    para$betaS = 0;
  }else{
    para$alpha = par$alpha;
    para$betaG = par$betaG;
    para$betaX = par$betaX;
    para$betaS = par$betaS;
  }

  dat = list(S=S,X=as.matrix(X),G=G,Y=Y,ns=ns);

  if(link=='logit'){
    if(modified){
      res = mpMLE.logit(dat,para);
    }else{
      res = mpMLE.logit(dat,para);
      para$par.initial = res$par;
      res0 = pMLE.logit(dat,para);
    }
  }
  if(link=='probit'){
    if(modified){
      res = mpMLE.probit(dat,para);
    }else{
      res = mpMLE.probit(dat,para);
      para$par.initial = res$par;
      res0 = pMLE.probit(dat,para);
    }
  }


  var.nam.logit = c('Intercept', 'SNP', paste('Covariate',1:n.X,sep='.'));
  if(n.S>1){
    var.nam.mpMLE = c('MAF', 'Intercept', 'Stratum', 'SNP', paste('Covariate',1:n.X,sep='.'));
  }else{
    var.nam.mpMLE = c('MAF', 'Intercept', 'SNP', paste('Covariate',1:n.X,sep='.'));
  }
  LOGIT = round(data.frame(EST=as.numeric(est.log),SE=as.numeric(se.log)),5);
  mpMLE=round(data.frame(EST=res$par,SE=res$se),5);
  row.names(LOGIT) = var.nam.logit;
  row.names(mpMLE) = var.nam.mpMLE;
  if(modified){
    res = list(LOGIT=LOGIT,mpMLE=mpMLE);
  }else{
    pMLE = round(data.frame(EST=res0$par,SE=res0$se),5);
    row.names(pMLE) = c(var.nam.mpMLE,paste('Lambda',1:n.S,sep='.'));
    res = list(LOGIT=LOGIT,mpMLE=mpMLE,pMLE=pMLE);
  }

  return(res);

}
