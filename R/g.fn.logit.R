g.fn.logit <-
function(G,theta){
  res = rep(NA,length(G));
  res[G==0] = (1-theta)^2;
  res[G==1] = 2*(1-theta)*theta;
  res[G==2] = theta^2;
  return(res);
}
