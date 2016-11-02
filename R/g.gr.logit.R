g.gr.logit <-
function(G,theta){
  res = rep(NA,length(G));
  res[G==0] = 2*(theta-1);
  res[G==1] = 2-4*theta;
  res[G==2] = 2*theta;
  return(res);
}
