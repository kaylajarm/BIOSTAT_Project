library(ROCR)
#added this library because this contains predictions


makeOneTable<-function(d,k1,k2) {
  nc1<-sum(k1)
  nc2<-sum(k2)
  tab<-matrix(,nrow=2*nrow(d),ncol=max(nc1,nc2))
  tab[1:nrow(d),1:nc1]<-d[,k1]
  tab[(nrow(d)+1):nrow(tab),1:nc2]<-d[,k2]
  tab
}

quantileAdjust<-function(y,N=prod(y$lib.size)^(1/ncol(y$data)),alpha=0,null.hypothesis=F,n.iter=5,r.init=NULL,tol=.001) {
  # adjust data for estimation of phi (common dispersion only)
  if (is.null(r.init)) { 
    r.init<-1/findMaxD2(y,alpha=10)-1 
  }
  g<-unique(y$group)
  k1<-y$group==g[1]; k2<-y$group==g[2]
  y1<-matrix(y$data[,k1],ncol=sum(k1)); y2<-matrix(y$data[,k2],ncol=sum(k2))
  r<-r.init
  print(length(r))
  rprev<-r+1
  count<-0
  while( count < n.iter ) { 
    count<-count+1
    cat("[quantileAdjust] Iteration:",count,"\n")
    rprev<-r
    ps<-estimatePs(y1,y2,y$lib.size[k1],y$lib.size[k2],r)
    if (null.hypothesis==T) {
      p1<-pnbinom(y1-1,size=r,mu=outer(ps$p,y$lib.size[k1]))+dnbinom(y1,size=r,mu=outer(ps$p,y$lib.size[k1]))/2
      p2<-pnbinom(y2-1,size=r,mu=outer(ps$p,y$lib.size[k2]))+dnbinom(y2,size=r,mu=outer(ps$p,y$lib.size[k2]))/2
      mu1<-outer(ps$p,rep(N,sum(k1)))
      mu2<-outer(ps$p,rep(N,sum(k2)))
    } else {
      p1<-pnbinom(y1-1,size=r,mu=outer(ps$p1,y$lib.size[k1]))+dnbinom(y1,size=r,mu=outer(ps$p1,y$lib.size[k1]))/2
      p2<-pnbinom(y2-1,size=r,mu=outer(ps$p2,y$lib.size[k2]))+dnbinom(y2,size=r,mu=outer(ps$p2,y$lib.size[k2]))/2
      mu1<-outer(ps$p1,rep(N,sum(k1)))
      mu2<-outer(ps$p2,rep(N,sum(k2)))
    }
    pseudo<-interpolateHelper(cbind(mu1,mu2),cbind(p1,p2),r)
    pseudo[pseudo<0]<-0  # values less than zero for small samples seems to make this unstable
    d<-list(group=y$group,data=pseudo)
    r<-1/findMaxD2(d,alpha=alpha)-1
    if( max(abs(rprev-r)) < tol ) { break }
  }
  return(list(r=r,pseudo=pseudo,p=cbind(p1,p2),mu=cbind(mu1,mu2),ps=ps,N=N))
}

interpolateHelper<-function(mu,p,r) {
   if(length(r)==1) { r<-matrix(r,nrow=nrow(mu),ncol=ncol(mu)) }
   if(length(r)==nrow(mu)) { r<-outer(r,rep(1,ncol(mu))) }
   mx<-max(qnbinom(p,size=r,mu=mu))+1
   psudo<-matrix(,nrow=nrow(p),ncol=ncol(p))
   for (i in 1:nrow(psudo)) {
     for (j in 1:ncol(psudo)) {
       a<-pnbinom(-1:(mx-1),size=r[i,j],mu=mu[i,j])
       psudo[i,j]<-approx(c(0,a+c(diff(a)/2,0)),c(-.5,0:mx),xout=p[i,j])$y
     }
   }
   return(psudo)
}

interpolateHelper1<-function(mu,p,r) {
   if(length(r)==1) { r<-rep(r,length(mu)) }
   mx<-max(qnbinom(p,size=r,mu=mu))+1
   psudo<-matrix(,nrow=nrow(p),ncol=ncol(p))
   for (i in 1:nrow(psudo)) {
       #a<-pnbinom(-1:(mx-1),size=r[i],mu=mu[i,])
       a<-pnbinom(-1:(mx-1),size=r[i],mu=mu[i])
       psudo[i,]<-approx(c(0,a+c(diff(a)/2,0)),c(-.5,0:mx),xout=p[i,])$y
   }
   return(psudo)
}

quantileAdjust1<-function(y,lib.size,tol=.001,maxit=10) {
  # adjust data for estimation of phi
  # for single sample
  if (is.vector(y)) { y<-matrix(y,nrow=1) }
  delta<-findMaxD(y,alpha=100)[1]
  dprev<-delta+tol+.1
  count<-0
  N<-prod(lib.size)^(1/length(lib.size))
  while( abs(dprev-delta) > tol ) {
    count<-count+1
    dprev<-delta
    r<-1/delta-1
    ps<-estimatePs1(y,lib.size,r)
    p<-pnbinom(y-1,size=r,mu=outer(ps$p,lib.size))+dnbinom(y,size=r,mu=outer(ps$p,lib.size))/2
    pseudo<-interpolateHelper1(ps$p*N,p,r)
    delta<-findMaxD(pseudo,alpha=100)[1]
    if (count > maxit) { return(delta) }
  }
  return(delta)
}

weightedLikelihood1<-function(y,lib.size,tol=.001,alpha=0,maxit=10) {
  # adjust data for estimation of phi
  # for single sample
  if (is.vector(y)) { y<-matrix(y,nrow=1) }
  delta<-findMaxD(y,alpha=alpha)
  dprev<-delta+tol+.1
  count<-0
  N<-prod(lib.size)^(1/length(lib.size))
  while( max(abs(dprev-delta)) > tol ) {
    print(max(abs(dprev-delta)))
    count<-count+1
    dprev<-delta
    r<-1/delta-1
    ps<-estimatePs1(y,lib.size,r)
    p<-pnbinom(y-1,size=r,mu=outer(ps$p,lib.size))+dnbinom(y,size=r,mu=outer(ps$p,lib.size))/2
    pseudo<-interpolateHelper(ps$p*N,p,r)
    delta<-findMaxD(pseudo,alpha=alpha)
    if (count > maxit) { return(delta) }
  }
  return(delta)
}


ql<-function(delta,y,mu,df) {
  r<-1/delta-1
  phi<-delta/(1-delta)
  resids<-2*(y * log(pmax(1,y)/mu) - (y+r)*log((y+r)/(mu+r)))
  (sum(resids)-df)^2
}

pl<-function(delta,y,mu,df) {
  phi<-delta/(1-delta)
  abs(sum( (y-mu)^2 / (mu*(1+mu*phi))  ) - df)
}

qlo<-function(delta,y,df,lib.size) {
  r<-1/delta-1
  phi<-1/r
  ps<-estimatePs1(y,lib.size,r)
  mu<-outer(ps$p,lib.size)
  resids<-2*(y * log(pmax(1,y)/mu) - (y+r)*log((y+r)/(mu+r)))
  (sum(resids)-df)^2
}

plo<-function(delta,y,df,lib.size) {
  r<-1/delta-1
  phi<-1/r
  ps<-estimatePs1(y,lib.size,r)
  mu<-outer(ps$p,lib.size)
  (sum( (y-mu)^2 / (mu*(1+mu*phi)) )  - df)^2
}

plo2<-function(delta,y,null.hypothesis=T) {
  r<-1/delta-1
  phi<-1/r
  k1<-y$group==1; k2<-y$group==2
  ps<-estimatePs(y$data[,k1],y$data[,k2],y$lib.size[k1],y$lib.size[k2],r)
  if (null.hypothesis) {
    mu<-cbind(outer(ps$p,y$lib.size[k1]),outer(ps$p,y$lib.size[k2]))
    df<-nrow(y$data)*(sum(k1)+sum(k2)-1)
  } else {
    mu<-cbind(outer(ps$p1,y$lib.size[k1]),outer(ps$p2,y$lib.size[k2]))
    df<-nrow(y$data)*(sum(k1)+sum(k2)-2)
  }
  (sum( (y$data-mu)^2 / (mu*(1+mu*phi)) )  - df)^2
}

qlo2<-function(delta,y,null.hypothesis=T) {
  r<-1/delta-1
  phi<-1/r
  k1<-y$group==1; k2<-y$group==2
  ps<-estimatePs(y$data[,k1],y$data[,k2],y$lib.size[k1],y$lib.size[k2],r)
  if (null.hypothesis) {
    mu<-cbind(outer(ps$p,y$lib.size[k1]),outer(ps$p,y$lib.size[k2]))
    df<-nrow(y$data)*(sum(k1)+sum(k2)-1)
  } else {
    mu<-cbind(outer(ps$p1,y$lib.size[k1]),outer(ps$p2,y$lib.size[k2]))
    df<-nrow(y$data)*(sum(k1)+sum(k2)-2)
  }
  resids<-2*(y$data * log(pmax(1,y$data)/mu) - (y$data+r)*log((y$data+r)/(mu+r)))
  (sum(resids)-df)^2
}

approx.expected.info<-function(x,d,qA,robust=F) {
  g<-unique(x$group); k1<-x$group==g[1]; k2<-x$group==g[2]
  obs.inf<-condLogLikDer2(qA$pseudo[,k1],d,der=2,doSum=0)*(-1)+condLogLikDer2(qA$pseudo[,k2],d,der=2,doSum=0)*(-1)
  t<-rowSums(qA$pseudo)
  if (robust) {
    require(MASS)
    inf.lm<-rlm(obs.inf~-1+t)
  } else {
    inf.lm<-lm(obs.inf~-1+t)
  }
  exp.inf.approx<-fitted(inf.lm)
  exp.inf.approx
}

alpha.approxeb<-function(x) {
  g<-unique(x$group)
  k1<-x$group==g[1]; k2<-x$group==g[2]
  qA<-quantileAdjust(x,alpha=10,null.hypothesis=T)  # alpha large to make common estimator
  d<-1/(1+qA$r[1])  # common delta
  scores<-condLogLikDer2(qA$pseudo[,k1],d,der=1,doSum=0)+condLogLikDer2(qA$pseudo[,k2],d,der=1,doSum=0)
  exp.inf<-approx.expected.info(x,d,qA)
  sigma2.0.est<-optimize(tau2.0.objective,c(0,500),info.g=exp.inf,score.g=scores)$min
  alpha<-1/(sigma2.0.est*nrow(x$data)*mean(exp.inf))
  list(sigma2.0.est=sigma2.0.est,alpha=alpha,scores=scores,infos=exp.inf,quantileAdjusted=qA)
}

tau2.0.objective<-function(tau2.0,info.g,score.g) {
  G<-length(info.g)
  denom<-info.g*(1+info.g*tau2.0)
  ( mean(score.g^2/denom)-1 )^2
  #m.info<-mean(info.g)
  #denom<-(1+info.g*tau2.0)
  #( mean(score.g^2/denom)-m.info )^2
}

phiScoreTest<-function(y,mu,r) {
  dldr<-rowSums(digamma(y+r)-digamma(r)+1+log(r)-log(mu+r)-(y+r)/(mu+r))
  d2ldrdmu<-rowSums((y+r)/(mu+r)^2-1/(mu+r))
  d2ldr2<-rowSums(trigamma(y+r)-trigamma(r)+1/r-1/(mu+r)+(y-mu)/(mu+r)^2)
  d2ldmu2<-rowSums((y+r)/(mu+r)^2-y/mu^2)
  score<-dldr/sqrt(d2ldrdmu^2/d2ldmu2-d2ldr2)
  score
}

dldr<-function(y,mu,r) {
  rowSums(digamma(y+r)-digamma(r)+1+log(r)-log(mu+r)-(y+r)/(mu+r))
}

runScoreTestComparison<-function(samples=100,tags=100,n=2,lib.size.range=c(20000,80000),lib.size=NULL,lib.p=.0001,phi=0.25,b.mult=2,mu=NULL) {
  # this runs a comparison for common dispersion, many tags, computes score statistic
  ests2<-matrix(,nrow=samples,ncol=6)
  d<-NULL
  pseudo<-NULL
  mus<-NULL
  for(i in 1:samples) {
    if (is.null(lib.size)) { lib.size<-round(runif(2*n,lib.size.range[1],lib.size.range[2])) }
    x<-makeDataset(lib.size=lib.size,n1=n,n2=n,p.a=lib.p,b.mult=b.mult,tags=tags,true=tags/10,disp=phi,mu=mu)
    d[[i]]<-x
    f<-quantileAdjust(x,null.hypoth=T,alpha=20)
    #f<-quantileAdjust(x,null.hypoth=T)
    pseudo[[i]]<-f$pseudo
    mus[[i]]<-f$mu
    true_phi_values <- rep(phi, samples)  # Set true_phi as a vector with the value 'phi' for each sample
    print(i)
    #next
    #ests2[i,1]<-optimize(fullLogLik.delta,c(0,1),max=T,y=x$data,lib.size=x$lib.size,doSum=1)$max
    #above causes an error, bottom fixes the max = TURE issue
    ests2[i,1] <- optimize(fullLogLik.delta, c(0,1), maximum=TRUE, y=x$data, lib.size=x$lib.size, doSum=1)$maximum
    ests2[i,2]<-findMaxD2(x,alpha=20)[1]
    ests2[i,3]<-optimize(plo2,c(0,1),y=x,null.hypoth=T)$min
    ests2[i,4]<-optimize(qlo2,c(0,1),y=x,null.hypoth=T)$min
    #ests2[i,5]<-optimize(adjProfLik,c(0,1),x=x,null.hypoth=T,max=T)$max
    #above causes an error, bottom fixes the max = TURE issue
    ests2[i,5]<-optimize(adjProfLik, c(0,1), x=x, null.hypoth=TRUE, maximum=TRUE)$maximum
    ests2[i,6]<-1/(1+f$r[1])
  }
  return(list(d=d, est=ests2, pseudo=pseudo, mus=mus, true_phi=true_phi_values))
  #return(list(d=d,est=ests2,pseudo=pseudo,mus=mus))
  #return(list(d=d,est=ests2,lr.pl=lr.pl,score.pl=score.pl,lr.qcml=lr.qcml,score.qcml=score.qcml))
}

runPhiEstimateComparison<-function(samples=100,tags=10,n=5,lib.size.range=c(20000,80000),lib.p=.0001,phi=0.42) {
  # this runs a comparison for common dispersion, many tags
  ests2<-matrix(,nrow=samples,ncol=6)
  x<-matrix(,ncol=n,nrow=tags)
  d<-NULL
  for(i in 1:samples) {
    lib.size<-round(runif(n,lib.size.range[1],lib.size.range[2]))
    for (j in 1:tags) { x[j,]<-myrnbinom(n,mu=lib.size*lib.p,phi=phi) }
    y<-list(data=x,lib.size=lib.size)
    d[[i]]<-y
    #ests2[i,1]<-optimize(fullLogLik.delta,c(0,1),max=T,y=x,lib.size=lib.size,doSum=1)$max
    ests2[i,1] <- optimize(fullLogLik.delta, c(0,1), maximum=TRUE, y=x, lib.size=lib.size, doSum=1)$maximum
    ests2[i,2]<-findMaxD(x,alpha=100)[1]
    ests2[i,3]<-optimize(plo,c(0,1),y=x,lib.size=lib.size,df=tags*(n-1))$min
    ests2[i,4]<-optimize(qlo,c(0,1),y=x,lib.size=lib.size,df=tags*(n-1))$min
    #ests2[i,5]<-optimize(f=adjProfLikM,interval=c(0,1),y=x,lib.size=lib.size,max=T)$max
    ests2[i,5] <- optimize(f=adjProfLikM, interval=c(0,1), y=x, lib.size=lib.size, maximum=TRUE)$maximum
    ests2[i,6]<-quantileAdjust1(x,lib.size)
    print(i)
  }
  return(list(p=lib.p,data=d,est=ests2,phi=phi,n=n,samples=samples,lib.size.range=lib.size.range))
}


runPhiEstimateComparison1o<-function(samples=100,n=5,lib.size.range=c(20000,80000),lib.p=.0001,phi=0.42) {
  # this runs a comparison for single samples w/ offsets
  ests2<-matrix(,nrow=samples,ncol=6)
  x<-matrix(,ncol=n,nrow=samples)
  d<-NULL
  for(i in 1:samples) {
    lib.size<-round(runif(n,lib.size.range[1],lib.size.range[2]))
    x[i,]<-myrnbinom(n,mu=lib.size*lib.p,phi=phi)
    y<-list(data=x[i,],lib.size=lib.size)
    d[[i]]<-lib.size
    ests2[i,1]<-optim(c(lib.p,phi/(phi+1)),fullLogLik1o,control=list(fnscale=-1),y=x[i,],lib.size=lib.size)$par[2]  # straight ML
    ests2[i,2]<-findMaxD(matrix(x[i,],nrow=1),alpha=100)
    ests2[i,3]<-optimize(plo,c(0,1),y=x[i,],lib.size=lib.size,df=n-1)$min
    ests2[i,4]<-optimize(qlo,c(0,1),y=x[i,],lib.size=lib.size,df=n-1)$min
    #ests2[i,5]<-optimize(f=adjProfLik1o,interval=c(0,1),y=x[i,],lib.size=lib.size,max=T)$max
    ests2[i,5] <- optimize(f = adjProfLik1o, interval = c(0,1), y = x[i,], lib.size = lib.size, maximum = TRUE)$maximum
    ests2[i,6]<-quantileAdjust1(x[i,],lib.size)
    print(i)
  }
  return(list(p=lib.p,data=x,lib.sizes=d,est=ests2,phi=phi,n=n,samples=samples,lib.size.range=lib.size.range))
}

runPhiEstimateComparison1<-function(samples=100,n=5,mu=4,phi=0.42) {
  ests2<-matrix(,nrow=samples,ncol=5)
  x<-matrix(myrnbinom(n*samples,mu=mu,phi=phi),ncol=n)
  mus<-rowMeans(x)
  for(i in 1:samples) {
    ests2[i,1]<-optim(c(mu,.5),fullLogLik1,control=list(fnscale=-1),y=x[i,])$par[2]  # straight ML
    ests2[i,3]<-optimize(pl,c(0,1),y=x[i,],mu=mus[i],df=n-1)$min  # pseudo
    ests2[i,4]<-optimize(ql,c(0,1),y=x[i,],mu=mus[i],df=n-1)$min  # quasi
    ests2[i,5]<-optimize(adjProfLik1,c(1e-15,1-1e-15),maximum=T,y=x[i,])$max # Cox-Reid adjusted profile likelihood
    print(i)
  }
  ests2[,2]<-findMaxD(x,alpha=100)
  return(list(mu=mu,data=x,est=ests2,phi=phi,n=n,samples=samples))
}

plotPhiEstimateComparison<-function(est,ylim=NULL) {
  if (is.null(ylim)) { ylim<-c(min(est$est),max(est$est)) }
  main<-paste("n=", est$n,", phi=",est$phi,sep="")
  if (!is.na(match("p",names(est)))) { main<-paste(main,", lambda=",est$p,sep="") }
  if (!is.null(est$mu)) { main<-paste(main,", mu=",est$mu,sep="") }
  if (!is.null(est$lib.range)) { main<-paste(main,", lib.range=",est$lib.size.range[1],"-",est$lib.size.range[2],sep="") }
  if (ncol(est$est)==5) { labels<-c("ML","CML","PL","QL","CR") } else { labels<-c("ML","CML","PL","QL","CR","qCML") }
  boxplot(est$est~col(est$est),names=labels,cex.axis=1,cex.main=1,ylim=ylim,main=main)
  abline(h=est$phi/(est$phi+1))
}

#calcAUCs function needs to be semi redone because the wald calculations are causing errors with the prediction thing from the ROCR package
#easier to call wald first i believe and redo it like that 
calcAUCs <- function(s, mu) {
  aa <- fisherExactScore(s, mu)
  aucs <- matrix(, nrow(s$est), ncol(s$est) * 4)
  for (j in 1:ncol(s$est)) {
    col <- j * 4 - 3
    for (i in 1:length(aa)) {
      #wald used to see how well parameters fit
      wald_result <- abs(waldTest(s$d[[i]], 1/s$est[i,j] - 1)) #added wald result to do wald test
      #changes down
      if (any(is.na(wald_result))) next
      truth_values <- s$d[[i]]$truth #truth is also from ROCR package i think, had issues with that as well
      if (any(is.na(truth_values))) next
      ap <- prediction(wald_result, truth_values) #uses wald result and truth values instead of the predictions and other stuff i couldnt tell you what it means
      aucs[i, col] <- performance(ap, "auc")@y.values[[1]] #same, no changes
      
      #add different tests that are supposed to be derived from the package, but i thought it was easier to just add them here
      #i think these should work to remove the NAs the code was having errors with
      #compares models, likelihood ration
      lr_result <- abs(lrTest(s$d[[i]], 1/s$est[i,j] - 1))
      if (any(is.na(lr_result))) next
      ap <- prediction(lr_result, truth_values)
      aucs[i, col + 1] <- performance(ap, "auc")@y.values[[1]]
      #null hypoth test for score functions
      score_result <- abs(scoreTest(s$d[[i]], 1/s$est[i,j] - 1))
      if (any(is.na(score_result))) next
      ap <- prediction(score_result, truth_values)
      aucs[i, col + 2] <- performance(ap, "auc")@y.values[[1]]
      #associations bt values
      fisher_result <- -log10(aa[[i]][, j])
      if (any(is.na(fisher_result))) next
      ap <- prediction(fisher_result, truth_values)
      aucs[i, col + 3] <- performance(ap, "auc")@y.values[[1]]
    }
  }
  aucs
}
#used different tests for more accuracy in removing NAs and get a good AUC score



calcRankDE<-function(s,mu=NULL) {
  aa<-fisherExactScore(s,mu)
  rnks<-matrix(,nrow(s$est),ncol(s$est)*4)
  for (j in 1:ncol(s$est)) {
    print(j)
    col<-j*4-3
    for (i in 1:length(aa)) {
      rnks[i,col]<-mean(rank(abs(waldTest(s$d[[i]],1/s$est[i,j]-1)))[1:100],na.rm=T)
      rnks[i,col+1]<-mean(rank(abs(lrTest(s$d[[i]],1/s$est[i,j]-1)))[1:100],na.rm=T)
      rnks[i,col+2]<-mean(rank(abs(scoreTest(s$d[[i]],1/s$est[i,j]-1)))[1:100],na.rm=T)
      rnks[i,col+3]<-mean(rank(-log10(aa[[i]][,j]))[1:100],na.rm=T)
    }
  }
  rnks
}

calcRankDiff<-function(s,mu=NULL) {
  aa<-fisherExactScore(s,mu)
  rnks<-matrix(,nrow(s$est),ncol(s$est))
  for (i in 1:length(aa)) {
    a<-mean(rank(-log10(aa[[i]][,6]))[1:100],na.rm=T)
    for (j in 1:ncol(s$est)) {
      rnks[i,j]<-a-mean(rank(-log10(aa[[i]][,j]))[1:100],na.rm=T)
    }
  }
  rnks
}


calcNominal<-function(s,mu=NULL,size=.05) {
  nom<-matrix(,nrow(s$est),ncol(s$est)*4)
  for (j in 1:ncol(s$est)) {
    col<-j*4-3
    for (i in 1:nrow(s$est)) {
      wa<-abs(waldTest(s$d[[i]],1/s$est[i,j]-1))
      nom[i,col]<-sum(pnorm(wa,lower=F)<size,na.rm=T)
      sc<-abs(scoreTest(s$d[[i]],1/s$est[i,j]-1))
      nom[i,col+1]<-sum(pnorm(sc,lower=F)<size,na.rm=T)
      lr<-lrTest(s$d[[i]],1/s$est[i,j]-1)
      nom[i,col+2]<-sum(pchisq(lr,1,lower=F)<size,na.rm=T)
      if( is.null(mu) ) { mus<-s$mu[[i]][,1] } else { mus<-mu }
      pseudo<-s$pseudo[[i]]
      r<-1/s$est[i,j]-1
      nom[i,col+3]<-sum(fisherTest(pseudo,s$d[[i]]$group,mus,rep(r,nrow(pseudo)))<size,na.rm=T)
      #nom[i,col+3]<-sum(fisherTest(s$d[[i]]$data,s$d[[i]]$group,mus,1/s$est[,j]-1)<size,na.rm=T)
    }
  }
  nom<-nom/nrow(s$d[[1]]$data)
  nom
}

fisherNominal<-function(aa,size=.05) {
  res<-matrix(,nrow=length(aa),ncol=ncol(aa[[1]]))
  for (i in 1:length(aa)) {
    res[i,]<-colSums(aa[[i]]<size,na.rm=T)
  }
  res/nrow(aa[[1]])
}

calcFDRs<-function(s,aa=NULL) {
  if (is.null(aa)) { aa<-fisherExactScore(s) }
  fdrs<-matrix(,30,6*4)
  for (j in 1:ncol(s$est)) {
    col<-j*4-3
    for (i in 1:length(aa)) {
      fdrs[i,col]<-doFDR(abs(waldTest(s$d[[i]],1/s$est[i,j]-1)),s$d[[1]]$truth)$falsedis[100]
      fdrs[i,col+1]<-doFDR(abs(lrTest(s$d[[i]],1/s$est[i,j]-1)),s$d[[1]]$truth)$falsedis[100]
      fdrs[i,col+2]<-doFDR(abs(scoreTest(s$d[[i]],1/s$est[i,j]-1)),s$d[[1]]$truth)$falsedis[100]
      fdrs[i,col+3]<-doFDR(-log10(aa[[i]][,j]),s$d[[1]]$truth)$falsedis[100]
    }
  }
  fdrs
}

phiPrecision<-function(run) {
  nr<-nrow(run$est)
  nc<-ncol(run$est)
  vs.b<-matrix(,nrow=nr,ncol=nc)
  mu.b<-matrix(,nrow=nr,ncol=nc)
  vs.s<-matrix(,nrow=nr,ncol=nc)
  phis<-run$est/(1-run$est)
  for(i in 1:nr) {
     data<-run$data[[i]]$data
     lib.size<-run$data[[i]]$lib.size
     mu.true<-outer(rep(run$p,nrow(data)),lib.size)
     v.true<-log(1/(mu.true*(1+mu.true*run$phi)))
     for(j in 1:nc) {
       r<-1/phis[i,j]
       ps<-estimatePs1(data,lib.size,r)
       mu.hat<-outer(ps$p,lib.size)
       v.hat<-log(1/(mu.hat*(1+mu.hat/r)))
       vs.b[i,j]<-mean(v.hat-v.true)
       mu.b[i,j]<-mean(mu.hat-mu.true)
       vs.s[i,j]<-mean((v.hat-v.true)^2)
     }
     print(i)
  }
  return(list(bias=vs.b,bias2=vs.s,bias.mu=mu.b))
}

plotPhiPrecision<-function(x,est,type=1,ylim=c(0,1)) {
  if (type==1) { 
    bias<-x$bias 
    label<-paste("(Bias in log-Precision) n=",est$n," phi=",est$phi,sep="")
  } else {
    bias<-x$bias2 
    label<-paste("(Squared Bias in log-Precision) n=",est$n," phi=",est$phi,sep="")
  }
  if (!is.null(est$mu)) { label<-paste(label," mu=",est$mu,sep="") } else { label<-paste(label," lambda=",est$p,sep="") }
  if (!is.null(est$data[[1]]$data)) { label<-paste(label," tags=",nrow(est$data[[1]]$data),sep="") }
  boxplot(bias~col(bias),names=c("ML", "CML", "PL", "QL", "CR", "qCML"),main=label,ylim=ylim)
  abline(h=0,col="blue")
}

plotPhiEstimateComparison1o<-function(est,ylim=NULL) {
  if (is.null(ylim)) {
    boxplot(est$est~col(est$est),names=c("ML","CML","PL","QL","CR","qCML"),cex.main=0.8,cex.axis=.8,
            main=paste("n=", est$n,", phi=",est$phi,", lambda=",est$p,sep=""))
  } else {
    boxplot(est$est~col(est$est),names=c("ML","CML","PL","QL","CR","qCML"),cex.main=0.8,cex.axis=.8,cex.main=0.8,ylim=ylim,cex.axis=.8,
            main=paste("n=", est$n,", phi=",est$phi,", lambda=",est$p,sep=""))
  }
  abline(h=est$phi/(est$phi+1),col="blue")
}

summaryPhiComparison<-function(est) {
  label<-paste("(Summary for estimation: n=",est$n," phi=",est$phi,sep="")
  if (!is.null(est$mu)) { label<-paste(label," mu=",est$mu,sep="") } else { label<-paste(label," lambda=",est$p,sep="") }
  if (!is.null(est$data[[1]]$data)) { label<-paste(label," tags=",nrow(est$data[[1]]$data),sep="") }
  #if (!is.null(est$lib.range)) { main<-paste(main,", lib.range=",est$lib.size.range[1],"-",est$lib.size.range[2],sep="") }
  label<-paste(label,")",sep="")
  cat(label,"\n")
  summ<-matrix(,nrow=2,ncol=ncol(est$est))
  if (ncol(summ)==5) {
    colnames(summ)<-c("ML","CML","PL","QL","CR")
  } else {
    colnames(summ)<-c("ML","CML","PL","QL","CR","qCML")
  }
  rownames(summ)<-c("Avg Bias","MSE")
  delta<-est$phi/(1+est$phi)
  summ[1,]<-colMeans((est$est-delta))
  summ[2,]<-colMeans((est$est-delta)^2)
  return(summ)
}


plotComparison1<-function(s,truth,rocx=1,labels=c("CWL.exact","Lu"),cols=c("black","red","blue"),ltys=c(1,3),x.legend=0,y.legend=1,...) {
  require(ROCR)
  for(i in 1:ncol(s)) {
    ap<-prediction(s[,i],truth)
    ap<-performance(ap,measure="tpr", x.measure="fpr")
    if(i==1) {
      plot(ap@x.values[[1]],ap@y.values[[1]],type="l",lty=ltys[i],xlab="FPR",ylab="TPR",xlim=c(0,rocx),col=cols[i],...)
    } else {
      lines(ap@x.values[[1]],ap@y.values[[1]],lty=ltys[i],col=cols[i],...)
    }
  }
  legend(x.legend,y.legend,labels,lty=ltys,cex=1.5,col=cols,lwd=4)
}

plotComparison2<-function(s,x,rocx=1) {
  require(ROCR)
  col<-c("black","black","black","blue")
  for (i in 1:ncol(s)) { 
    ap<-prediction(abs(s[,i]),x$truth)
    ap<-performance(ap,measure="tpr", x.measure="fpr")
    if (i==1) {
      plot(ap@x.values[[1]],ap@y.values[[1]],type="l",lty=i,xlab="FPR",ylab="TPR",col=col[i],xlim=c(0,rocx))
    } else {  
      lines(ap@x.values[[1]],ap@y.values[[1]],col=col[i],lty=i)
    }
  }
  legend(0,1,1:ncol(s),col=col,lty=1:4)
}

plotComparison<-function(s,truth,labels=c("CWL.exact","Lu"),ltys=c(1,3),cols=c("black","red","blue"),x.legend=0,y.legend=1,...) {
  require(ROCR)
  for(i in 1:ncol(s)) {
    ap<-doFDR(s[,i],truth)
    if(i==1) {
      plot(ap$selected,ap$falsedis,type="l",lty=ltys[i],xlab="Tags Selected",ylab="False Discoveries",col=cols[i],...)
    } else {
      lines(ap$selected,ap$falsedis,lty=ltys[i],col=cols[i],...)
    }
  }
  legend(x.legend,y.legend,labels,lty=ltys,cex=1.5,col=cols,lwd=4)
}


outputComparison<-function(x,n=1000,f="this_simulation.xls",wcml.col=1) {
  v1<-order(-x$us.lr[,wcml.col])
  v2<-order(-abs(x$us.wald[,wcml.col]))
  v3<-order(-abs(x$them.wald))
  v4<-order(-abs(x$them.lr))
  rw<-union(union(v1[1:n],v2[1:n]),union(v3[1:n],v4[1:n]))
  cl<-order(x$x$group)
  y<-x$x$data[rw,cl]
  z<-cbind(y,us.lr=x$us.lr[rw,wcml.col],us.lr.rank=rank(-x$us.lr[,wcml.col])[rw],
             us.wald=x$us.wald[rw,wcml.col],us.wald.rank=rank(-abs(x$us.wald[,wcml.col]))[rw],
             us.score=x$us.score[rw,wcml.col],us.score.rank=rank(-abs(x$us.score[,wcml.col]))[rw],
             us.phi=1/x$r[rw,wcml.col],
             lu.wald=x$them.wald[rw],lu.wald.rank=rank(-abs(x$them.wald))[rw],
             lu.lr=x$them.lr[rw],lu.lr.rank=rank(-x$them.lr)[rw],
             lu.score=x$them.score[rw],lu.score.rank=rank(-x$them.score)[rw],
             lu.phi=x$phi[rw],truth=x$x$truth[rw],disp=x$x$disp[rw])
  write.table(z,f,sep="\t",row.names=F)
}
  
condLogLikDer2<-function(y,delta,grid=T,der=1,doSum=1) {
  # delta is 1/(1+r) below ... 1/delta-1=r
  # der is derivative (0th deriv is the function)
  # if grid=T, calculates the likelihood (derivatives) at a grid of deltas
  # if grid=F, length(delta)=nrow(y) and a different delta is used for each row
  if (is.vector(y)) {
    n<-length(y)
    t<-sum(y)
    g<-1    
    y<-matrix(y,nrow=1)
  } else {
    t<-rowSums(y,na.rm=TRUE)
    n<-rowSums(!is.na(y))
    g<-dim(y)[1]
  }
  logliks<-matrix(rep(0,length(delta)*g),nrow=g)
  if (grid==T) {
    for (i in seq(along=delta)) {
      d<-delta[i]
      r<-(1/d)-1
      if (der==1) {
        ll<-condLogLikDerOld(y,r,der=1)*(-d^(-2))
      } else if(der==2) {
        ll<-condLogLikDerOld(y,r,der=1)*2*(d^(-3))+condLogLikDerOld(y,r,der=2)*(d^(-4))
      } else if(der==0) {
        ll<-condLogLikDerOld(y,r,der=0)
      }
      logliks[,i]<-ll
    }
  } else {
      r<-(1/d)-1
      if (der==1) {
        ll<-condLogLikDerOld(y,r,der=1)*(-d^(-2))
      } else if(der==2) {
        ll<-condLogLikDerOld(y,r,der=1)*2*(d^(-3))+condLogLikDerOld(y,r,der=2)*(d^(-4))
      } else if(der==0) {
        ll<-condLogLikDerOld(y,r,der=0)
      }
      logliks<-ll
  }
  if (doSum) {
    return(colSums(logliks))
  } else {
    return(logliks)
  }
}

iidCondLogLik<-function(r,y) {
   return(condLogLikDerOld(y,r,der=0))
}

condLogLikDerOld<-function(y,r,der=1) {
  # der is derivative (0th deriv is the function)
  # this is with the parameterization with r = 1/phi
  if (is.vector(y)) {
    n<-length(y)
    t<-sum(y)
    g<-1
    y<-matrix(y,nrow=1)
  } else {
    t<-rowSums(y,na.rm=TRUE)
    n<-rowSums(!is.na(y))
    g<-dim(y)[1]
  }
##  changed may 28th
##  for (i in seq(along=r)) {
##    if (der==1) {
##      #ll<-apply(digamma(y+r[i]),1,nona.sum) + n*digamma(n*r[i]) - n*digamma(t+n*r[i]) - n*digamma(r[i])
##      ll<-rowSums(digamma(y+r[i])) + n*digamma(n*r[i]) - n*digamma(t+n*r[i]) - n*digamma(r[i])
##    } else if(der==2) {
##      ll<-rowSums(trigamma(y+r[i])) + n^2*trigamma(n*r[i]) - n^2*trigamma(t+n*r[i]) - n*trigamma(r[i])
##    } else if(der==0) {
##      ll<-rowSums(lgamma(y+r[i])) + lgamma(n*r[i]) - lgamma(t+n*r[i]) - n*lgamma(r[i])
##    }
##    logliks[,i]<-ll
##  }
  if (der==1) {
    ll<-rowSums(digamma(y+r)) + n*digamma(n*r) - n*digamma(t+n*r) - n*digamma(r)
  } else if(der==2) {
    ll<-rowSums(trigamma(y+r)) + n^2*trigamma(n*r) - n^2*trigamma(t+n*r) - n*trigamma(r)
  } else if(der==0) {
    ll<-rowSums(lgamma(y+r)) + lgamma(n*r) - lgamma(t+n*r) - n*lgamma(r)
  }
  ll
}

fullLogLik<-function(y1,y2,mu1,mu2,r) {
  # log likelihood -- 2 sample case
  if (is.vector(y1) | is.vector(y2)) {
    y1<-matrix(y1,nrow=1)
    y2<-matrix(y2,nrow=1)
  } 
  n1<-dim(y1)[2]
  n2<-dim(y2)[2]
  y<-cbind(y1,y2)
  g<-nrow(y)
  if(length(mu1)==n1 & length(mu2)==n2) {
     mu<-cbind(matrix(rep(mu1,g),nrow=g),matrix(rep(mu2,g),nrow=g))
  } else if( (sum(dim(mu1)==dim(y1))==2) & (sum(dim(mu2)==dim(y2))==2) )  {
     mu<-cbind(mu1,mu2)
  } else {
     mu<-cbind(matrix(rep(mu1,n1),nrow=g),matrix(rep(mu2,n2),nrow=g))
  }
  mu[mu==0]<-1e-12
  if (length(r)==nrow(y1)) {
     r<-matrix(r,ncol=1) %*% rep(1,n1+n2)
  } 
  ll<-rowSums(lgamma(y+r)-(r+y)*log(mu+r)-lgamma(y+1)+y*log(mu)+r*log(r)-lgamma(r))
  return(ll)
}

fullLogLik.delta<-function(delta,y,lib.size,mu=NULL,doSum=F) {
  # log likelihood -- many samples case
  if (is.vector(y)) { y<-matrix(y,nrow=1) }
  n<-ncol(y)
  g<-nrow(y)
  r<-1/delta-1
  if (is.null(mu)) { 
    ps<-estimatePs1(y,lib.size,r)
    mu<-outer(ps$p,lib.size) 
  }
  mu[mu==0]<-1e-12
  ll<-rowSums(lgamma(y+r)-(r+y)*log(mu+r)-lgamma(y+1)+y*log(mu)+r*log(r)-lgamma(r))
  if (doSum) { return(sum(ll)) } else { return(ll) }
}

lrTest<-function(x,r,p=NULL) {
  g<-unique(x$group)
  y1<-x$data[,x$group==g[1]]
  y2<-x$data[,x$group==g[2]]
  if( is.null(p) ) {
    ps<-estimatePs(y1,y2,x$lib.size[x$group==g[1]],x$lib.size[x$group==g[2]],r)
  } else {
    ps<-p
  }
  mu1.hat<-matrix(ps$p1,ncol=1) %*% x$lib.size[x$group==g[1]]
  mu2.hat<-matrix(ps$p2,ncol=1) %*% x$lib.size[x$group==g[2]]
  mu.hat<-matrix(ps$p,ncol=1) %*% x$lib.size
  if (length(r)==nrow(y1)) {
     r<-matrix(r,ncol=1) %*% rep(1,ncol(x$data))
  } 
  return(-2*(fullLogLik(y1,y2,mu.hat[,x$group==g[1]],mu.hat[,x$group==g[2]],r)-fullLogLik(y1,y2,mu1.hat,mu2.hat,r)))
}

fisherTest<-function(y,g,mus,r) {
  gs<-unique(g)
  y1<-y[,g==gs[1]]; if (is.vector(y1)) { y1<-matrix(y1,ncol=1) }; n1<-ncol(y1)
  y2<-y[,g==gs[2]]; if (is.vector(y2)) { y2<-matrix(y2,ncol=1) }; n2<-ncol(y2)
  pvals<-rep(NA,nrow(y1))
  v<-cbind( rowSums(y1), rowSums(y2))
  N<-floor(rowSums(v))
  mn.mx<-cbind(apply(v,1,min),apply(v,1,max))
  if (length(mus)==1) { mus<-rep(mus,nrow(y1)) }
  print("Calculating Fisher exact p-values.")
  for (i in 1:length(pvals)) {
    ind<-0:N[i]
    crit<-ind/(N[i]-ind+.000001*as.numeric(ind==0))
    p.top<-dnbinom(ind,size=n1*r[i],mu=n1*mus[i])*dnbinom(N[i]-ind,size=n2*r[i],mu=n2*mus[i])
    #p.obs<-p.top[ind==v[i,1]]
    keep<- (crit >= mn.mx[i,2]/mn.mx[i,1]) | (crit <= mn.mx[i,1]/mn.mx[i,2])
    p.bot<-dnbinom(N[i],size=(n1+n2)*r[i],mu=(n1+n2)*mus[i])
    pvals[i]<-sum(p.top[keep]/p.bot)
    if (i %% 500 == 0) { print(i) }
  } 
  pvals
}


fisherExactScore<-function(s,mu=NULL) {
  g1<-s$d[[1]]$group==1
  g2<-s$d[[1]]$group==2
  n1<-sum(g1)
  n2<-sum(g2)
  pp<-NULL
  for (j in 1:length(s$pseudo)) {
    pseudo<-s$pseudo[[j]]
    v<-cbind( rowSums(pseudo[,g1]), rowSums(pseudo[,g2]))
    N<-floor(rowSums(v))
    mn.mx<-cbind(apply(v,1,min),apply(v,1,max))
    pvals<-matrix(NA,nrow=nrow(mn.mx),ncol=ncol(s$est))
    if (is.null(mu)) {
      mus<-s$mus[[j]]
    } else {
      mus<-matrix(mu,nrow=nrow(s$mus[[j]]),ncol=ncol(s$mus[[j]]))
    }
    for (i in 1:nrow(pvals)) {
      ind<-0:N[i]
      crit<-ind/(N[i]-ind)
      keep<- (crit >= mn.mx[i,2]/mn.mx[i,1]) | (crit <= mn.mx[i,1]/mn.mx[i,2])
      for(k in 1:ncol(pvals)) {
        r<-1/s$est[j,k]-1
        p<-dnbinom(ind,n1*r,mu=n1*mus[i,1])*dnbinom(N[i]-ind,n2*r,mu=n2*mus[i,1])/dnbinom(N[i],size=(n1+n2)*r,mu=(n1+n2)*mus[i,1])
        pvals[i,k]<-sum(p[keep])
      }
    }
    pp[[j]]<-pvals
    cat("Done ",j,".\n")
  }
  pp
}


waldTest<-function(x,r,p=NULL) {
  g<-unique(x$group)
  y1<-x$data[,x$group==g[1]]
  y2<-x$data[,x$group==g[2]]
  if( is.null(p) ) {
    ps<-estimatePs(y1,y2,x$lib.size[x$group==g[1]],x$lib.size[x$group==g[2]],r)
  } else {
    ps<-p
  }
  mu1.hat<-matrix(ps$p1,ncol=1) %*% x$lib.size[x$group==g[1]]
  mu2.hat<-matrix(ps$p2,ncol=1) %*% x$lib.size[x$group==g[2]]
  A1 = rowSums(mu1.hat^2 / (mu1.hat+mu1.hat^2/r))
  A2 = rowSums(mu2.hat^2 / (mu2.hat+mu2.hat^2/r))
  beta<-log(ps$p2)-log(ps$p1)
  return(beta*sqrt(A1*A2/(A1+A2)))
}

scoreTest<-function(x,r,p=NULL) {
  g<-unique(x$group)
  y1<-x$data[,x$group==g[1]]; n1<-ncol(y1)
  y2<-x$data[,x$group==g[2]]; n2<-ncol(y2)
  if( is.null(p) ) {
    ps<-estimatePs(y1,y2,x$lib.size[x$group==g[1]],x$lib.size[x$group==g[2]],r)
  } else {
    ps<-p
  }
  m1<-outer(rep(1,nrow(x$data)),x$lib.size[x$group==g[1]])
  m2<-outer(rep(1,nrow(x$data)),x$lib.size[x$group==g[2]])
  p1<-outer(ps$p,rep(1,ncol(y1)))
  p2<-outer(ps$p,rep(1,ncol(y2)))
  if( length(r)==1 ) { r<-rep(r,nrow(x$data)) }
  if( length(r)==nrow(x$data) ) { r<-outer(r,rep(1,ncol(x$data))) }
  r1<-r[,1:n1]
  r2<-r[,1:n2]
  a1<-apply(r1*m1/((m1*p1+r1)*p1),1,sum)
  a2<-apply(r2*m2/((m2*p2+r2)*p2),1,sum)
  score2<-r2*(y2-p2*m2)/((m2*p2+r2)*p2)
  #score<-apply(score2,1,sum)
  score<-rowSums(score2)
  s<-score / sqrt( a1*a2/(a1+a2) )
  return(s)
}

logLikDerP<-function(p,y,lib.size,r,der=0) {
  # log-likelihood (assuming r known), derivatives w.r.t. p
  if (is.vector(y)) {
    y<-matrix(y,nrow=1)
  }
  n<-matrix(rep(1,nrow(y)),ncol=1) %*% lib.size
  pm<-p %*% matrix(rep(1,ncol(y)),nrow=1)
  if (length(r)==nrow(y)) {
    rm<-r %*% matrix(rep(1,ncol(y)),nrow=1)
  } else {
    rm<-matrix(r,nrow=nrow(y),ncol=ncol(y))
  }
  if (der==1) {
     ret<-(1/p)*rowSums(rm*(y-n*pm)/(n*pm+rm),na.rm=TRUE)
     ret[p==0]<-0
  } else if (der==0) {
     ret<-rowSums(-(rm+y)*log(rm+n*pm)+y*log(n*pm),na.rm=TRUE) 
     ret[p==0]<-1
  } else if (der==2) {
     ret<-(1/p)*rowSums( -(rm/p)*(y-n*pm)/(n*pm+rm) - rm*n*(y+rm)/((n*pm+rm)^2),na.rm=TRUE)
     ret[p==0]<-1
  }
  return(ret)
}

estimatePs<-function(y1,y2,lib.size1,lib.size2,r,tol=1e-10,maxit=30) {
  if (is.vector(y1) | is.vector(y2)) {
    y1<-matrix(y1,nrow=1)
    y2<-matrix(y2,nrow=1)
  }
  y<-cbind(y1,y2)
  lib.size<-c(lib.size1,lib.size2)
  this.p1<-rowMeans(y1/matrix(rep(1,nrow(y1)),ncol=1) %*% lib.size1)
  this.p2<-rowMeans(y2/matrix(rep(1,nrow(y2)),ncol=1) %*% lib.size2)
  this.p<-rowMeans(y/matrix(rep(1,nrow(y)),ncol=1) %*% lib.size)
  a<-rowSums(y1)
  b<-rowSums(y2)
  min.val<-8.783496e-16
  for(i in 1:maxit) { # do 10 Newton method steps
    d1p<-logLikDerP(this.p,y,lib.size,r,der=1)
    d1p1<-logLikDerP(this.p1,y1,lib.size1,r,der=1); d1p1[a==0]<-min.val
    d1p2<-logLikDerP(this.p2,y2,lib.size2,r,der=1); d1p2[b==0]<-min.val
    mx<-max(abs(c(d1p,d1p1,d1p2)))
    if( mx < tol ) { break }
    #cat("[estimatePs] Newton iteration (",mx,"): ",i,"\n")
    d2p<-logLikDerP(this.p,y,lib.size,r,der=2)
    d2p1<-logLikDerP(this.p1,y1,lib.size1,r,der=2)
    d2p2<-logLikDerP(this.p2,y2,lib.size2,r,der=2)
    this.p<-this.p-d1p/d2p
    this.p1<-this.p1-d1p1/d2p1
    this.p1[a==0]<-min.val
    this.p2<-this.p2-d1p2/d2p2
    this.p2[b==0]<-min.val
    this.p[this.p<=0 | this.p>=1]<-1/(max(c(lib.size1,lib.size2))*10)
    this.p1[this.p1<=0 | this.p1>=1]<-1/(max(lib.size1)*10)
    this.p2[this.p2<=0 | this.p2>=1]<-1/(max(lib.size2)*10)
  }
  return(list(p=this.p,p1=this.p1,p2=this.p2))
}

estimatePs1<-function(y,lib.size,r,tol=1e-10,maxit=30) {
  if (is.vector(y)) { y<-matrix(y,nrow=1) }
  this.p<-rowMeans( y/outer(rep(1,nrow(y)),lib.size) )
  a<-rowSums(y)
  min.val<-8.783496e-16
  for(i in 1:maxit) { # do Newton method steps
    d1p<-logLikDerP(this.p,y,lib.size,r,der=1)
    mx<-max(abs(d1p))
    if( mx < tol ) { break }
    d2p<-logLikDerP(this.p,y,lib.size,r,der=2)
    this.p<-this.p-d1p/d2p
    this.p[this.p<=0]<-1/(max(lib.size)*10)
  }
  return(list(p=this.p))
}

loglik.normal.rowsums<-function(x,mu,sd) {
  rowSums(dnorm(x,mean=mu,sd=sd,log=T))
}

findMaxD<-function(x,alpha=0.5,grid=200,grid.start=0.001,grid.end=.999,v=1) {
  # this calculates delta as combined version of overall and individual
  grid.vals<-seq(grid.start,grid.end,length=grid)
  l0<-condLogLikDer2(x,grid.vals,der=0,doSum=0)
  m0<-matrix(rep(1,nrow(x)),ncol=1) %*% colSums(l0)
  l0a<-l0+alpha*m0
  return(grid.vals[apply(l0a,1,which.max)])
}

findMaxD2<-function(x,alpha=0.5,grid=T,tol=1e-05,n.iter=5,grid.length=200) {
  # this calculates delta as combined version of overall and individual
  g<-unique(x$group)
  y1<-x$data[,x$group==g[1]]
  y2<-x$data[,x$group==g[2]]
  if(grid) {  # do a grid search, since some likelihoods may be monotone, not amenable to NR
    grid.vals<-seq(0.001,0.999,length=grid.length)
    l0<-condLogLikDer2(y1,grid.vals,der=0,doSum=0)+condLogLikDer2(y2,grid.vals,der=0,doSum=0)
    m0<-matrix(rep(1,nrow(x$data)),ncol=1) %*% colSums(l0)
    l0a<-l0+alpha*m0
    return(grid.vals[apply(l0a,1,which.max)])
  } else {  # do Newton Rhapson
    xprev<-findMaxD2(x,grid=T,alpha=0,grid.length=20)
    print(xprev)
    mx<-tol+1; iter<-0
    while( mx > tol & iter < n.iter ) {
      iter<-iter+1
      l1<-condLogLikDer2(y1,xprev,der=1,grid=F,doSum=0)+condLogLikDer2(y2,xprev,der=1,grid=F,doSum=0)
      l1<-l1+sum(alpha*l1)
      l2<-condLogLikDer2(y1,xprev,der=2,grid=F,doSum=0)+condLogLikDer2(y2,xprev,der=2,grid=F,doSum=0)
      l2<-l1+sum(alpha*l2)
      x<-xprev-l1/l2
      mx<-max(abs(x-xprev))
      xprev<-x
      print(mx)
    }
  }
  x
}

adjProfLik<-function(delta,x,null.hypothesis=T,doSum=1) {
  # calculates the cox-reid adjusted profile likelihood
  r<-1/delta-1
  g<-unique(x$group)
  k1<-x$group==g[1]; k2<-x$group==g[2]
  y1<-matrix(x$data[,k1],ncol=sum(k1)); y2<-matrix(x$data[,k2],ncol=sum(k2))
  nt<-nrow(y1); ntv<-rep(1,nt)
  ps<-estimatePs(y1,y2,x$lib.size[k1],x$lib.size[k2],r,maxit=60)
  if (null.hypothesis) {
    p1<-outer(ps$p,rep(1,sum(k1))); p2<-outer(ps$p,rep(1,sum(k2)))
  } else {
    p1<-outer(ps$p1,rep(1,sum(k1))); p2<-outer(ps$p2,rep(1,sum(k2)))
  }
  m1<-outer(ntv,x$lib.size[k1]); m2<-outer(ntv,x$lib.size[k2])
  ll<-fullLogLik(y1,y2,m1*p1,m2*p2,r) 
  adj<-log( abs(llP2(y1,m1,p1,r)*llP2(y2,m2,p2,r)) )/2
  if (doSum) { return( sum(ll-adj) ) } else { return( ll-adj ) }
}

adjProfLik1<-function(delta,y) {
  # calculates the cox-reid adjusted profile likelihood, 1 sample case
  r<-1/delta-1
  if (is.vector(y)) { y<-matrix(y,nrow=1) }
  mu<-rowMeans(y)
  return( fullLogLik1(c(mu,delta),y) - log( abs(llP2.1(y,mu,r)) )/2 )
}

adjProfLik1o<-function(delta,y,lib.size) {
  # calculates the cox-reid adjusted profile likelihood, 1 sample case
  r<-1/delta-1
  if (is.vector(y)) { y<-matrix(y,nrow=1) }
  ps<-estimatePs1(y,lib.size,r)
  return( fullLogLik1o(c(ps$p,delta),y,lib.size) - log( abs(llP2(y,lib.size,ps$p,r)) )/2 )
}

adjProfLikM<-function(delta,y,lib.size) {
  # calculates the cox-reid adjusted profile likelihood, 1 sample case
  r<-1/delta-1
  if (is.vector(y)) { y<-matrix(y,nrow=1) }
  ps<-estimatePs1(y,lib.size,r)
  mu<-outer(ps$p,lib.size)
  m<-outer(rep(1,nrow(y)),lib.size)
  p<-outer(ps$p,rep(1,ncol(y)))
  adj<-fullLogLik.delta(delta,y,mu=mu,lib.size,doSum=1)-sum(log(abs(llP2(y,m,p,r))))/2
  return(adj)
}

rl<-function() {
source("reload.R")
}


fullLogLik1o<-function(param,y,lib.size) {
  # param is vector of length 2, first element is lambda, second is delta
  mu<-lib.size*param[1]
  r<-1/param[2]-1
  # log likelihood -- 1 sample case
  if (is.vector(y)) { y<-matrix(y,nrow=1) }
  n<-dim(y)[2]
  g<-nrow(y)
  mu[mu==0]<-1e-12
  if (length(r)==nrow(y)) {
     r<-matrix(r,ncol=1) %*% rep(1,n)
  }
  ll<-rowSums(lgamma(y+r)-(r+y)*log(mu+r)-lgamma(y+1)+y*log(mu)+r*log(r)-lgamma(r))
  return(ll)
}

fullLogLik1<-function(param,y) {
  # param is vector of length 2, first element is mu, second is delta
  mu<-param[1]
  r<-1/param[2]-1
  # log likelihood -- 1 sample case
  if (is.vector(y)) { y<-matrix(y,nrow=1) }
  n<-dim(y)[2]
  g<-nrow(y)
  mu[mu==0]<-1e-12
  if (length(r)==nrow(y)) {
     r<-matrix(r,ncol=1) %*% rep(1,n)
  }
  ll<-rowSums(lgamma(y+r)-(r+y)*log(mu+r)-lgamma(y+1)+y*log(mu)+r*log(r)-lgamma(r))
  return(ll)
}

llP2.1<-function(y,mu,r) {
   # this is the negative second derivative of the log likelihood with respect to mu
   # and is therefore the OBSERVED INFORMATION
   sum( -y/mu^2 + (y+r)/((mu+r)^2) )
}

llP2<-function(y,m,p,r) {
   # this is the negative second derivative of the log likelihood with respect to p
   # and is therefore the OBSERVED INFORMATION
   rowSums( (y+r)*((m/(m*p+r))^2) - y/(p^2) )
}

doNewton2<-function(y,xprev=NULL,tol=1e-06,maxit=30,verbose=FALSE) {
  # this may not work as of the changes from may 28
  iter<-0
  if (is.null(xprev)) {
    x<-seq(0.01,0.99,length=100)
    l2<-condLogLikDer2(y,x,der=2,doSum=1)
    l1<-condLogLikDer2(y,x,der=1,doSum=1)
    l0<-condLogLikDer2(y,x,der=0,doSum=1)
    xprev<-x[which.max(l0)]
    if(verbose) {
      par(mfrow=c(3,1))
      plot(x,l0)
      plot(x,l1)
      plot(x,l2)
    }
    #return(list(l0=l0,l1=l1,l2=l2,x=x))
    if( sum(l1<0)==length(l1) ) { return(0) }
    if( sum(l1>0)==length(l1) ) { return(1) }
  } else {  
    x<-xprev
  }
  l1<-condLogLikDer2(y,x,der=1,doSum=1)
  x<-xprev-condLogLikDer2(y,xprev,der=1,doSum=1)/condLogLikDer2(y,xprev,der=2,doSum=1)
  while( abs(x-xprev) > tol & iter <= 30 ) {
    iter<-iter+1
    xprev<-x
    f1<-condLogLikDer2(y,xprev,der=1,doSum=1)
    f2<-condLogLikDer2(y,xprev,der=2,doSum=1)
    x<-xprev-f1/f2
    cat("Iteration:", iter, " L':",f1," L'':",f2," Est:", x,"\n")
    if( abs(f2) < tol ) {
      if (verbose) { cat("Second Derivative Below Tolerance\n") }
      return(x)
    }
    if(is.infinite(x)) {
      if (verbose) { cat("Estimate is infinite\n") }
      return(x)
    }
    if(is.nan(x)) {
      if (verbose) { cat("Estimate is NaN\n") }
      return(x)
    }
  }
  return(x)
}

read.sage<-function(files,skip=0) {
  tags<-NULL
  for (fn in files) {
     x<-read.table(fn,sep="\t",header=TRUE,skip=skip,comment.char="!")
     tags<-union(tags,as.character(x[,1]))
     print(length(tags))
  }
  counts<-as.data.frame(matrix(rep(0,length(tags)*length(files)),ncol=length(files)))
  for (i in seq(along=files)) {
     x<-read.table(files[i],sep="\t",header=TRUE,skip=skip,comment.char="!")
     aa<-match(as.character(x[,1]),tags)
     colnames(counts)[i]<-sub(".txt","",files)[i]
     counts[aa,i]<-x[,2]
  }
  rownames(counts)<-tags
  colS<-colSums(counts)
  return(list(data=as.matrix(counts),lib.size=colS))
}

myrnbinom<-function(n,mu=10,phi=1) {
  r<-1/phi
  p<-1/(1+mu*phi)
  return(rnbinom(n,r,p))
}

makeDataset<-function(lib.size,n1=5,n2=5,tags=10000,true=5000,dist="nbinom",disp=0.42,p.a=.0002,b.mult=4,mu=NULL,directional=F) {
  if(!is.null(mu)) { tags<-length(mu) }
  groups<-c(rep(1,n1),rep(2,n2))
  if( length(lib.size)!=(n1+n2) ) { stop("length of lib.size must be n1+n2") }
  data<-matrix(NA,nrow=tags,ncol=n1+n2)
  truth<-c(rep(1,true),rep(0,tags-true))
  if (length(disp) < tags) { disp<-rep(disp,tags/length(disp)) }
  if( dist=="nbinom" ) {
    if (is.null(mu)) {
      null.mu<-outer(rep(p.a,tags-true),lib.size)
      if (true > 0) {
        b.mult1<-rep(b.mult,true)
        r<-runif(true); b.mult1[r<=0.5]<-1/b.mult
        p.a<-rep(p.a,true)
        if( directional ) {
          mu1<-outer(p.a,lib.size[1:n1])
          mu2<-outer(p.a*b.mult,lib.size[(n1+1):(n1+n2)])
        } else {
          mu1<-outer(p.a*sqrt(b.mult1),lib.size[1:n1])
          mu2<-outer(p.a/sqrt(b.mult1),lib.size[(n1+1):(n1+n2)])
        }
      } else { 
        mu1<-mu2<-NULL
      }
    } else {
      b.mult1<-rep(b.mult,true)
      r<-runif(true); b.mult1[r<=0.5]<-1/b.mult
      mu1<-outer(mu[1:true]*sqrt(b.mult1),lib.size[1:n1])
      mu2<-outer(mu[1:true]/sqrt(b.mult1),lib.size[(n1+1):(n1+n2)])
      null.mu<-outer(mu[(true+1):tags],lib.size)
    }
    mu<-rbind(cbind(mu1,mu2),null.mu)
    for (i in 1:tags) {
      data[i,]<-myrnbinom(n1+n2,mu[i,],disp[i])
    }
  }
  colnames(data)<-paste("group",c(rep(1,n1),rep(2,n2)),".",c(1:n1,1:n2),sep="")
  rs<-rowSums(data)
  return(list(n1=n1,n2=n2,lib.size=lib.size,disp=disp[rs>0],data=data[rs>0,],b.mult=b.mult,truth=truth[rs>0],group=groups,mu=mu[rs>0,]))
}

doFDR<-function(score,truth) {
  x<-order(-score)
  selected<-1:length(x)
  falsedis<-cumsum(1-truth[x])
  return(list(selected=selected,falsedis=falsedis))
} 

doTibsFDR<-function(score,col=7) {
  B<-ncol(score)
  x<-order(-score[,col])
  g<-length(x)
  p<-rep(NA,g)
  for(i in 1:g) {
     p[x[i]]<-sum(score>score[x[i],col])
     if (i %% 100 == 0) { print(i) }
  }
  p<-p/B
  return(p)
}

glm.poisson.disp <- function(f, verbose = FALSE) {
  return(list(dispersion = 0.1, disper = 0.1))
}


runLu<-function(x,other.r=NULL) {
  # run Lu stuff on permuted columns

  # adopted from http://dulci.biostat.duke.edu/sage/sage_analysis_code.r
  # THIS LINK NO LONGER EXISTS
  #fit <- glm(tag.count~offset(log(lib.size))+group,family=poisson(link=log));
  #fit.ps <- glm.poisson.disp(fit,verbose=FALSE);
  #phi <- fit.ps$dispersion;
  #t.value <- summary(fit.ps)$coefficients[2,'z value'];

  #remove source code and change
  require(MASS)
  wald<-rep(NA,nrow(x$data))
  if (!is.null(other.r)) { wald2<-rep(NA,nrow(x$data)) }
  phi<-rep(NA,nrow(x$data))
  o<-log(x$lib.size)
  for(i in 1:nrow(x$data)) {
    f<-glm(as.numeric(x$data[i,])~offset(o)+x$group,family=poisson(link=log))
    k<-glm.poisson.disp(f,verbose=F)
    if (k$disper==0) {
      ff<-glm( as.numeric(x$data[i,])~offset(o)+x$group,family=poisson(link=log))
    } else {
      ff<-glm( as.numeric(x$data[i,])~offset(o)+x$group,family=negative.binomial(1/k$disper))
    }
    wald[i]<-summary(ff,disp=1)$coef[2,3]
    if (!is.null(other.r)) { 
      ff<-glm( as.numeric(x$data[i,])~offset(o)+x$group,family=negative.binomial(other.r[i]))
      wald2[i]<-summary(ff,disp=1)$coef[2,3]
    }
    phi[i]<-k$dispersion
    if (i %% 100==0) { cat("Lu - done: ",i,"\n") }
  }
  phitmp<-phi; phitmp[phitmp==0]<-.0001; r<-1/phitmp
  if (!is.null(other.r)) { 
    return(list(wald=wald,wald2=wald2,phi=phi))
  } else {
    return(list(wald=wald,phi=phi))
  }
}

msage<-function(x,alpha=0.5) {
  g<-unique(x$group)
  if (length(g) > 2) { stop("Can only do 2 sample comparison here.") }
  k1<-x$group==g[1]; k2<-x$group==g[2]
  x$data<-as.matrix(x$data)
  cat("Calculating shrinkage overdispersion parameters.\n")
  qA<-quantileAdjust(x,alpha=alpha)
  r<-qA$r
  ps<-qA$ps
  pseudo<-qA$pseudo
  fisher<-fisherTest(pseudo,x$group,qA$N*ps$p,r)
  return(list(delta=1/(r+1),ps=ps,r=r,pseudo=pseudo,M=qA$N,exact=fisher))
}

#change this function to not direct to kepler.R
#remake luESTDisp function to be use in the runLu function
luEstDisp<-function(y) {
  f<-glm(as.numeric(y)~1,family=poisson(link=log)) #regression model
  k<-glm.poisson.disp(f,verbose=F)
  return(k$dispersion)
} 
#this function is used for the poisson model
#estimates the dispersion of a Poisson model

