#' @importFrom stats as.formula binomial predict sd
NULL
#' Fit the doubly interval-censored AFT model with linear regression model
#'
#' Fit inverse weighted linear regression with doubly interval-censored data
#'
#' @param L left-censoring time, having 0 if left-censored.
#' @param R right-censoring time, having \code{Inf} if right-censored.
#' @param T exactly observed time.
#' @param delta censoring indicator, 1: observed; 0: interval-censored.
#' @param x X matrix of baseline covariates.
#' @param estimation estimating method of partly interval censored, if estimation="DR", doubly robust estimator is estimated.
#' @param wttype weight estimating method, default is "Param", Beran's nonparametric KM estimating method as "Beran", and  Ishwaran's random survival forests KM estimating method as "Ishwaran".
#' @param hlimit bandwidth value, default is NULL.
#' @param id cluster id. If the data does not have clustered structure, set \code{id=NULL}.
#' @param index index of cluster weight, default is 1
#' @param maxit maximum number of iteration for the log-rank estimator, default is 100.
#' @param tol tolerance of iteration for the log-rank estimator, default is 1e-3.
#'
#' @return \code{dcwls} returns a data frame containing at least the following components:
#' \itemize{
#'   \item \code{coefficients}: regression estimator.
#'   \item \code{se}: standard error estimates for \code{est}.
#'   \item \code{pvalue}: p-value.
#'   \item \code{lower bd}: lower bound of coefficients under 95% confidence level.
#'   \item \code{upper bd}: upper bound of coefficients under 95% confidence level.
#' }
#'
#' @details
#' see Kim et al., (2023+) for detailed method explanation.
#'
#' @references
#' 
#' Beran, R. (1981). Nonparametric Regression with Randomly Censored Survival Data. Technical Report, Univ.California, Berkeley.
#' 
#' Ishwaran, H., U. B. Kogalur, E. H. Blackstone, and M. S. Lauer (2008). Random survival forests. The annals of applied statistics. 2 (3), 841–860.
#' 
#' Gorfine, M., Y. Goldberg, and Y. Ritov (2017). A quantile regression model for failure-time data with time- dependent covariates. Biostatistics. 18 (1), 132–146.
#' 
#' Chiou, S. H., Kang, S. and Yan, J. (2015). Rank-based estimating equations with general weight for accelerated failure time models: an induced smoothing approach. Statistics in Medicine 34(9): 1495–-1510.
#' 
#' Pan, C. (2021). PICBayes: Bayesian Models for Partly Interval-Censored Data. R package. https://CRAN.R-project.org/package=PICBayes.
#' 
#' Kim, Y., Choi, T. (2023+). Inverse weighted linear regression with partially interval-censored data.
#' 
#'
#'
#' @examples
#' \dontrun{
#' # Simulations
#' n=200; x1=runif(n,-1.2,1.7); x2=rbinom(n,1,0.6)
#' T = 1.7+x1+x2+rnorm(n)*(1-0.1*x2)
#' L=runif(n,-2.8,1.9); R=L+runif(n,4.2,8.1)
#' Y=pmin(R,pmax(T,L))
#' delta=case_when(
#'  T<L ~ 3,
#'  T>R ~ 2,
#'  TRUE ~ 1 #observed
#')
#'L=L; R=R; T=T; delta=delta; x=cbind(x1,x2); tau=0.3
#'dcwls(L,R,T,delta,x)
#'dcwls(L,R,T,delta,x,estimation = "DR")
#'dcwls(L,R,T,delta,x,wttype = "Ishwaran")
#'dcwls(L,R,T,delta,x,wttype = "Beran",hlimit = 0.9)
#' }
#' @export
#'
#'
#'


dcwls=function(L,R,T,delta,x,estimation=NULL,wttype="Param",hlimit=NULL,id=NULL,index=NULL,maxit=100,tol=1e-3){
  library(extRemes)
  library(MASS)
  library(tidyverse)
  library(survival)
  library(tidyverse)
  library(extRemes)
  library(quantreg)
  
  wtfunc=function(L,R,T,delta){
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y)
    kml = survfit(Surv(-Y,delta==3)~1)
    kmr = survfit(Surv(Y,delta==2)~1)
    ww = rep(0,n)
    
    for (i in 1:n) {
      if(delta[i]==1){
        sl = approx( c(0, -kml$time, 100), c(1,1-kml$surv,0), xout=Y[i])$y
        sr = approx( c(0, kmr$time, 100), c(1, kmr$surv, 0), xout=Y[i])$y
        ww[i] = 1/pmax( sr-sl, 1e-3)
      }
    }
    ww
  }
  
  Rwtfunc=function(L,R,T,delta){
    
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y)
    kml = survfit(Surv(-Y,delta==3)~1)
    kmr = survfit(Surv(Y,delta==2)~1)
    wr = rep(0,n)
    for (i in 1:n) {
      if(delta[i]==1){
        sl = approx( c(0, -kml$time, 100), c(1,1-kml$surv,0), xout=Y[i])$y
        sr = approx( c(0, kmr$time, 100), c(1, kmr$surv, 0), xout=Y[i])$y
        wr[i] = 1/pmax( sr, 1e-3)
      }
    }
    wr
  }
  
  Berwtfunc = function(L,R,T,x,delta, h=NULL) {
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y); y=Y
    ker = dnorm(outer(x[,1],x[,1],"-")/h)
    Wnj = ker / rowSums(ker)
    sr = sl = srl= rep(0,n)
    denomr = rowSums(outer(y,y,"<=")*(Wnj))
    denoml = rowSums(outer(y,y,">=")*(Wnj))
    for (i in 1:n) {
      if(delta[i]==1){
        y0 = y[i]
        etar = 1*(y<=y0 & delta==2)
        etal = 1*(y>=y0 & delta==3)
        nom = Wnj[,i]
        sr = prod((1 - nom/denomr)^etar)
        sl = 1-prod((1 - nom/denoml)^etal)
        srl[i] = 1/pmax( sr-sl, 1e-3)
      }
    }
    srl
  }
  
  Ishrfwtfunc = function(L,R,T,x,delta) {
    library(randomForestSRC)
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y)
    statusl=ifelse(delta==3,0,1)
    statusr=ifelse(delta==2,0,1)
    dt=data.frame(L=L,R=R,statusl=statusl,statusr=statusr)
    
    kml.obj <- rfsrc(Surv(L, statusl) ~ ., data=dt)
    kml <- get.brier.survival(kml.obj, cens.model="rfsrc")
    survl=kml$surv.aalen; survl[is.na(survl)]=0; survl
    
    kmr.obj <- rfsrc(Surv(R, statusr) ~ ., data=dt)
    kmr <- get.brier.survival(kmr.obj, cens.model="rfsrc")
    survr=kmr$surv.aalen; survr[is.na(survr)]=0; survr
    
    ww= rep(0,n)
    for (i in 1:n) {
      if(delta[i]==1){
        sl = approx( c(0, (kml$time), 100), c(1,survl,0), xout=Y[i])$y
        sr = approx( c(0, (kmr$time), 100), c(1, survr, 0), xout=Y[i])$y
        ww[i] = 1/pmax( sr-sl, 1e-3)
      }
    }
    ww[is.na(ww)]=0
    ww
  }
  
  DCls=function(L,R,T,x,delta,ww,eta){
    
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y)
    as.numeric(lm(Y~x, weights = ww*eta)$coef) #intc, beta1, beta2
    # lm((Y)~x1+x2, weights = ww, data = d)$coef #intc, beta1, beta2
  }
  
  Efunc=function(L,R,T,x,delta,beta,ww,eta,cluster){
    
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y);
    xx = as.matrix(cbind(1,x)); p=ncol(xx)
    res = as.numeric(Y - xx%*%beta)
    U = as.vector( t(ww*xx*eta)%*%(res) )
    U/cluster
  }
  
  DREfunc=function(L,R,T,x,delta,beta,ww,wr,eta,cluster){
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y);
    wl=(ww-wr)/(ww*wr); wl[is.nan(wl)]=0; n=length(Y);
    xx = as.matrix(cbind(1,x)); p=ncol(xx)
    res = as.numeric(Y - xx%*%beta)
    U = as.vector( t(ww*xx*eta)%*%(res) )
    UR=matrix(0,p,1)
    UL=matrix(0,p,1)
    
    for(i in 1:n){
      yind=Y>=Y[i]
      denom=sum(yind*eta )
      if(delta[i]==2){
        indr=Y>=Y[i]
        dNir = Y<=Y[i]
        dMr=as.numeric(dNir-( (yind/denom) *dNir))
        Rft=as.numeric(t(xx*wr*dMr*eta)%*%( res *indr))
        UR=UR+as.numeric((Rft/n))
      }
      if(delta[i]==3){
        indl=Y<=Y[i]
        dNil = Y>=Y[i]
        dMl=as.numeric(dNil-( ((1-yind)/(n+1-denom)) *dNil))
        Lft=as.numeric(t(xx*wl*dMl*eta)%*%( res *indl))
        UL=UL+((Lft/n))
      }
    }
    UDR=as.numeric((U+UR+UL))
    UDR/cluster
  }
  
  Amat=function(x,ww,eta,cluster){
    n=nrow(x);
    xx = as.matrix(cbind(1,x)); p=ncol(xx)
    (t(ww*eta*xx)%*%xx)/cluster
  }
  
  Mmat=function(L,R,T,x,delta,beta,ww,eta,cluster){
    
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y);
    xx = as.matrix(cbind(1,x)); p=ncol(xx)
    delta2=ifelse(delta==2,1,0); delta3=ifelse(delta==3,1,0)
    res = as.numeric(Y - xx%*%beta) #nx1
    M1=ww*res*xx
    Dr=Dl=D1r=D1l=NULL
    M4=M44=NULL
    M5=M55=NULL
    
    for(i in 1:n){
      
      yind=Y>=Y[i]
      denom=sum(yind*eta ) 
      
      D1r=(colSums(ww*yind*xx*res*eta)/denom)
      Dr=rbind(Dr,D1r)
      
      M4=(delta2[i]*yind[i]*D1r*eta[i])/(sum(yind[i]*eta ) )
      M44=rbind(M44,M4)
      
      D1l=(colSums(ww*(1-yind)*xx*res*eta)/(n+1-denom))
      Dl=rbind(Dl,D1l)
      
      M5=(delta3[i]*(1-yind[i])*D1l*eta[i])/(n+1-sum(yind[i]*eta ) )
      M55=rbind(M55,M5)
    }
    M2=delta2*Dr
    M3=delta3*Dl
    (M=(M1-M2-M3+M44+M55)*eta)
    (t(M) %*% (M))/cluster
  }
  
  up_Sigma=function(Y,Afunc, Mfunc, cluster){
    n=length(Y)
    invA = solve(Afunc)
    newSigma = (( (invA) %*% Mfunc %*% (invA) ) )
    newSigma/cluster
  }
  
  
  # # # update 'Gamma' matrix in sandwich variance (Zeng and Lin, 2008)
  Gfunc2 = function(L,R,T,x,delta,ww,eta,cluster,beta,B=100) {
    n=length(L);
    library(MASS)
    Shat = t(replicate(B,{
      id = sample(n,n,replace = TRUE)
      Efunc(L=L[id],R=R[id],T=T[id],x=x[id,],delta=delta,ww=ww[id],eta=eta[id],cluster=cluster,beta = beta)
    }))
    Sigma = cov(Shat) * (cluster)
    Sigma
  }
  
  Gfunc3= function(L,R,T,x,delta,ww,eta,id,cluster,beta,B=100) {
    n=length(L)
    library(MASS)
    Shat = t(replicate(B,{
      tabid=as.vector(table(id))
      idx = as.vector(unlist(lapply(tabid, function(x) sample(x=x,size=x,replace = TRUE))))
      Efunc(L=L[idx],R=R[idx],T=T[idx],x=x[idx,],delta=delta,ww=ww[idx],eta=eta[idx],cluster=cluster,beta = beta)
    }))
    Sigma = cov(Shat) * (cluster)
    Sigma
  }
  
  L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y);
  if(is.null(id)){eta=rep(1,n); cluster=n}
  else{ci=rep(c(table(id)),c(table(id))); wi=(1/ci); eta=(wi^(index)); cluster=length(table(id))}
  
  if(wttype=="Param"){ww=wtfunc(L=L,R=R,T=T,delta=delta);}
  else if(wttype=="Ishwaran"){ww=Ishrfwtfunc(L=L,R=R,T=T,delta=delta,x=x);}
  else if(wttype=="Beran" & is.null(hlimit)==F){ww=Berwtfunc(L=L,R=R,T=T,delta=delta,x=x,h=hlimit);}
  xx = as.matrix(cbind(1,x)); p = ncol(xx)
  old_beta = init = beta = DCls(L=L,R=R,T=T,delta=delta,x=x,ww=ww,eta=eta)
  
  if(is.null(estimation)){
    new_beta = BB::dfsane(par=beta,fn=Efunc,L=L,R=R,T=T,x=x,delta=delta,ww=ww,eta=eta,cluster=cluster,control=list(trace=FALSE))$par ##IPCW
  }else if(estimation=="DR"){
    wr=Rwtfunc(L=L,R=R,T=T,delta=delta)
    new_beta = BB::dfsane(par=beta,fn=DREfunc,L=L,R=R,T=T,x=x,delta=delta,ww=ww,wr=wr,eta=eta,cluster=cluster,control=list(trace=FALSE))$par ##AIPCW
  }
  A=Amat(x=x,ww=ww,eta=eta,cluster=cluster); M=Mmat(L=L,R=R,T=T,x=x,delta=delta,beta=new_beta,ww=ww,eta=eta,cluster=cluster)
  new_Sigma = diag(up_Sigma(Y=Y,Afunc=A,Mfunc=M,cluster=cluster))
  se=sqrt(new_Sigma)
  res=data.frame(est=new_beta,
                 se=se,
                 pvalue = 1 - pnorm(abs(new_beta/se)),
                 lb = new_beta-1.96*se, ub = new_beta+1.96*se)
  colnames(res)=c("coefficients","se","pvalue","lower bd","upper bd")
  rownames(res)[1]="Intercept"
  round((res), 6)
}
