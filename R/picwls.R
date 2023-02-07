#' @importFrom stats as.formula binomial predict sd
NULL
#' Fit the partly interval-censored AFT model with linear regressions
#'
#' Fit inverse weighted linear regression with partially interval-censored data
#'
#' @param L left-censoring time, having 0 if left-censored.
#' @param R right-censoring time, having \code{Inf} if right-censored.
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
#' @return \code{picwls} returns a data frame containing at least the following components:
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
#' Beran, R. (1981). Nonparametric Regression with Randomly Censored Survival Data. Technical Report, Univ.California, Berkeley.
#' 
#' Ishwaran, H., U. B. Kogalur, E. H. Blackstone, and M. S. Lauer (2008). Random survival forests. The annals of applied statistics. 2 (3), 841–860.
#' 
#' Gorfine, M., Y. Goldberg, and Y. Ritov (2017). A quantile regression model for failure-time data with time- dependent covariates. Biostatistics. 18 (1), 132–146.
#' 
#' Chiou, S. H., Kang, S. and Yan, J. (2015). Rank-based estimating equations with general weight for accelerated failure time models: an induced smoothing approach. Statistics in Medicine 34(9): 1495–-1510.
#' 
#' Zeng, D. and D. Lin (2008). Efficient resampling methods for nonsmooth estimating functions. Biostatistics 9 (2), 355–363.
#' 
#' Pan, C. (2021). PICBayes: Bayesian Models for Partly Interval-Censored Data. R package. https://CRAN.R-project.org/package=PICBayes.
#' 
#' Kim, Y., Choi, T., Park, S., Choi, S. and Bandyopadhyay, D. (2023+). Inverse weighted quantile regression with partially interval-censored data.
#'
#' @examples
#' \dontrun{
#' # Simulations
#' set.seed(111)
#' n = 200
#' x1 = runif(n,-1,1)
#' x2 = rbinom(n,1,0.43)
#' x = cbind(x1,x2)
#' T = 2 + x1 + x2 + rnorm(n)
#' U = (1 - 0.25*x1)*runif(n, -6, 5)
#' V = U + (1 - 0.1*x2)*runif(n, 6, 20)
#' U = exp(dplyr::case_when(TRUE ~ T, T>V ~ V, T<U ~ -Inf))
#' V = exp(dplyr::case_when(TRUE ~ T, T>V ~ Inf, T<U ~ U))
#' delta = ifelse(U==V, 1, 0)
#' tau=0.3
#' picwls(L=V,R=U,delta=delta,x=x)
#' picwls(L=V,R=U,delta=delta,x=x,estimation = "DR")
#' 
#' 
#' # Data example
#' library(PICBayes)
#' data("mCRC")
#' d = with(data.frame(mCRC), data.frame(U = ifelse(y==0,R,L),
#'                                       V = ifelse(y==2,L,R),
#'                                       # Cluster weighted data
#'                                       id=(rep(c(table(SITE)),c(table(SITE)))),
#'                                       # Treatment arm: 0 = FOLFIRI alone, 1 = Panitumumab + FOLFIRI.
#'                                       x1= case_when(TRT_C == 0 ~ 0, #Pan et al data
#'                                                     TRT_C == 1 ~ 1),
#'                                       # Tumor KRAS mutation status: 0 = wild-type, 1 = mutant.
#'                                       x2= case_when(KRAS_C == 0 ~ 1,
#'                                                     KRAS_C == 1 ~ 0),
#'                                       delta = case_when(IC == 0 ~ 1,
#'                                                         IC == 1 ~ 0)
#'));
#'L=d$U;R=d$V; delta=d$delta
#'L=(log(d$U));R=log(d$V); delta=d$delta
#'x = cbind(d$x1,d$x2); id=d$id;  tau=0.1;
#'picwls(L,R,delta,x=x)
#'picwls(L,R,delta,x=x,hlimit=0.1,wttype="Beran")
#'picwls(L,R,delta,x=x,estimation = "DR")
#'picwls(L,R,delta,x=x,id=id,index = 1)
#'
#'
#' }
#' @export
#'
#'

library(extRemes)

picwls=function(L,R,delta,x,estimation=NULL,wttype="Param",hlimit=NULL,id=NULL,index=NULL,maxit=100,tol=1e-3){
  library(extRemes)
  library(MASS)
  library(tidyverse)
  library(survival)
  library(tidyverse)
  library(extRemes)
  library(quantreg)
  
  wtpicft=function(L,R,delta){
    
    L = pmax(L,1e-8); R = pmax(R,1e-8); n=length(L)
    kml = survfit(Surv(L) ~ 1)
    kmr = survfit(Surv(R) ~ 1)
    ww = rep(0,n)
    
    for (i in 1:n) {
      if(delta[i]==1){
        sl=approx(c(0,kml$time,100),c(1,kml$surv,0), xout=L[i])$y
        sr=approx(c(0,kmr$time,100),c(1,kmr$surv,0), xout=R[i])$y
        ww[i] = 1/pmax((1-(sr-sl)),0.001)
      }
    }
    ww
  }
  
  
  Rwtpicft=function(L,R,delta){
    
    L = pmax(L,1e-8); R = pmax(R,1e-8); n=length(L)
    kml = survfit(Surv(L) ~ 1)
    kmr = survfit(Surv(R) ~ 1)
    ww = rep(0,n)
    
    for (i in 1:n) {
      if(delta[i]==1){
        sl=approx(c(0,kml$time,100),c(1,kml$surv,0), xout=L[i])$y
        sr=approx(c(0,kmr$time,100),c(1,kmr$surv,0), xout=R[i])$y
        ww[i] = 1/pmax((sr),0.001)
      }
    }
    ww
  }
  
  Berwtpicfunc = function(L,R,x,delta, h=NULL) {
    L = pmax(L,1e-8); R = pmax(R,1e-8); Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(L); y=Y
    ker = dnorm(outer(x[,1],x[,1],"-")/h)
    Wnj = ker / rowSums(ker)
    sr = sl = srl= rep(0,n)
    denomr = rowSums(outer(y,y,">=")*(Wnj))
    denoml = rowSums(outer(y,y,"<=")*(Wnj))
    for (i in 1:n) {
      if(delta[i]==1){
        y0 = y[i]
        etal = 1*(y>=y0 & delta==0)
        etar = 1*(y<=y0 & delta==0)
        nom = Wnj[,i]
        sr = prod((1 - nom/denomr)^etar)
        sl = 1-prod((1 - nom/denoml)^etal)
        srl[i] = 1/pmax(1-(sr-sl),0.001)
      }
    }
    srl[is.na(srl)]=0
    srl
  }
  
  Ishrfwtpicfunc = function(L,R,x,delta) {
    library(randomForestSRC)
    L = pmax(L,1e-8); R = pmax(R,1e-8); Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y); 
    statusl=ifelse(delta==0,0,1)
    statusr=ifelse(delta==0,0,1)
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
        ww[i] = 1/pmax((1-(sr-sl)),0.001)
      }
    }
    ww[is.na(ww)]=0
    ww
  }
  
  PICls=function(L,R,x,delta,ww,eta){
    
    L = pmax(L,1e-8); R = pmax(R,1e-8); Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y); 
    as.numeric(lm(Y~x, weights = ww*eta)$coef) #intc, beta1, beta2
    # lm((Y)~x1+x2, weights = ww, data = d)$coef #intc, beta1, beta2
  }
  
  Efunc=function(L,R,x,delta,beta,ww,eta,cluster){
    
    L = pmax(L,1e-8); R = pmax(R,1e-8); Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y); 
    xx = as.matrix(cbind(1,x)); p=ncol(xx)
    res = as.numeric(Y - xx%*%beta)
    U = as.vector( t(ww*xx*eta)%*%(res) )
    U/cluster
  }
  
  DREfunc=function(L,R,x,delta,beta,ww,wr,eta,cluster){
    L = pmax(L,1e-8); R = pmax(R,1e-8); Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y); 
    wl=(ww-wr)/(ww*wr); wl[is.nan(wl)]=0; n=length(Y);
    xx = as.matrix(cbind(1,x)); p=ncol(xx)
    res = as.numeric(Y - xx%*%beta)
    U = as.vector( t(ww*xx*eta)%*%(res) )
    UR=matrix(0,p,1)
    UL=matrix(0,p,1)
    
    for(i in 1:n){
      yind=Y>=Y[i]
      denom=sum(yind*eta )
      if(delta[i]==0){
        indr=Y>=Y[i]
        dNir = Y<=Y[i]
        dMr=as.numeric(dNir-( (yind/denom) *dNir))
        Rft=as.numeric(t(xx*wr*dMr*eta)%*%( res *indr))
        UR=UR+as.numeric((Rft/n))
        
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
  
  Mmat=function(L,R,x,delta,beta,ww,eta,cluster){
    
    L = pmax(L,1e-8); R = pmax(R,1e-8); Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y); 
    xx = as.matrix(cbind(1,x)); p=ncol(xx)
    delta0=ifelse(delta==0,1,0)
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
      
      M4=(delta0[i]*yind[i]*D1r*eta[i])/(sum(yind[i]*eta ) )
      M44=rbind(M44,M4)
      
      D1l=(colSums(ww*(1-yind)*xx*res*eta)/(n+1-denom))
      Dl=rbind(Dl,D1l)
      
      M5=(delta0[i]*(1-yind[i])*D1l*eta[i])/(n+1-sum(yind[i]*eta ) )
      M55=rbind(M55,M5)
    }
    M2=delta0*Dr
    M3=delta0*Dl
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
  
  L = pmax(L,1e-8); R = pmax(R,1e-8); Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y); 
  if(is.null(id)){eta=rep(1,n); cluster=n}
  else{ci=rep(c(table(id)),c(table(id))); wi=(1/ci); eta=(wi^(index)); cluster=length(table(id))}
  if(wttype=="Param"){ww=wtpicft(L=L,R=R,delta=delta);}
  else if(wttype=="Ishwaran" & n==sum(delta==1)){print("Use another weight estimating method.")}
  else if(wttype=="Ishwaran"){ww=Ishrfwtpicfunc(L=L,R=R,delta=delta,x=x);}
  else if(wttype=="Beran" & n==sum(delta==1)){print("Use another weight estimating method.")}
  else if(wttype=="Beran" & is.null(hlimit)==F){ww=Berwtpicfunc(L=L,R=R,delta=delta,x=x,h=hlimit);}
  xx = as.matrix(cbind(1,x)); p = ncol(xx)
  old_beta = init = beta = PICls(L=L,R=R,delta=delta,x=x,ww=ww,eta=eta)
  
  if(is.null(estimation)){
    new_beta = BB::dfsane(par=beta,fn=Efunc,L=L,R=R,x=x,delta=delta,ww=ww,eta=eta,cluster=cluster,control=list(trace=FALSE))$par ##IPCW
  }else if(estimation=="DR"){
    wr=Rwtpicft(L=L,R=R,delta=delta)
    new_beta = BB::dfsane(par=beta,fn=DREfunc,L=L,R=R,x=x,delta=delta,ww=ww,wr=wr,eta=eta,cluster=cluster,control=list(trace=FALSE))$par ##AIPCW
  }
  A=Amat(x=x,ww=ww,eta=eta,cluster=cluster); M=Mmat(L=L,R=R,x=x,delta=delta,beta=new_beta,ww=ww,eta=eta,cluster=cluster)
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
