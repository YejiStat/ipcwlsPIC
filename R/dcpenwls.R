#' @importFrom stats as.formula binomial predict sd
NULL
#' Fit the doubly interval-censored AFT model with linear regression model
#'
#' Fit inverse weighted linear regression with doubly interval-censored data
#'
#' @param L left-censoring time, having 0 if left-censored.
#' @param R right-censoring time, having \code{Inf} if right-censored.
#' @param T exactly observed time.
#' @param delta censoring indicator, 1: exactly observed; 2: left-censored; 3: right-censored;
#' @param x X matrix of baseline covariates.
#' @param estimation estimating method of partly interval censored, if estimation="DR", doubly robust estimator is estimated.
#' @param wttype weight estimating method, default is "KM", Beran's nonparametric KM estimating method as "Beran", and  Ishwaran's random survival forests KM estimating method as "Ishwaran".
#' @param type penalized estimating method, default is "wls", lasso penalized estimating method as "lasso", elasticnet penalized estimating method as "elasticnet", SCAD penalized estimating method as "SCAD"
#' @param hlimit bandwidth value, default is NULL.
#' @param nonzero.index index of nonzero position, default is NULL.
#' @param lamb.len the length of grid of lambdas, default is NULL.
#' @param id cluster id. If the data does not have clustered structure, set \code{id=NULL}.
#' @param index index of cluster weight, default is 1
#' @param maxit maximum number of iteration for the log-rank estimator, default is 100.
#' @param tol tolerance of iteration for the log-rank estimator, default is 1e-3.
#'
#' @return \code{dcpenwls} returns a data frame containing at least the following components:
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
#' Kim, Y., Choi, S. (2023+). On weighted-least squares regression with partially interval-censored data.
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
#'dcpenwls(L,R,T,delta,x=x,hlimit=0.1,wttype="KM")
#' }
#' @export
#'



dcpenwls=function(L,R,T,delta,x,estimation=NULL,beta0,type,wttype="KM",hlimit=0.1,id=NULL,index=NULL,nonzero.index=NULL,lamb.len=200,maxit=100,tol=1e-10){
  library(tidyverse)
  library(extRemes)
  library(glmnet)
  library(ncvreg)
  library(MASS)
  library(survival)
  library(quantreg)
  library(randomForestSRC)
  
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
    n=length(Y);
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
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y)
    status=ifelse(delta==1,1,0)
    dt=data.frame(L=L,R=R,status=status)
    
    kml.obj <- rfsrc(Surv(L, status) ~ ., data=dt)
    kml <- get.brier.survival(kml.obj, cens.model="rfsrc")
    survl=kml$surv.aalen; survl[is.na(survl)]=0; survl
    
    kmr.obj <- rfsrc(Surv(R, status) ~ ., data=dt)
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
  
  
  DCls=function(L,R,T,x,delta,ww,eta,type,lambda=NULL,old_beta=NULL){
    options(warn=-1)
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y); p=ncol(x);
    w=ww*eta; wtx = (w*x);  wtY = (w*((Y)));
    if(type=="nonpenalty"){
      as.numeric(lm(Y~x, weights = ww*eta)$coef)[-1] #beta1, beta2
    }else if(type=="oracle"){
      as.numeric(lm(Y~x, weights = ww*eta)$coef)[-1] #beta1, beta2
    }else if(type=="elasticnet"){
      # elastic.model <- glmnet(x=wtx, y=(wtY), alpha  = 0.5,lambda = lambda, penalty.factor = (1/old_beta))
      elastic.model <- cv.glmnet(x=wtx, y=(wtY), alpha=0.5,lambda=(exp(1))^seq(-2, 3, length = lamb.len)); plot(elastic.model)
      best_lambda <- elastic.model$lambda.min
      Elasticnet_coef <- coef(elastic.model, s=best_lambda)
      list(beta=(as.numeric(round(Elasticnet_coef,3)))[-1],
           lambda=(elastic.model$lambda),
           bestlambda=elastic.model$lambda.min,
           cvm=elastic.model$cvm, cvsd=elastic.model$cvsd,
           cvup=elastic.model$cvup, cvlo=elastic.model$cvlo
      )
    }else if(type=="MCP"){
      (as.numeric(ncvreg(wtx, wtY, penalty="MCP", lambda = lambda, penalty.factor=(1/old_beta))$beta)[-1])  #beta1, beta2
    }else if(type=="lasso"){
      options(warn=-1)
      lasso.model <- cv.glmnet(x=wtx, y=(wtY), alpha=1,lambda=(exp(1))^seq(-2, 3, length = lamb.len)); plot(lasso.model)
      best_lambda <- lasso.model$lambda.min
      Lasso_coef <- coef(lasso.model, s=best_lambda)
      list(beta=(as.numeric(round(Lasso_coef,3)))[-1],
           lambda=(lasso.model$lambda),
           bestlambda=lasso.model$lambda.min,
           cvm=lasso.model$cvm, cvsd=lasso.model$cvsd,
           cvup=lasso.model$cvup, cvlo=lasso.model$cvlo
      )
    }else if(type=="SCAD"){
      (as.numeric(ncvreg(wtx, wtY, penalty="SCAD", lambda = lambda, penalty.factor=(1/old_beta))$beta)[-1])  #beta1, beta2
    }
  }
  
  
  
  BICft = function(L,R,T,x,ww,eta,delta,beta,beta0=beta0){
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y)
    p=ncol(x)
    sglamb=sglambft(L=L,R=R,T=T,x=x,ww=ww,eta=eta,delta=delta,beta=beta)
    df = sum(beta[beta!=0])
    n * log(sglamb/n) + df * log(n)
  }
  
  sglambft = function(L,R,T,x,ww,eta,delta,beta){
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y)
    p=ncol(x)
    res = as.numeric((ww*eta*Y) - (ww*eta*x)%*%beta)
    R <- rank(res)
    (sum(R*res)*2)/n
  }
  
  Betafunclamb=function(L,R,T,x,ww,eta,delta,beta,type,beta0=beta0){
    lamb=(exp(1))^seq(-2, 3, length = lamb.len); opt=NULL; penbeta=penbeta2=NULL; p=ncol(x)
    x=as.matrix(x)
    old_beta = DCls(L=L,R=R,T=T,delta=delta,x=x,ww=ww,eta=eta,type="nonpenalty")
    for (i in 1:length(lamb)) {
      new_beta=DCls(L=L,R=R,T=T,x=x,delta=delta,ww=ww,eta=eta,lambda=lamb[i],type=type,old_beta = old_beta); 
      opt=rbind(opt,BICft(L=L,R=R,T=T,delta=delta,x=x,ww=ww,eta=eta,beta=new_beta,beta0=beta0))
      # if(type=="elasticnet"|type=="lasso"){old_beta = PICls(U=U,V=V,delta=delta,x=x,ww=ww,eta=eta,type="nonpenalty")}
      # else{old_beta=(new_beta)}
      old_beta=(new_beta)
      new_beta=matrix(new_beta,nrow=1)
      penbeta=rbind(penbeta,new_beta)
      penbeta2=rbind(cbind(penbeta2),cbind(new_beta,lamb[i]))
    }
    # pen_beta=penbeta[which.min(abs(opt)),]; lambda=penbeta2[which.min(abs(opt)),(p+1)]
    pen_beta=penbeta[which.min((opt)),]; lambda=penbeta2[which.min((opt)),(p+1)]
    list(pen_beta=pen_beta, lambda=lambda, pen_beta2=penbeta2)
  }
  
  Amat=function(type,x,ww,eta,cluster){
    n=nrow(x);
    p=ncol(x)
    (t(x)%*%diag(eta)%*%x)/cluster 
  }
  
  Mmat=function(L,R,T,x,delta,beta,ww,eta,cluster,type){
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y)
    p=ncol(x)
    delta2=ifelse(delta==2,1,0); delta3=ifelse(delta==3,1,0)
    M11=ww*Y*x
    xbeta = as.numeric(x%*%beta)
    M12=(ww*   xbeta )*x
    Dr=Dl=D1r=D1l=NULL
    M4=M44=NULL
    M5=M55=NULL
    
    for(i in 1:n){
      yind=Y>=Y[i]
      denom=sum(yind*eta ) 
      D1r=(colSums(ww*yind*x*Y*eta - ww*   yind*x*xbeta*eta   )/denom)
      Dr=rbind(Dr,D1r)
      M4=(delta2[i]*yind[i]*D1r*eta[i])/(sum(yind[i]*eta ) )
      M44=rbind(M44,M4)
      D1l=(colSums(ww*(1-yind)*x*Y*eta - ww*   (1-yind)*x*xbeta*eta)/(n+1-denom))
      Dl=rbind(Dl,D1l)
      M5=(delta3[i]*(1-yind[i])*D1l*eta[i])/(n+1-sum(yind[i]*eta ) )
      M55=rbind(M55,M5)
    }
    M2=delta2*Dr
    M3=delta3*Dl
    (M=(M11-M12-M2-M3+M44+M55)*eta)
    (t(M) %*% (M))/cluster
  }
  
  up_Sigma=function(Y,Afunc, Mfunc, cluster){
    n=length(Y)
    invA = solve(Afunc)
    newSigma = (( (invA) %*% Mfunc %*% (invA) ) )
    newSigma/cluster
  }
  
  L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y)
  if(is.null(id)){eta=rep(1,n); cluster=n}
  else{ci=rep(c(table(id)),c(table(id))); wi=(1/ci); eta=(wi^(index)); cluster=length(table(id))}
  if(wttype=="KM"){ww=wtfunc(L=L,R=R,T=T,delta=delta);}
  else if(wttype=="Ishwaran"){ww=Ishrfwtfunc(L=L,R=R,T=T,delta=delta,x=x);}
  else if(wttype=="Beran" & is.null(hlimit)==F){ww=Berwtfunc(L=L,R=R,T=T,delta=delta,x=x,h=hlimit);}
  p=ncol(x)
  init = beta = DCls(L=L,R=R,T=T,delta=delta,x=x,ww=ww,eta=eta,type="nonpenalty")
  
  if(type=="wlse"){
    # new_beta = Betafunc(L=L,R=R,T=T,delta=delta,x=x,ww=ww,eta=eta,type=type,beta = init)
    new_beta = DCls(L=L,R=R,T=T,delta=delta,x=x,ww=ww,eta=eta,type="nonpenalty")
  }else if(type=="oracle"){
    new_beta = DCls(L=L,R=R,T=T,delta=delta,x=x,ww=ww,eta=eta,type="oracle")
    new_beta[-c(nonzero.index)]=0
  }
  else if(type=="elasticnet"|type=="lasso"){
    new_beta = DCls(L=L,R=R,T=T,delta=delta,x=x,ww=ww,eta=eta,type=type)$beta
    betalamb = DCls(L=L,R=R,T=T,delta=delta,x=x,ww=ww,eta=eta,type=type)$lambda
  }
  else {
    new_beta = Betafunclamb(L=L,R=R,T=T,x,delta,ww,eta,beta,type,beta0=beta0)$pen_beta
    betalamb = Betafunclamb(L=L,R=R,T=T,x,delta,ww,eta,beta,type,beta0=beta0)$pen_beta2
  }
  Afunc=A=Amat(type=type,x=x,ww=ww,eta=eta,cluster=cluster); Mfunc=M=Mmat(L=L,R=R,T=T,delta=delta,x=x,beta=new_beta,ww=ww,eta=eta,cluster=cluster,type=type)
  new_Sigma = diag(up_Sigma(Y=Y,Afunc=A,Mfunc=M,cluster=cluster))
  se=sqrt(new_Sigma)
  
  biasabs=abs(new_beta-beta0); 
  amad_beta = mean((biasabs)); 
  cor=sum(abs(new_beta[beta0!=0])>tol); incor=sum(abs(new_beta[beta0==0])>tol);
  
  dat=list(res=data.frame(
    est=new_beta,
    se=se,
    pvalue = 1 - pnorm(abs(new_beta/se)),
    lb = new_beta-1.96*se, ub = new_beta+1.96*se
  ),
  cor=cor, incor=incor, amad_beta=amad_beta)
  colnames(dat$res)=c("coefficients","se","pvalue","95% lower bd","95% upper bd")
  dat$res=round((dat$res), 6)
  dat
}


