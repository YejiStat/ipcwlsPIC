#' @importFrom stats as.formula binomial predict sd
NULL
#' Fit the interval-censored AFT model with linear regression model
#' 
#' Fit inverse weighted linear regression with partially interval-censored data
#'
#' @param U left-censoring time, having 0 if left-censored.
#' @param V right-censoring time, having \code{Inf} if right-censored.
#' @param delta censoring indicator, 1: exactly observed; 2: left-censored; 3: right-censored; 4: interval-censored.
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
#' @return \code{picpenwls} returns a data frame containing at least the following components:
#' \itemize{
#'   \item \code{coefficients}: regression estimator.
#'   \item \code{se}: standard error estimates for \code{est}.
#'   \item \code{pvalue}: p-value.
#'   \item \code{95% lower bd}: lower bound of coefficients under 95% confidence level.
#'   \item \code{95% upper bd}: upper bound of coefficients under 95% confidence level.
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
#' Kim, Y., Choi, T. (2023+). On weighted-least squares regression with partially interval-censored data.
#' 
#' Pan, C. (2021). PICBayes: Bayesian Models for Partly Interval-Censored Data. R package. https://CRAN.R-project.org/package=PICBayes.
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
#' delta = ifelse(U==V, 1, 4)
#' tau=0.3
#' picpenwls(U=log(U),V=log(V),delta,x=x,h=0.1,wttype="KM")
#' 
#' 
#' # Data example
#' library(PICBayes)
#' library(tidyverse)
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
#'U=(log(d$U));V=log(d$V); delta=d$delta
#'x = cbind(d$x1,d$x2); id=d$id;  tau=0.1;
#'picpenwls(U=U,V=V,delta,x=x,h=0.1,wttype="KM")
#' }
#' @export
#'
#'


picpenwls=function(U,V,delta,x,estimation=NULL,beta0,type="wls",wttype="KM",hlimit=NULL,id=NULL,index=NULL,nonzero.index=NULL,lamb.len=200,maxit=100,tol=1e-2){
  library(tidyverse)
  library(extRemes)
  library(glmnet)
  library(ncvreg)
  library(MASS)
  library(survival)
  library(quantreg)
  library(randomForestSRC)
  
  wtpicft=function(U,V,delta){
    
    U = pmax(U,1e-8); V = pmax(V,1e-8); Y=pmax(ifelse(delta==3,V,U),1e-8); n=length(U);
    kml = survfit(Surv(-Y,delta==2)~1)
    kmr = survfit(Surv(Y,delta==3)~1)
    
    ww = rep(0,n)
    
    for (i in 1:n) {
      if(delta[i]==1){
        sl=approx(c(0,kml$time,100),c(1,kml$surv,0), xout=U[i])$y
        sr=approx(c(0,kmr$time,100),c(1,kmr$surv,0), xout=V[i])$y
        ww[i] = 1/pmax((1-(sr-sl)),0.001)
      }
    }
    ww[is.na(ww)]=0;
    ww
  }
  
  
  Vwtpicft=function(U,V,delta){
    
    U = pmax(U,1e-8); V = pmax(V,1e-8); Y=pmax(ifelse(delta==3,V,U),1e-8); n=length(U);
    kml = survfit(Surv(-Y,delta==2)~1)
    kmr = survfit(Surv(Y,delta==3)~1)
    
    ww = rep(0,n)
    
    for (i in 1:n) {
      if(delta[i]==1){
        sl=approx(c(0,kml$time,100),c(1,kml$surv,0), xout=U[i])$y
        sr=approx(c(0,kmr$time,100),c(1,kmr$surv,0), xout=V[i])$y
        ww[i] = 1/pmax(sr,0.001)
      }
    }
    ww[is.na(ww)]=0;
    ww
  }
  
  Berwtpicfunc = function(U,V,x,delta, h=NULL) {
    U = pmax(U,1e-8); V = pmax(V,1e-8); Y=pmax(ifelse(delta==3,V,U),1e-8); n=length(U); y=Y
    ker = dnorm(outer(x[,1],x[,1],"-")/h) #x1: continuous variable
    Wnj = ker / rowSums(ker)
    sr = sl = srl= rep(0,n)
    denomr = rowSums(outer(y,y,">=")*(Wnj))
    denoml = rowSums(outer(y,y,"<=")*(Wnj))
    for (i in 1:n) {
      if(delta[i]==1){
        y0 = y[i]
        etal = 1*(y>=y0 & delta!=1)
        etar = 1*(y<=y0 & delta!=1)
        nom = Wnj[,i]
        sr = prod((1 - nom/denomr)^etar)
        sl = 1-prod((1 - nom/denoml)^etal)
        srl[i] = 1/pmax(1-(sr-sl),0.001)
      }
    }
    srl[is.na(srl)]=0
    srl
  }
  
  Ishrfwtpicfunc = function(U,V,x,delta) {
    U = pmax(U,1e-8); V = pmax(V,1e-8); Y=pmax(ifelse(delta==3,V,U),1e-8); n=length(U);
    status=ifelse(delta==1,1,0)
    dt=data.frame(U=U,V=V,status=status)
    
    kml.obj <- rfsrc(Surv(U, status) ~ ., data=dt)
    kml <- get.brier.survival(kml.obj, cens.model="rfsrc")
    survl=kml$surv.aalen; survl[is.na(survl)]=0; survl
    
    kmr.obj <- rfsrc(Surv(V, status) ~ ., data=dt)
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
  
  
  PICls=function(U,V,x,delta,ww,eta,type,lambda,old_beta=NULL){
    options(warn=-1)
    U = pmax(U,1e-8); V = pmax(V,1e-8); Y=pmax(ifelse(delta==3,V,U),1e-8); n=length(U); p=ncol(x);
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
    }
    else if(type=="MCP"){
      (as.numeric(ncvreg(wtx, wtY, penalty="MCP", lambda = lambda, penalty.factor=(1/old_beta))$beta)[-1])  #beta1, beta2
    }
    else if(type=="lasso"){
      # (as.numeric(ncvreg(wtx, wtY, penalty="lasso", lambda = lambda, penalty.factor=(1/old_beta))$beta)[-1])  #beta1, beta2
      # elastic.model <- glmnet(x=wtx, y=(wtY), alpha  = 1,lambda = lambda, penalty.factor = (1/old_beta))
      # (as.numeric(round(elastic.model$beta,3)))
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
  
  BICft = function(U,V,x,delta,ww,eta,beta,beta0=beta0){
    U = pmax(U,1e-8); V = pmax(V,1e-8); Y=pmax(ifelse(delta==3,V,U),1e-8); n=length(U);
    p=ncol(x)
    sglamb=sglambft(U=U,V=V,x=x,ww=ww,eta=eta,delta=delta,beta=beta)
    df = sum(beta0[beta0!=0])
    n * log(sglamb/n) + df * log(n)
  }
  
  sglambft = function(U,V,x,delta,ww,eta,beta){
    U = pmax(U,1e-8); V = pmax(V,1e-8); Y=pmax(ifelse(delta==3,V,U),1e-8); n=length(U);
    p=ncol(x)
    res = as.numeric((ww*eta*Y) - (ww*eta*x)%*%beta)
    R <- rank(res)
    (sum(R*res)*2)/n
  }
  
  Betafunclamb=function(U,V,x,delta,ww,eta,beta,type,lamb_len=lamb.len,beta0=beta0){
    lamb=(exp(1))^seq(-2, 3, length = lamb_len); opt=NULL; penbeta=penbeta2=NULL; p=ncol(x)
    old_beta = PICls(U=U,V=V,delta=delta,x=x,ww=ww,eta=eta,type="nonpenalty")
    for (i in 1:length(lamb)) {
      new_beta=PICls(U=U,V=V,delta=delta,x=x,type=type,ww=ww,eta=eta,lambda=lamb[i], old_beta=old_beta);
      opt=rbind(opt,BICft(U=U,V=V,delta=delta,x=x,ww=ww,eta=eta,beta=new_beta,beta0=beta0))
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
  
  Mmat=function(U,V,x,delta,beta,ww,eta,cluster,type){
    
    U = pmax(U,1e-8); V = pmax(V,1e-8); Y=pmax(ifelse(delta==3,V,U),1e-8); n=length(U);
    p=ncol(x)
    delta0=ifelse(delta!=1,1,0)
    res = as.numeric(Y - x%*%beta) #nx1
    M1=ww*res*x
    Dr=Dl=D1r=D1l=NULL
    M4=M44=NULL
    M5=M55=NULL
    
    for(i in 1:n){
      yind=Y>=Y[i]
      denom=sum(yind*eta ) 
      D1r=(colSums(ww*yind*x*res*eta)/denom)
      Dr=rbind(Dr,D1r)
      M4=(delta0[i]*yind[i]*D1r*eta[i])/(sum(yind[i]*eta ) )
      M44=rbind(M44,M4)
      D1l=(colSums(ww*(1-yind)*x*res*eta)/(n+1-denom))
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
  
  U = pmax(U,1e-8); V = pmax(V,1e-8); Y=pmax(ifelse(delta==3,V,U),1e-8); n=length(U);
  if(is.null(id)){eta=rep(1,n); cluster=n}
  else{ci=rep(c(table(id)),c(table(id))); wi=(1/ci); eta=(wi^(index)); cluster=length(table(id))}
  if(wttype=="KM"){ww=wtpicft(U=U,V=V,delta=delta);}
  else if(wttype=="Ishwaran"){ww=Ishrfwtpicfunc(U=U,V=V,delta=delta,x=x);}
  else if(wttype=="Beran"){ww=Berwtpicfunc(U=U,V=V,delta=delta,x=x,h=hlimit);}
  p=ncol(x)
  init = beta = PICls(U=U,V=V,delta=delta,x=x,ww=ww,eta=eta,type="nonpenalty")
  
  if(type=="wlse"){
    new_beta = PICls(U=U,V=V,delta=delta,x=x,ww=ww,eta=eta,type="nonpenalty")
  }else if(type=="oracle"){
    new_beta = PICls(U=U,V=V,delta=delta,x=x,ww=ww,eta=eta,type="oracle")
    new_beta[-c(nonzero.index)]=0
  }
  else if(type=="elasticnet"|type=="lasso"){
    new_beta = PICls(U=U,V=V,delta=delta,x=x,ww=ww,eta=eta,type=type)$beta
    betalamb = PICls(U=U,V=V,delta=delta,x=x,ww=ww,eta=eta,type=type)$lambda
  }
  else {
    new_beta = Betafunclamb(U,V,x,delta,ww,eta,beta,type,lamb_len=lamb.len,beta0=beta0)$pen_beta
    betalamb = Betafunclamb(U,V,x,delta,ww,eta,beta,type,lamb_len=lamb.len,beta0=beta0)$pen_beta2
  }
  Afunc=A=Amat(type=type,x=x,ww=ww,eta=eta,cluster=cluster); Mfunc=M=Mmat(U=U,V=V,delta=delta,x=x,beta=new_beta,ww=ww,eta=eta,cluster=cluster,type=type)
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



