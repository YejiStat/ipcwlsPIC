#' @importFrom stats as.formula binomial predict sd
NULL
#' Fit the partly interval-censored AFT model with linear regressions
#'
#' Fit inverse weighted linear regression with partly interval-censored data
#'
#' @param L left-censoring time, having 0 if left-censored.
#' @param R right-censoring time, having \code{Inf} if right-censored.
#' @param delta censoring indicator, 0: exactly observed; 1: left-censored; 2: right-censored; 3: interval-censored;
#' @param x X matrix of baseline covariates.
#' @param beta0 true parameter values, including intercept, with the default being 1 for all parameters."
#' @param type penalized estimating method, default is "wlse", or "lasso", "alasso", "elastic", "scad" penalties are applicable.
#' @param wttype weight estimating method, default is "KM", Beran's nonparametric KM estimating method as "Beran", and  Ishwaran's random survival forests KM estimating method as "Ishwaran".
#' @param hlimit bandwidth value, default is 0.1
#' @param contx1.pos position of the continuous covariate of variable x used in the kernel of the Beran method. The default is 1.
#' @param contx2.pos position of the same or another continuous covariate of variable x used in the kernel of the Beran method. The default is 1.
#' @param selection Whether to use variable selection method; when using "wlse" type, the variable selection option must be set to FALSE, with the default being FALSE.
#' @param id cluster id. If the data does not have clustered structure, must set \code{id=NULL}.
#' @param index index of cluster weight, default is 1
#' @param maxit maximum time value of event time T or L and R, default is 100.
#' @param tol tolerance of the minimum threshold for the calculated weights, default is 1e-3.
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
#' see Kim et al., (2024+) for detailed method explanation.
#'
#' @references
#' Beran, R. (1981). Nonparametric Regression with Randomly Censored Survival Data. Technical Report, Univ.California, Berkeley.
#' 
#' Ishwaran, H., U. B. Kogalur, E. H. Blackstone, and M. S. Lauer (2008). Random survival forests. The annals of applied statistics. 2 (3), 841–860.
#' 
#' Pak, D., Langohr, K., Ning, J., Cort ́es Mart ́ınez, J., G ́omez Melis, G., and Shen, Y. (2020). Modeling the coronavirus disease 2019 incubation period: impact on quarantine policy. Mathematics, 8(9):1631.
#' 
#' Kim, Y., Park, S., Choi, S. (2024+). On weighted-least squares regression with doubly interval-censored COVID-19 study.
#'
#' @examples
#' \dontrun{
#' # Simulations
#' set.seed(111)
#' n = 200
#' x1 = runif(n,-1,1)
#' x2 = rbinom(n,1,0.43)
#' x = cbind(x1,x2)
#' T = 4 + x1 + x2 + rnorm(n)
#' U = 4 + (1 - 0.25*x1)*runif(n, -6, 5)
#' V = U + (1 - 0.1*x2)*runif(n, 6, 20)
#' U = exp(dplyr::case_when(TRUE ~ T, T>V ~ V, T<U ~ -Inf))
#' V = exp(dplyr::case_when(TRUE ~ T, T>V ~ Inf, T<U ~ U))
#' delta = ifelse(U==V, 0, 3)
#' picwls(L=L,R=R,delta=delta,x=x)
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
#'                                       delta = case_when(IC == 0 ~ 0,
#'                                                         IC == 1 ~ 3)
#'));
#'L=d$U;R=d$V; delta=d$delta
#'L=(log(d$U));R=log(d$V); delta=d$delta
#'x = cbind(d$x1,d$x2); id=d$id;
#'picwls(L=L,R=R,delta=delta,x=x)
#'picwls(L=L,R=R,delta=delta,x=x,id=id,index=1)
#'
#' }
#' @export
#'
#'

picpenwls=function(L,R,delta,x,beta0=rep(1,ncol(x)),type="wlse",wttype="KM",hlimit=0.1,contx1.pos=1,contx2.pos=1,selection=FALSE,id=NULL,index=1,maxit=1000,tol=1e-3){
  library(extRemes)
  library(MASS)
  library(tidyverse)
  library(survival)
  library(quantreg)
  library(ncvreg)
  library(glmnet)
  
  
  wtft = function(L,R,T=NULL,estimation=NULL,delta){
    
    if(sum(delta==3)!=0){ 
      #pic
      L = pmax(L,1e-8); R = pmax(R,1e-8); n=length(L)
      deltaL = ifelse(delta==3|delta==1,1,0)
      deltaR = ifelse(delta==3|delta==2,1,0)
      kml = survfit(Surv(-L,deltaL==0)~1)
      kmr = survfit(Surv(R,deltaR==0)~1)
      
    }else if(sum(delta==2)!=0 & sum(delta==1)!=0){ 
      #dc
      L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmin(R,T) );n=length(Y);
      deltaL = ifelse(delta==1,1,0)
      deltaR = ifelse(delta==2,1,0)
      kml = survfit(Surv(-Y,deltaL==1)~1)
      kmr = survfit(Surv(Y,deltaR==1)~1)
      
    }else if(sum(delta==2)!=0 & sum(delta==1)==0){  
      #rc
      R=pmax(R,1e-8); Y=pmin(R,T); n=length(Y)
      deltaR = ifelse(delta==2,1,0)
      kmr = survfit(Surv(Y,deltaR==1)~1)
      
    }else if(sum(delta==1)!=0){ 
      #lc
      L = pmax(L,1e-8); Y=pmax(L,T); n=length(Y)
      deltaL = ifelse(delta==1,1,0)
      kml = survfit(Surv(-Y,deltaL==1)~1)
    }
    
    ww = rep(0,n)
    for (i in 1:n) {
      if(delta[i]==0){
        
        if(sum(delta==0)==n){
          ww[i] = 1
        }else if(sum(delta==0)!=0){ 
          sl=approx(c(0,-kml$time,maxit),c(1,1-kml$surv,0), xout=L[i])$y
          sr=approx(c(0,kmr$time,maxit),c(1,kmr$surv,0), xout=R[i])$y
          ww[i] = 1/pmax(1-(sr-sl), tol)
        }else if(sum(delta==2)!=0 & sum(delta==1)!=0 & is.null(estimation)){
          sl = approx( c(0, -kml$time, maxit), c(1, 1-kml$surv,0), xout=Y[i])$y
          sr = approx( c(0, kmr$time, maxit), c(1, kmr$surv, 0), xout=Y[i])$y
          ww[i] = 1/pmax( sr-sl, tol)
        }else if(sum(delta==2)!=0 & sum(delta==1)!=0 & estimation=="DR"){ 
          sl = approx( c(0, -kml$time, maxit), c(1, 1-kml$surv,0), xout=Y[i])$y
          sr = approx( c(0, kmr$time, maxit), c(1, kmr$surv, 0), xout=Y[i])$y
          ww[i] = 1/pmax( sr, tol)
        }else if(sum(delta==2)!=0 & sum(delta==1)==0){  
          sr = approx( c(0, kmr$time, maxit), c(1, kmr$surv, 0), xout=Y[i])$y
          ww[i] = 1/pmax(sr, tol)
        }else if(sum(delta==1)!=0){ 
          sl = approx( c(0, -kml$time, maxit), c(1, 1-kml$surv,0), xout=Y[i])$y
          ww[i] = 1/pmax(1-sl, tol)
        }
      }
    }
    ww
  }
  
  
  Berfunc = function(L,R,T=NULL,x,delta) {
    
    if(sum(delta==0)==n){
      Y=T; n=length(Y); y=Y
      
    }else if(sum(delta==3)!=0){ 
      #pic
      L = pmax(L,1e-8); R = pmax(R,1e-8); Y=pmax(ifelse(delta==3,R,L),1e-8); n=length(L); y=Y
      
    }else if(sum(delta==2)!=0 & sum(delta==1)!=0){ 
      #dc
      L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmin(R,T) );n=length(Y); y=Y
      deltaL = ifelse(delta==1,1,0)
      deltaR = ifelse(delta==2,1,0)
      
    }else if(sum(delta==2)!=0 & sum(delta==1)==0){  
      #rc
      R=pmax(R,1e-8); Y=pmin(R,T); n=length(Y); y=Y
      deltaR = ifelse(delta==2,1,0)
      
    }else if(sum(delta==1)!=0){ 
      #lc
      L = pmax(L,1e-8); Y=pmax(L,T); n=length(Y); y=Y
      deltaL = ifelse(delta==1,1,0)
      
    }
    
    ker = dnorm(outer(x[,contx1.pos],x[,contx2.pos],"-")/hlimit)
    Wnj = ker / rowSums(ker)
    if(sum(delta==3)!=0){
      denomr = rowSums(outer(y,y,">=")*(Wnj))
      denoml = rowSums(outer(y,y,"<=")*(Wnj))
    }else{
      denomr = rowSums(outer(y,y,"<=")*(Wnj))
      denoml = rowSums(outer(y,y,">=")*(Wnj))
    }
    
    sr = sl = ww= rep(0,n)
    for (i in 1:n) {
      if(delta[i]==0){
        
        y0 = y[i]; nom = Wnj[,i]
        
        if(sum(delta==0)==n){
          ww[i] = 1
          
        }else if(sum(delta==3)!=0){ 
          etar = 1*(y<=y0 & delta!=0)
          etal = 1*(y>=y0 & delta!=0)
          sr = prod((1 - nom/denomr)^etar)
          sl = 1-prod((1 - nom/denoml)^etal)
          ww[i] = 1/pmax(1-(sr-sl), tol)
          
        }else if(sum(delta==2)!=0 & sum(delta==1)!=0){ 
          etar = 1*(y<=y0 & deltaR==1)
          etal = 1*(y>=y0 & deltaL==1)
          sr = prod((1 - nom/denomr)^etar)
          sl = 1-prod((1 - nom/denoml)^etal)
          ww[i] = 1/pmax( sr-sl, tol)
          
        }else if(sum(delta==2)!=0 & sum(delta==1)==0){  
          etar = 1*(y<=y0 & deltaR==1)
          sr = prod((1 - nom/denomr)^etar)
          ww[i] = 1/pmax(sr, tol)
          
        }else if(sum(delta==1)!=0){ 
          etal = 1*(y>=y0 & deltaL==1)
          sl = 1-prod((1 - nom/denoml)^etal)
          ww[i] = 1/pmax(1-sl, tol)
          
        }
      }
    }
    ww[is.na(ww)]=0
    ww
  }
  
  
  
  Ishfunc = function(L,R,T=NULL,x,delta) {
    
    if(sum(delta==3)!=0){ 
      #pic
      L = pmax(L,1e-8); R = pmax(R,1e-8); n=length(L)
      deltaL = ifelse(delta==3|delta==1,0,1)
      deltaR = ifelse(delta==3|delta==2,0,1)
      dt=data.frame(L=L,R=R,statusl=deltaL,statusr=deltaR,x=x,xx=1)
      
    }else if(sum(delta==2)!=0 & sum(delta==1)!=0){ 
      #dc
      L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmin(R,T) );n=length(Y);
      deltaL = ifelse(delta==1,0,1)
      deltaR = ifelse(delta==2,0,1)
      dt=data.frame(L=L,R=R,statusl=deltaL,statusr=deltaR,x=x,xx=1)
      
    }else if(sum(delta==2)!=0 & sum(delta==1)==0){  
      #rc
      R=pmax(R,1e-8); Y=pmin(R,T); n=length(Y); y=Y
      deltaR = ifelse(delta==2,0,1)
      dt=data.frame(R=R,statusr=deltaR,x=x,xx=1)
      
    }else if(sum(delta==1)!=0){ 
      #lc
      L = pmax(L,1e-8); Y=pmax(L,T); n=length(Y); y=Y
      deltaL = ifelse(delta==1,0,1)
      dt=data.frame(L=L,statusl=deltaL,x=x,xx=1)
    }
    
    
    if(sum(delta==0)==n){
      survl=survr=0
      
    }else if(sum(delta==3)!=0){ 
      kml.obj <- rfsrc(Surv(-L, statusl==1) ~ .-L-statusl-R-statusr, data=dt)
      # kml.obj <- predict(rfsrc(Surv(L, statusl) ~ xx, data=dt))
      kml <- get.brier.survival(kml.obj, cens.model="rfsrc")
      survl=kml$surv.aalen; survl[is.na(survl)]=0; survl
      
      kmr.obj <- rfsrc(Surv(R, statusr==0) ~ .-L-statusl-R-statusr, data=dt)
      # kmr.obj <- predict(rfsrc(Surv(R, statusr) ~ xx, data=dt))
      kmr <- get.brier.survival(kmr.obj, cens.model="rfsrc")
      survr=kmr$surv.aalen; survr[is.na(survr)]=0; survr
      
    }else if(sum(delta==2)!=0 & sum(delta==1)!=0){ 
      kml.obj <- rfsrc(Surv(-L, statusl==1) ~ .-L-statusl-R-statusr, data=dt)
      # kml.obj <- predict(rfsrc(Surv(L, statusl) ~ xx, data=dt))
      kml <- get.brier.survival(kml.obj, cens.model="rfsrc")
      survl=kml$surv.aalen; survl[is.na(survl)]=0; survl
      
      kmr.obj <- rfsrc(Surv(R, statusr==0) ~ .-L-statusl-R-statusr, data=dt)
      # kmr.obj <- predict(rfsrc(Surv(R, statusr) ~ xx, data=dt))
      kmr <- get.brier.survival(kmr.obj, cens.model="rfsrc")
      survr=kmr$surv.aalen; survr[is.na(survr)]=0; survr
      
    }else if(sum(delta==2)!=0 & sum(delta==1)==0){
      kmr.obj <- rfsrc(Surv(R, statusr==0) ~ .-R-statusr, data=dt)
      # kmr.obj <- predict(rfsrc(Surv(R, statusr) ~ xx, data=dt))
      kmr <- get.brier.survival(kmr.obj, cens.model="rfsrc")
      survr=kmr$surv.aalen; survr[is.na(survr)]=0; survr
      
    }else if(sum(delta==1)!=0){ 
      kml.obj <- rfsrc(Surv(L, statusl==0) ~ .-L-statusl, data=dt)
      # kml.obj <- predict(rfsrc(Surv(L, statusl) ~ xx, data=dt))
      kml <- get.brier.survival(kml.obj, cens.model="rfsrc")
      survl=kml$surv.aalen; survl[is.na(survl)]=0; survl
    }
    
    
    ww= rep(0,n)
    for (i in 1:n) {
      if(delta[i]==0){
        
        if(sum(delta==0)==n){
          ww[i] = 1
          
        }else if(sum(delta==3)!=0){ 
          sl = approx( c(0, (kml$event.info$time.interest), maxit), c(1, survl, 0), xout=L[i])$y
          sr = approx( c(0, (kmr$event.info$time.interest), maxit), c(1, survr, 0), xout=R[i])$y
          ww[i] = 1/pmax(1-(sr-sl),tol)
          
        }else if(sum(delta==2)!=0 & sum(delta==1)!=0){ 
          sl = approx( c(0, (kml$event.info$time.interest), maxit), c(1, survl, 0), xout=Y[i])$y
          sr = approx( c(0, (kmr$event.info$time.interest), maxit), c(1, survr, 0), xout=Y[i])$y
          ww[i] = 1/pmax( sr-sl, tol)
          
        }else if(sum(delta==2)!=0 & sum(delta==1)==0){ 
          sr = approx( c(0, (kmr$event.info$time.interest), maxit), c(1, survr, 0), xout=Y[i])$y
          ww[i] = 1/pmax(sr, tol)
          
        }else if(sum(delta==1)!=0){ 
          sl = approx( c(0, (kml$event.info$time.interest), maxit), c(1, survl, 0), xout=Y[i])$y
          ww[i] = 1/pmax(1-sl, tol)
          
        }
      }
    }
    ww[is.na(ww)]=0
    ww
  }
  
  PICls=function(L,R,x,delta,ww,eta){
    L = pmax(L,1e-8); R = (pmax(R,1e-8)); Y=pmax(ifelse(delta==0,L,R),1e-8); n=length(L);
    if(type=="oracle"){
      as.numeric(lm(Y~x[,(beta0!=0)], weights = ww*eta)$coef) #beta1, beta2
    }else{
      as.numeric(lm(Y~x, weights = ww*eta)$coef) #beta1, beta2
    }
  }
  
  Betafunc = function(L,R,x,delta,ww,eta,lambda=NULL,type){
    L = pmax(L,1e-8); R = (pmax(R,1e-8)); Y=pmax(ifelse(delta==0,L,R),1e-8); n=length(L);
    p=ncol(x); nonzero_p=sum(beta0!=0)
    init=as.numeric(lm(Y~x, weights = ww*eta)$coef)
    if (type=="wlse") {
      beta=init
      betas=NULL; lambda=NULL; best_lambda=NULL; cve=NULL; cvse=NULL
    }else if(type=="lasso"){
      beta=lassopen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta)$beta
      betas = lassopen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta)$betas
      lambda=lassopen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta)$lambda
      best_lambda=lassopen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta)$best_lambda
      cve=lassopen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta)$cve
      cvse=lassopen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta)$cve
    }else if(type=="alasso"){
      init=init[-1]
      beta=alassopen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta,init=init)$beta
      betas = alassopen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta,init=init)$betas
      lambda=alassopen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta,init=init)$lambda
      best_lambda=alassopen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta,init=init)$best_lambda
      cve=alassopen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta,init=init)$cve
      cvse=alassopen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta,init=init)$cve
    }else if(type=="elasticnet"){
      beta=elasticpen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta)$beta
      betas = elasticpen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta)$betas
      lambda=elasticpen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta)$lambda
      best_lambda=elasticpen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta)$best_lambda
      cve=elasticpen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta)$cve
      cvse=elasticpen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta)$cve
    }else if(type=="scad"){
      init=init[-1]
      beta=scadpen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta,init=init)$beta
      betas = scadpen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta,init=init)$betas
      lambda=scadpen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta,init=init)$lambda
      best_lambda=scadpen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta,init=init)$best_lambda
      cve=scadpen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta,init=init)$cve
      cvse=scadpen(L=L,R=R,x=x,delta=delta,ww=ww,eta=eta,init=init)$cve
    }
    else if(type=="oracle"){
      beta=as.numeric(lm(Y~x, weights = ww*eta)$coef)[-1]
      betas=NULL; lambda=NULL; best_lambda=NULL; cve=NULL; cvse=NULL
    }
    list(beta=beta, betas=betas, lambda=lambda, cve=cve, cvse=cvse, best_lambda=best_lambda)
  }
  
  lassopen=function(L,R,x,delta,ww,eta){
    L = pmax(L,1e-8); R = (pmax(R,1e-8)); Y=pmax(ifelse(delta==0,L,R),1e-8); n=length(L);
    p=ncol(x); nonzero_p=sum(beta0!=0)
    X=(ww*eta*x); Y=(ww*eta*Y)
    
    # Fitting a lasso regression model.
    lasso_model <- glmnet(X, Y, alpha = 1)  # alpha = 1: Lasso
    
    # Selecting the optimal lambda (Lambda) using cross-validation.
    cv_model <- cv.glmnet(X, Y, alpha = 1)
    
    # Checking the optimal lambda value.
    best_lambda <- cv_model$lambda.min
    
    # Estimating the optimal model.
    best_model <- glmnet(X, Y, alpha = 1, lambda = best_lambda)
    
    # Regression coefficients.
    est=coef(best_model)[-1]
    
    list(beta=coef(best_model)[-1],
         betas = cv_model$glmnet.fit$beta,
         lambda=cv_model$glmnet.fit$lambda,
         best_lambda=cv_model$lambda.min,
         cve=cv_model$cvm, cvse=cv_model$cvsd,
         cvup=cv_model$cvup, cvlo=cv_model$cvlo
    )
  }
  
  alassopen=function(L,R,x,delta,ww,eta,init){
    L = pmax(L,1e-8); R = (pmax(R,1e-8)); Y=pmax(ifelse(delta==0,L,R),1e-8); n=length(L);
    p=ncol(x); nonzero_p=sum(beta0!=0)
    X=(ww*eta*x); Y=(ww*eta*Y)
    
    # Fitting an alasso regression model.
    scad_model <- ncvreg(X, Y, penalty = "lasso", penalty.factor = (1/abs(init)))
    cv_model <- cv.ncvreg(X, Y, penalty = "lasso", penalty.factor = (1/abs(init)))
    best_lambda <- cv_model$lambda.min
    best_model <- ncvreg(X, Y, penalty = "lasso", lambda = best_lambda, penalty.factor = (1/abs(init)))
    # Regression coefficients.
    est=(coef(best_model)[-1])
    
    list(beta=coef(best_model)[-1],
         betas = cv_model$fit$beta,
         lambda=cv_model$fit$lambda,
         best_lambda=cv_model$lambda.min,
         cve=cv_model$cve, cvse=cv_model$cvse,
         cvup=cv_model$cve+1.96*cv_model$cvse, cvlo=cv_model$cve-1.96*cv_model$cvse
    )
  }
  
  elasticpen=function(L,R,x,delta,ww,eta){
    L = pmax(L,1e-8); R = (pmax(R,1e-8)); Y=pmax(ifelse(delta==0,L,R),1e-8); n=length(L);
    p=ncol(x); nonzero_p=sum(beta0!=0)
    X=(ww*eta*x); Y=(ww*eta*Y)
    
    # Fitting an elasticnet regression model.
    elastic_net_model <- glmnet(X, Y, alpha = 0.5)  # alpha = 0.5는 Elastic Net을 나타냄
    
    # Selecting the optimal lambda (Lambda) using cross-validation.
    cv_model <- cv.glmnet(X, Y, alpha = 0.5)
    length(cv_model$lambda)
    
    # Checking the optimal lambda value.
    best_lambda <- cv_model$lambda.min
    
    # Estimating the optimal model.
    best_model <- glmnet(X, Y, alpha = 0.5, lambda = best_lambda)
    
    # Regression coefficients.
    est=coef(best_model)[-1]
    
    list(beta=coef(best_model)[-1],
         betas = cv_model$glmnet.fit$beta,
         lambda=cv_model$glmnet.fit$lambda,
         best_lambda=cv_model$lambda.min,
         cve=cv_model$cvm, cvse=cv_model$cvsd,
         cvup=cv_model$cvup, cvlo=cv_model$cvlo
    )
  }
  
  scadpen=function(L,R,x,delta,ww,eta,init){
    L = pmax(L,1e-8); R = (pmax(R,1e-8)); Y=pmax(ifelse(delta==0,L,R),1e-8); n=length(L);
    p=ncol(x); nonzero_p=sum(beta0!=0)
    X=(ww*eta*x); Y=(ww*eta*Y)
    
    # Fitting a SCAD regression model.
    scad_model <- ncvreg(X, Y, penalty = "SCAD", penalty.factor = (1/abs(init)))
    cv_model <- cv.ncvreg(X, Y, penalty = "SCAD", penalty.factor = (1/abs(init)))
    best_lambda <- cv_model$lambda.min
    best_model <- ncvreg(X, Y, penalty = "SCAD", lambda = best_lambda, penalty.factor = (1/abs(init)))
    
    # Regression coefficients.
    est=(coef(best_model)[-1])
    
    list(beta=coef(best_model)[-1],
         betas = cv_model$fit$beta,
         lambda=cv_model$fit$lambda,
         best_lambda=cv_model$lambda.min,
         cve=cv_model$cve, cvse=cv_model$cvse,
         cvup=cv_model$cve+1.96*cv_model$cvse, cvlo=cv_model$cve-1.96*cv_model$cvse
    )
  }
  
  Amat=function(x,ww,eta,cluster){
    n=nrow(x);
    xx = as.matrix(cbind(1,x)); p=ncol(xx)
    (t(ww*eta*xx)%*%xx)/cluster
  }
  
  Mmat=function(L,R,x,delta,beta,ww,eta,cluster,type){
    L = pmax(L,1e-8); R = (pmax(R,1e-8)); Y=pmax(ifelse(delta==0,L,R),1e-8); n=length(L);
    xx = as.matrix(cbind(1,x)); p=ncol(xx)
    delta0=ifelse(delta!=0,1,0)
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
  
  L = pmax(L,1e-8); R = (pmax(R,1e-8)); Y=pmax(ifelse(delta==0,L,R),1e-8); n=length(L);
  if(is.null(id)){eta=rep(1,n); cluster=n}
  else if(index==0){eta=rep(1,n); cluster=n}
  else{ci=rep(c(table(id)),c(table(id))); wi=(1/ci); eta=(wi^(index)); cluster=length(table(id))}
  if(wttype=="KM"){ww=wtft(L=L,R=R,delta=delta);}
  if(wttype=="Ishwaran"){ww=Ishfunc(L=L,R=R,delta=delta,x=x);}
  if(wttype=="Beran"){ww=Berfunc(L=L,R=R,delta=delta,x=x);}
  # p=ncol(x); init = beta = PICls(L=L,R=R,delta=delta,x=x,ww=ww,eta=eta)
  
  if(type=="wlse"&selection==TRUE){
    new_beta = Betafunc(L=L,R=R,delta=delta,x=x,ww=ww,eta=eta,type=type)$beta[-1]
  }else if(type=="wlse"){
    new_beta = Betafunc(L=L,R=R,delta=delta,x=x,ww=ww,eta=eta,type=type)$beta
    Afunc=A=Amat(x=x,ww=ww,eta=eta,cluster=cluster); 
    Mfunc=M=Mmat(L=L,R=R,x=x,delta=delta,beta=new_beta,ww=ww,eta=eta,cluster=cluster,type=type)
    new_Sigma = diag(up_Sigma(Y=Y,Afunc=A,Mfunc=M,cluster=cluster))
    se=sqrt(new_Sigma)
  }else if(type=="oracle"){
    new_beta = Betafunc(L=L,R=R,delta=delta,x=x,ww=ww,eta=eta,type=type)$beta
    new_beta[beta0==0]=0
  }else {
    new_beta = Betafunc(L=L,R=R,delta=delta,x=x,ww=ww,eta=eta,type=type)$beta
  }
  
  if(selection!=FALSE){
    biasabs=abs(new_beta-beta0); amad_beta = mean((biasabs)); 
    new_beta_selected_vars <- new_beta != 0
    cor <- sum(new_beta_selected_vars[beta0!=0]); incor <- sum(new_beta_selected_vars) - cor
    mrme_beta <- (t(x%*%biasabs)%*%(x%*%biasabs))/(nrow(x));
    betas = Betafunc(L=L,R=R,delta=delta,x=x,ww=ww,eta=eta,type=type)$betas
    lambda = Betafunc(L=L,R=R,delta=delta,x=x,ww=ww,eta=eta,type=type)$lambda
    best_lambda = Betafunc(L=L,R=R,delta=delta,x=x,ww=ww,eta=eta,type=type)$best_lambda
    cve = Betafunc(L=L,R=R,delta=delta,x=x,ww=ww,eta=eta,type=type)$cve
  }
  
  if(type=="wlse"& selection==FALSE){
    dat=list(res=data.frame(
      est=new_beta,
      se=se,
      pvalue = 1 - pnorm(abs(new_beta/se)),
      lb = new_beta-1.96*se, ub = new_beta+1.96*se
    ),
    ww=ww
    )
    colnames(dat$res)=c("coefficients","se","pvalue","lower bd","upper bd")
  }else{
    dat=list(res=data.frame(est=new_beta), 
             # ww=ww,
             cor=cor, incor=incor, amad_beta=amad_beta, mrme_beta=mrme_beta,
             betas=betas, lambda=lambda, best_lambda=best_lambda, cve=cve)
    colnames(dat$res)=c("coefficients")
  }
  dat
}


Etable_nonpen=function(est,se,beta0){
  bias=colMeans(est)-beta0
  ase=apply(est, 2, sd)
  ese=colMeans(se)
  cpl=t( t(est)-1.96*t(se) < beta0 ) #95% coverage probability
  cpu=t( t(est)+1.96*t(se) > beta0 ) #95% coverage probability
  cp=colMeans( cpl*cpu )
  (round(cbind(bias,ase,ese,cp), 3))
}


Etable_pen=function(type,x,est,beta0,AMADest,MRMEest,corest,incorest){
  p=ncol(x)
  biasabs=abs(colMeans(est)-beta0)
  mrme_beta <- (t(x%*%biasabs)%*%(x%*%biasabs))/(nrow(x)); 
  amad_mean=mean(AMADest); amad_sd=sd(AMADest)
  mrme_mean=mean(MRMEest); mrme_sd=sd(MRMEest)
  cor_mean=mean(corest); cor_sd=sd(corest)
  incor_mean=mean(incorest); incor_sd=sd(incorest)
  (round(data.frame(mrme_beta,
                    amad_mean,amad_sd,
                    mrme_mean,mrme_sd,
                    cor_mean, cor_sd,
                    incor_mean,incor_sd), 3))
}
