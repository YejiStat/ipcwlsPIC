# Introduction to `ipcwlsPIC` R package


## Introduction
`ipcwlsPIC` is the R package to fit a new inverse-probability censoring weighted (IPCW) estimating procedure for censored linear regression when the data are partially interval-censored that include doubly-censored (DC) data and partly interval-censored (PIC) and possibly correlated within the same cluster.
Let $T$ and $X$ be the event time of interest and its related $p$-vector of covariates, respectively.
Our main objective is to estimate 
the $p$-dimensional linear coefficient vector ${\boldsymbol{\beta}}_0$
in the following linear linear regression model:

$$T_i = {\bf x}_i^T {\boldsymbol{\beta}}_0 + \epsilon_i,\quad i=1, \ldots ,n, $$

where $e_i$ is the random error.
When the data are subject to partially interval-censoring, 
left and right endpoints of the censoring time, $L$ and $R$,
are observed instead of $T$ such that $T\in(L,R)$.
Note that double-censoring  can also  be viewed as 
a special case of partly interval-censoring, 
i.e., $T$ is left-censored if $L=0$ and right-censored if $R=\infty$. 


## Description
This R package `ipcwlsPIC` implements an inverse-probability censoring weighted (IPCW) procedure for weighted linear regression for (cluster-correlated) partially interval-censored data, which includes both double-censoring and partially interval-censoring.

Vignettes is available in [here](http://htmlpreview.github.io/?https://github.com/YejiStat/ipcwlsPIC/blob/main/vignettes/ipcwlsPIC.html).


## Usages 
```{r}
library(PICBayes)
data("mCRC")
d = with(data.frame(mCRC), data.frame(U = ifelse(y==0,R,L),
                                      V = ifelse(y==2,L,R),
                                      # Cluster weighted data
                                      id=(rep(c(table(SITE)),c(table(SITE)))),
                                      # Treatment arm: 0 = FOLFIRI alone, 1 = Panitumumab + FOLFIRI.
                                      x1= case_when(TRT_C == 0 ~ 0, #Pan et al data
                                                    TRT_C == 1 ~ 1),
                                      # Tumor KRAS mutation status: 0 = wild-type, 1 = mutant.
                                      x2= case_when(KRAS_C == 0 ~ 1,
                                                    KRAS_C == 1 ~ 0),
                                      delta = case_when(IC == 0 ~ 1,
                                                        IC == 1 ~ 4)
));
L=(log(d$U));R=log(d$V); delta=d$delta
x = cbind(d$x1,d$x2); id=d$id;
ipcwlsPIC::picpenwls(L=L,R=R,delta=delta,x=x,beta0 = c(1,1,1),type="wlse")$res
#>   coefficients        se    pvalue   lower bd  upper bd
#> 1   4.34668364 0.2039375 0.0000000  3.9469661 4.7464012
#> 2   0.10503791 0.1652651 0.2625276 -0.2188816 0.4289575
#> 3   0.08319163 0.1980515 0.3372243 -0.3049893 0.4713725
ipcwlsPIC::picpenwls(L=L,R=R,delta=delta,x=x,beta0 = c(1,1,1),type="wlse",index = 1)$res
#>   coefficients        se    pvalue   lower bd  upper bd
#> 1   4.34668364 0.2039375 0.0000000  3.9469661 4.7464012
#> 2   0.10503791 0.1652651 0.2625276 -0.2188816 0.4289575
#> 3   0.08319163 0.1980515 0.3372243 -0.3049893 0.4713725
ipcwlsPIC::picpenwls(L=L,R=R,delta=delta,x=x,beta0 = c(1,1,1),type="wlse",wttype="Beran",hlimit=0.1,id=id,index = 1)$res
#>   coefficients        se     pvalue   lower bd  upper bd
#> 1    3.7209189 0.2772746 0.00000000  3.1774607 4.2643771
#> 2    0.1262165 0.2639139 0.31623655 -0.3910548 0.6434879
#> 3    0.4026830 0.2889063 0.08168615 -0.1635733 0.9689393
```


## References


* Beran, R. (1981). Nonparametric Regression with Randomly Censored Survival Data. Technical Report, Univ.California, Berkeley.

* Ishwaran, H., U. B. Kogalur, E. H. Blackstone, and M. S. Lauer (2008). Random survival forests. The annals of applied statistics. 2 (3), 841–860.

* Pan, C. (2021). 
PICBayes: Bayesian Models for Partly Interval-Censored Data. R package. 
https://CRAN.R-project.org/package=PICBayes.

* Kim, Y., Choi, S. (2023+). 
Kim, Y., Choi, T. (2023+). On weighted-least squares regression with partially interval-censored data.
*Scheduled to aaa*.
