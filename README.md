# Introduction to `ipcwlsPIC` R package


## Introduction
```{R}
`ipcwlsPIC` is the R package to fit a new inverse-probability censoring weighted (IPCW) estimating procedure for censored linear regression when the data are partially interval-censored that include doubly-censored (DC) data and partly interval-censored (PIC) and possibly correlated within the same cluster.
Let $T$ and $X$ be the event time of interest and its related $p$-vector of covariates, respectively.
Our main objective is to estimate 
the $p$-dimensional linear coefficient vector ${\boldsymbol{\beta}}_0$
in the following linear linear regression model:
$$
T_i = {\bf x}_i^T {\boldsymbol{\beta}}_0 + e_i,\quad i=1, \ldots ,n, 
$$
where $e_i$ is the random error.
When the data are subject to partially interval-censoring, 
left and right endpoints of the censoring time, $L$ and $R$,
are observed instead of $T$ such that $T\in(L,R)$.
Note that double-censoring  can also  be viewed as 
a special case of partly interval-censoring, 
i.e., $T$ is left-censored if $L=0$ and right-censored if $R=\infty$. 
```

## Description
This R package `ipcwlsPIC` implements an inverse-probability censoring weighted (IPCW) procedure for censored quantile regression for (cluster-correlated) partially interval-censored data, which includes both double-censoring and partially interval-censoring.

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
                                                        IC == 1 ~ 0)
));
L=(log(d$U));R=log(d$V); delta=d$delta
x = cbind(d$x1,d$x2); id=d$id;  tau=0.1;
ipcwlsPIC::picwls(L,R,delta,x=x)
#>           coefficients       se   pvalue  lower bd upper bd
#> Intercept     3.762554 0.289753 0.000000  3.194637 4.330470
#> 2             0.154889 0.251553 0.269037 -0.338156 0.647933
#> 3             0.386410 0.281427 0.084870 -0.165188 0.938008
ipcwlsPIC::picwls(L,R,delta,x=x,id=id,index = 1)
#>           coefficients       se   pvalue  lower bd upper bd
#> Intercept     3.762554 0.270516 0.000000  3.232342 4.292765
#> 2             0.154889 0.250467 0.268157 -0.336028 0.645805
#> 3             0.386410 0.278151 0.082384 -0.158765 0.931585
ipcwlsPIC::picwls(L,R,delta,x=x,wttype="Beran",hlimit=0.1,id=id,index = 1)
#>           coefficients       se   pvalue  lower bd upper bd
#> Intercept     3.720919 0.277275 0.000000  3.177461 4.264377
#> 2             0.126217 0.263914 0.316237 -0.391055 0.643488
#> 3             0.402683 0.288906 0.081686 -0.163573 0.968939
```


## References


* Beran, R. (1981). Nonparametric Regression with Randomly Censored Survival Data. Technical Report, Univ.California, Berkeley.

* Ishwaran, H., U. B. Kogalur, E. H. Blackstone, and M. S. Lauer (2008). Random survival forests. The annals of applied statistics. 2 (3), 841–860.

* Gorfine, M., Y. Goldberg, and Y. Ritov (2017). A quantile regression model for failure-time data with time-
dependent covariates. Biostatistics. 18 (1), 132–146.

* Pan, C. (2021). 
PICBayes: Bayesian Models for Partly Interval-Censored Data. R package. 
https://CRAN.R-project.org/package=PICBayes.

* Kim, Y., Choi, T. (2023+). 
Inverse weighted linear regression with partially interval-censored data.
*Submitted to JKSS*.
