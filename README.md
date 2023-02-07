# Introduction to `ipcwqrPIC` R package


## Introduction
`ipcwlsPIC` is the R package to fit a new inverse-probability censoring weighted (IPCW) estimating procedure for censored linear regression when the data are partially interval-censored that include doubly-censored (DC) data and partly interval-censored (PIC) and possibly correlated within same cluster.
Let $T$ and $X$ be the event time of interest and its related $p$-vector of covariates, respectively.
Our main objective is to estimate 
the $p$-dimensional quantile coefficient vector ${\boldsymbol{\beta}}_0(\tau)$
for some $\tau \in[\tau_L,\tau_R]\subset (0, 1)$ 
in the following linear regression model:

$$T_i = {\bf x}_i^T {\boldsymbol{\beta}}_0(\tau) + e_i(\tau),\quad i=1, \ldots ,n, $$

where $e_i(\tau)$ is the random error 
whose $\tau$ th quantile conditional on 
${\bf x}_i$ equals 0. 
When the data are subject to partially interval-censoring, 
left and right endpoints of the censoring time, $L$ and $R$,
are observed instead of $T$ such that $T\in(L,R)$.
Note that double-censoring  can also  be viewed as 
a special case of partly interval-censoring, 
i.e., $T$ is left-censored if $L=0$ and right-censored if $R=\infty$. 


## Description
This R package `ipcwlsPIC` implements an inverse-probability censoring weighted (IPCW) procedure for censored quantile regression for (cluster-correlated) partially interval-censored data, which includes both double-censoring and partially interval-censoring.

Vignettes is available in [here](http://htmlpreview.github.io/?https://github.com/YejiStat/ipcwqrPIC/blob/main/vignettes/ipcwlsPIC.html).


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
```


## References


* Beran, R. (1981). Nonparametric Regression with Randomly Censored Survival Data. Technical Report, Univ.California, Berkeley.

* Ishwaran, H., U. B. Kogalur, E. H. Blackstone, and M. S. Lauer (2008). Random survival forests. The annals of applied statistics. 2 (3), 841–860.

* Gorfine, M., Y. Goldberg, and Y. Ritov (2017). A quantile regression model for failure-time data with time-
dependent covariates. Biostatistics. 18 (1), 132–146.

* Chiou, S. H., Kang, S. and Yan, J. (2015). Rank-based estimating equations with general weight for accelerated failure time models: an induced smoothing approach. Statistics in Medicine 34(9): 1495–-1510.

* Zeng, D. and D. Lin (2008). Efficient resampling methods for nonsmooth estimating functions. Biostatistics 9 (2), 355–363.

* Pan, C. (2021). 
PICBayes: Bayesian Models for Partly Interval-Censored Data. R package. 
https://CRAN.R-project.org/package=PICBayes.

* Kim, Y., Choi, T., Park, S., Choi, S. and Bandyopadhyay, D. (2023+). 
Inverse weighted quantile regression with partially interval-censored data.
*Submitted to SMMR*.
