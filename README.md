# NCETS

NCOR(Negative Control Outcome Regression),is an R package for causal effect estimation of environmental exposure on health outcome. Just taking advantage of a pre-exposure outcome as an auxiliary variable, the NCOR can obtain an unbiased and robust causal effect estimation of exposure on outcome. 


# Installation
It is easy to install the development version of NCOR package using the 'devtools' package. The typical install time on a "normal" desktop computer is less than one minute.

```
# install.packages("devtools")
library(devtools)
install_github("yuyy-shandong/NCOR")
```


# Usage
There are two main functions in NCOR package. One is ncor_ind which could eliminate the unmeasured confounders and estimate causal effect for individual data. And the other one is ncor_summary which eliminate the unmeasured confounders and estimate causal effect for summary data.
You can find the instructions by '?ncor_ind'  and  '?ncor_summary'.

library(NCOR)

?ncor_ind

?ncor_summary


# Example


```
 NN <- 10
 data_total <- NULL
 for(i in 1:NN){
 c0 <- 0.5
 dia <- 0.5
 diff_c1 <- 0
 r <- 0.5
 c1 <- rnorm(1,0.5,1)
 c2 <- rnorm(1,0.5,1)
 N <- 100
 Sigma <- matrix(c(1,r,r^2,r,1,r,r^2,r,1),3,3)
 u <- MASS::mvrnorm(n=N, mu=rep(0,3), Sigma=Sigma)
 colnames(u) <- c('u1','u2','u3')
 u <- data.frame(u)
 u1 <- u$u1
 u2 <- u$u2
 u3 <- u$u3

 eps_y1 <- rnorm(N,0,0.1)
 eps_y3 <- rnorm(N,0,0.1)
 eps_x2 <-  rnorm(N,0,0.1)
 y1 <- c1*u1 + eps_y1
 x2 <- c2*u2 + eps_x2
 y3 <- c0*x2 + c1*u3+diff_c1*u3 + dia*y1 + eps_y3
 data_cenre <- data.frame(x2,u1,u2,u3,y1,y3,i)
 data_total <- rbind(data_total,data_cenre)

 }

 model <- ncor_ind (data=data_total, pre_outc_name='y1',expo_namem = 'x2',
                   post_outc_name='y3',centre='i',method='IVW',boot_no=1000)

 model

```


# Development
This R package is developed by Yuanyuan Yu, HongKai Li and Fuzhong Xue.
