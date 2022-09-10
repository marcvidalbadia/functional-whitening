#SIMULATION in Vidal & Aguilera (2022)
source("whiten.R")
library(fda)
library(matrixcalc)
library(MASS)
library(ggplot2)

witt <- c("PCA","PCA-cor","ZCA","ZCA-cor","Cholesky")
Tab1 <- matrix(0,5,10)
rownames(Tab1) <- c("PCA","PCA-cor","ZCA","ZCA-cor","CHOL")
colnames(Tab1) <- c("ISR-B-13",
                    "ISR-B-25",
                    "trCROSS-B-13",
                    "trCROSS-B-25",
                    "CROSS2-B-13",
                    "CROSS2-B-25",
                    "trCORR-B-13",
                    "trCORR-B-25",
                    "CORR2-B-13",
                    "CORR2-B-25")

#------------------------------------------------------------------------------
#RUN SIMULATION
#Results are provided bellow
#Select parameters:

iter <- 180 #iterations
n.curves <- 1000 #number of curves
#--
Tabs <- list() #data storage
datasets <- list()

for (n.iter in 1:iter) { 
  #set.seed(1234)
  D <- 50 #grid
  x_predict <- seq(1,D,len=D)
  l <- 15
  SE <- function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2) 
  covi <- function(X, Y) outer(X, Y, SE, l) 
  COV <- covi(x_predict, x_predict)
  
  #set.seed(NULL)
  n <- n.curves 
  t <- seq(0, D, length.out = D)
  values1 <- MASS::mvrnorm(n, sqrt(2)*sin(2*pi*2*t/D), COV) 
  datas <- list(values1)
  
  for (i in 1:n) values1[i,] <- values1[i,] + rnorm(D, mean = 0, sd = 0.4)

  datas[[2]] <- rbind(values1); 
  data <- datas[[2]]
  arg <- 1:ncol(data)
  
  pp <- 0 #penalty parameter
  datasets[[n.iter]] <- data
  #TABLE
    for (basisdim in c(13,25)) {
      bbasis <- create.bspline.basis(rangeval=c(min(arg),max(arg)), nbasis=basisdim)
      fdx <- Data2fd(t(data), argvals=arg, bbasis, lambda=pp)
      if(basisdim == 13) {
        for (i in witt) {
          
          j <- which(i==witt)
          wfobj <- whiten.fd(fdx, proc=i, pre.center = T)
          Tab1[j,1]  <- mean((datas[[1]] - t(eval.fd(arg,fd(wfobj$ISR,bbasis))) )^2)
          Tab1[j,3]  <- wfobj$trCROSS
          Tab1[j,5]  <- wfobj$CROSS2
          Tab1[j,7]  <- wfobj$trCORR
          Tab1[j,9]  <- wfobj$CORR2
          
        }
      } else if (basisdim == 25) {#level basisdim
        for (i in witt) {
          
          j <- which(i==witt)
          wfobj <- whiten.fd(fdx, proc=i, pre.center = T)
          Tab1[j,2]  <- mean((datas[[1]] - t(eval.fd(arg,fd(wfobj$ISR,bbasis))) )^2)
          Tab1[j,4]  <- wfobj$trCROSS
          Tab1[j,6]  <- wfobj$CROSS2
          Tab1[j,8]  <- wfobj$trCORR
          Tab1[j,10] <- wfobj$CORR2
          
        }
      } 
  Tabs[[n.iter]] <- Tab1
    }
  print(n.iter)
}
#END SIMULATION
#------------------------------------------------------------------------------
