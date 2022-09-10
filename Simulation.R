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
  #set.seed(NULL)
  D <- 50 #grid
  x_predict <- seq(1,D,len=D)
  l <- 15
  SE <- function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2) 
  covi <- function(X, Y) outer(X, Y, SE, l) 
  COV <- covi(x_predict, x_predict)
  
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
      } else if (basisdim == 25) {
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
#Visualize results

#Table 1
files <- c("Tabs.RData","Tabs-500.RData","Tabs-1000.RData")
load(paste("results/",files[1],sep="")) # <- Select 1 for n=180; 2, n=500; 3, n=1000
Tabs.mean <- Reduce("+",Tabs) / length(Tabs)
round(Tabs.mean[c(3,1,4,2,5),c(1,2,3,4,7,8,5,6,9,10)], 3)
#--
#Plot Figure 1
Tabs <- array(unlist(Tabs), dim = c(5,10,iter))
Score <- c()
for (j in seq(1,9,2)) for (i in 1:5) Score <- append(Score,Tabs[i,j,])
for (j in seq(2,10,2)) for (i in 1:5) Score <- append(Score,Tabs[i,j,])
df <- as.data.frame(Score)

df$Basis <- c(rep("basis 13",iter*25),rep("basis 25",iter*25))
df$Basis <- factor(df$Basis, levels=c('basis 13','basis 25'))

h <- iter*5
df$Measure <- rep(c(rep('ISR',h),rep('trCROSS',h),rep('CROSS2',h),
             rep('trCORR',h),rep('CORR2',h)),2)
df$Measure <- factor(df$Measure, levels=c('ISR','trCROSS','trCORR','CROSS2','CORR2'))
df$Operator <- rep(c(rep('PCA',iter),rep('PCA-cor',iter),rep('ZCA',iter),
                     rep('ZCA-cor',iter),rep('Cholesky',iter)), 5)
df$Operator <- factor(df$Operator, levels=c('ZCA','PCA','ZCA-cor','PCA-cor','Cholesky'))

ggplot(data = df, aes(x=Operator, y=Score)) +
geom_boxplot(aes(fatten = NULL, fill = Basis), position=position_dodge(.9),
             outlier.size = 0.08, color = "darkgrey", lwd = 0.1) +
facet_wrap( ~ Measure, scales="free") +
scale_fill_brewer(palette = 1)
