#SIMULATION
source("whiten.R")
library(fda)
library(matrixcalc)

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
