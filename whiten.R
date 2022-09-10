whiten.fd  <- function (fdx, 
                        proc = c("PCA",
                                 "PCA-cor",
                                 "ZCA",
                                 "ZCA-cor",
                                 "Cholesky"),
                        pre.center = FALSE,
                        post.center = FALSE) {
 
  #################################################################
  # whiten.fd: R function to whiten functional data
  # fdx............a functional data object created with the fda package
  # proc...........whitening procedure
  # pre.center.....center the data before whitening
  # post.center....center the data after whitening
  #################################################################
  
  PD = function(U) return( sweep(U, 2, sign(diag(U)), "*") )
  
  if (!(inherits(fdx, "fd"))) 
    stop("Argument FD  not a functional data object. See fda package")
  
  anc <- mean.fd(fdx)$coefs
  if (pre.center) fdx <- center.fd(fdx);
  
  a <- fdx$coefs
  nrep <- dim(a)[2]
  phi <- fdx$basis
  nbasis <- fdx$basis$nbasis
  typebasis <- fdx$basis$type
  G <- inprod(phi, phi)
  GChol <- chol(G)
  iGChol <- solve(GChol)
  arg <- fdx$fdnames$time
  
  covc <- tcrossprod(a)/nrep
  C2 <- GChol %*% covc %*% t(GChol)
  C2 <- (C2 + t(C2))/2
  dC2 <- eigen(C2, symmetric = TRUE)
  
  if (proc == "PCA") {
    
    U <- PD(dC2$vectors)
    w <- diag(1/sqrt(dC2$values)) %*% t(U)
    
  } else if (proc == "PCA-cor") {
    
    v <- diag(C2)
    R <- cov2cor(C2)
    eR <- eigen(R, symmetric = TRUE)
    varphi <- PD(eR$vectors)
    w <- diag(1/sqrt(eR$values)) %*% t(varphi) %*% diag(1/sqrt(v))
    
  } else if (proc == "ZCA") {
    
    U <- PD(dC2$vectors)
    w <- U %*% diag(1/sqrt(dC2$values)) %*% t(U)
    
  } else if (proc == "ZCA-cor") {
    
    v <- diag(C2)
    R <- cov2cor(C2)
    eR <- eigen(R, symmetric = TRUE)
    varphi <- PD(eR$vectors)
    w <- varphi %*% diag(1/sqrt(eR$values)) %*% t(varphi) %*% diag(1/sqrt(v))
    
  } else if (proc == "Cholesky") {
    
    w <- chol(solve(C2))
    
  } 
  
  #Mapping coefficients
  wa <- iGChol %*% w %*% GChol %*% a
  
  #Identity validation (uncomment to check)
  #print(diag(JChol%*%crossprod(t(wa))%*%t(JChol))/nrep)
  
  #Whitened functional sample & post-centering
  wfdx <- fd(wa,phi)
  if (post.center) wfdx <- center.fd(wfdx);
  
  #ISR
  nb <- phi$nbasis
  W.in <- diag(1/sqrt(diag(abs(C2))))
  ISR <- (iGChol %*% W.in %*% w %*% GChol %*% a) + matrix(rep(anc,nrep),nb,nrep)
  
  #--
  cross <- tcrossprod(a,wfdx$coefs)/nrep
  CC <- GChol %*% cross %*% t(GChol)
  CC <- (CC + t(CC))/2
  CR <- cor(t(GChol%*%wfdx$coefs),t(GChol%*%a))
  CR <- (CR + t(CR))/2
  
  #Trace & Compression indexes
  trCROSS <- matrixcalc::hilbert.schmidt.norm(CC) 
  CROSS2 <- max(abs(diag(CC%*%t(CC))))
  
  trCORR <- matrixcalc::hilbert.schmidt.norm(CR) 
  CORR2 <- max(abs(diag(CR%*%t(CR))))

  WFDX <- list(wfdx, ISR, trCROSS, CROSS2, trCORR, CORR2)
  names(WFDX) <- c("wfdx", "ISR", "trCROSS", "CROSS2", "trCORR", "CORR2")
  return(WFDX)
}
