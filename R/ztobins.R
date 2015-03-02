ztobins <- function (zmat, n.association.status = 3, n.bins = 120, nulltype = 0, ...) {
  
  if (requireNamespace("asd", quietly = TRUE))
  {
    M    <- dim(zmat)[1] 
    n.studies <- dim(zmat)[2] 
    
    binned.z.mat <- array(dim=c(M,n.studies))
    dimnames(binned.z.mat) <- dimnames(zmat)
    
    pdf.binned.z  <- array(0,dim=c(n.studies, n.bins-1,n.association.status))
    dimnames(pdf.binned.z) <- list(sprintf("Study %i",1:n.studies),
                                   sprintf("bin %i",1:(n.bins-1)),
                                   if (n.association.status==2) c("H:0","H:1") else
                                     c("H:-1","H:0","H:1"))
    
    par(mfrow=c(n.studies,1))
    
    for(i in 1:n.studies)
    {
      
      expr <- expression("asd"(zmat[,i], bre = n.bins,nulltype = nulltype,
                                        main = dimnames(zmat)[[2]][i], ...)$mat)
      
      lfout <- tryCatch(eval(expr),error= function(e) e)
      
      if (is(lfout,"error"))
      {
        message("locfdr() produced the following error:")
        stop(lfout)
      }
      
      x <- lfout[,"x"]
      
      bins  <- ceiling((zmat[,i] - min(zmat[,i]))/(x[2] - x[1]))
      bins[bins == 0]   <- 1
      bins[bins ==  n.bins] <- n.bins-1
      binned.z.mat[,i] <- bins
      
      f0    <- lfout[,"f0theo"]
      p1f1  <- lfout[,"p1f1"]
      
      if (n.association.status == 2)
      {
        pdf.binned.z[i,,1]  <- f0 / sum(f0)
        pdf.binned.z[i,,2]  <- p1f1 / sum(p1f1)
      } else
        if (n.association.status == 3)
        {
          pdf.binned.z[i,,2]  <- f0 / sum(f0)
          pdf.binned.z[i,,1]  <- ifelse(x < 0,p1f1,rep(0, n.bins-1))
          pdf.binned.z[i,,3]  <- ifelse(x > 0,p1f1,rep(0, n.bins-1))
          pdf.binned.z[i,,1]  <- pdf.binned.z[i,,1] / sum(pdf.binned.z[i,,1])
          pdf.binned.z[i,,3]  <- pdf.binned.z[i,,3] / sum(pdf.binned.z[i,,3])
        } else stop("Invalide number of hypothesis states.")
    }
    par(mfrow=c(1,1))
    return(list(pdf.binned.z=pdf.binned.z,binned.z.mat=binned.z.mat))
  }
  else
  {
    stop("locfdr package is needed for this function to work. Please install it.")
  }
}