bkde <- function(x,kernel="normal",canonical=F,bandwidth,
                 gridsize=401,range.x,truncate=T)

# Last changed: 16/06/95

{
   # Rename common variables

   n <- length(x)
   M <- gridsize
   
   # Set canonical scaling factors   
   
   if (kernel=="normal") del0 <- (1/(4*pi))^(1/10)
   if (kernel=="box") del0 <- (9/2)^(1/5)
   if (kernel=="epanech") del0 <- 15^(1/5)
   if (kernel=="biweight") del0 <- 35^(1/5)
   if (kernel=="triweight") del0 <- (9450/143)^(1/5)

   # Set default bandwidth

   if (missing(bandwidth))
   {
      if (canonical) 
      {
         bandwidth <- (243/(35*n))^(1/5)*sqrt(var(x))
      }
      else
      {
         bandwidth <- del0*(243/(35*n))^(1/5)*sqrt(var(x))
      }
   }
   h <- bandwidth  
 
   # Set kernel support values

   if (canonical)
   {
      if (kernel=="normal") {tau <- 4*del0} 
      else {tau <- del0}
   }
   else
   {
      if (kernel=="normal") {tau <- 4} 
      else {tau <- 1}
   }

   if (missing(range.x))
   {
      range.x <- c(min(x)-tau*h,max(x)+tau*h)  
   }
   a <- range.x[1]
   b <- range.x[2]

   # Set up grid points and bin the data
  
   gpoints <- seq(a,b,length=M)
   gcounts <- linbin(x,gpoints,truncate)

   # Compute kernel weights 

   L <- min(floor(tau*h*(M-1)/(b-a)),M)
   lvec <- (0:L)
   delta  <- (b-a)/(h*(M-1))
   if (canonical==F) del0 <- 1
   if (kernel=="normal")
   {
      kappa <- dnorm(lvec*delta/del0)/(n*h*del0)
   }
   else if (kernel=="box")
   {
      kappa <- 0.5*dbeta(0.5*(lvec*delta/del0+1),1,1)/(n*h*del0)
   }
   else if (kernel=="epanech")
   {
      kappa <- 0.5*dbeta(0.5*(lvec*delta/del0+1),2,2)/(n*h*del0)
   }
   else if (kernel=="biweight")
   {
      kappa <- 0.5*dbeta(0.5*(lvec*delta/del0+1),3,3)/(n*h*del0)
   }
   else if (kernel=="triweight") 
   {
      kappa <- 0.5*dbeta(0.5*(lvec*delta/del0+1),4,4)/(n*h*del0)
   }

   # Now combine weight and counts to obtain estimate

   P <- 2^(ceiling(log(M+L)/log(2)))
   kappa <- c(kappa,rep(0,P-2*L-1),kappa[(L+1):2])
   gcounts <- c(gcounts,rep(0,P-M))
   kappa <- fft(kappa)
   gcounts <- fft(gcounts)
   return(list(x=gpoints,y=(Re(fft(kappa*gcounts,T))/P)[1:M])) 
}


bkde.surv <- function(x,cens,kernel="normal",bandwidth,
                 gridsize=401,range.x=c(0,max(x)),truncate=T)

# Last changed: 19/03/98

{
   # Set up jump vector

   
   s.fit <- surv.fit(x,cens)
   jumps <- -c(0,diff(s.fit$surv))

   # Rename common variables

   n <- length(x)
   M <- gridsize
   
   # Set canonical scaling factors   
   
   if (kernel=="normal") del0 <- (1/(4*pi))^(1/10)
   if (kernel=="box") del0 <- (9/2)^(1/5)
   if (kernel=="epanech") del0 <- 15^(1/5)
   if (kernel=="biweight") del0 <- 35^(1/5)
   if (kernel=="triweight") del0 <- (9450/143)^(1/5)

   # Set default bandwidth

   if (missing(bandwidth))
   {
      bandwidth <- del0*(243/(35*n))^(1/5)*sqrt(var(x))
   }
   h <- bandwidth  
 
   # Set kernel support values

   if (kernel=="normal") {tau <- 4} 
   else {tau <- 1}

   # Set default evaluation grid.

   if (missing(range.x))
   {
      range.x <- c(0,max(x)+tau*h)  
   }
   a <- range.x[1]
   b <- range.x[2]

   # Set up grid points and bin the data
  
   gpoints <- seq(a,b,length=M)

   counts <- rlbin(x,jumps,gpoints,truncate)

   jcounts <- counts$ycounts

   # Compute kernel weights 

   L <- min(floor(tau*h*(M-1)/(b-a)),M)
   lvec <- (0:L)
   delta  <- (b-a)/(h*(M-1))
   del0 <- 1
   if (kernel=="normal")
   {
      kappa <- dnorm(lvec*delta/del0)/(n*h*del0)
   }
   else if (kernel=="box")
   {
      kappa <- 0.5*dbeta(0.5*(lvec*delta/del0+1),1,1)/(n*h*del0)
   }
   else if (kernel=="epanech")
   {
      kappa <- 0.5*dbeta(0.5*(lvec*delta/del0+1),2,2)/(n*h*del0)
   }
   else if (kernel=="biweight")
   {
      kappa <- 0.5*dbeta(0.5*(lvec*delta/del0+1),3,3)/(n*h*del0)
   }
   else if (kernel=="triweight") 
   {
      kappa <- 0.5*dbeta(0.5*(lvec*delta/del0+1),4,4)/(n*h*del0)
   }

   # Now combine weight and counts to obtain estimate

   P <- 2^(ceiling(log(M+L)/log(2)))
   kappa <- c(kappa,rep(0,P-2*L-1),kappa[(L+1):2])
   jcounts <- c(jcounts,rep(0,P-M))
   kappa <- fft(kappa)
   jcounts <- fft(jcounts)
   return(list(x=gpoints,y=(Re(fft(kappa*jcounts,T))/P)[1:M],
               bandwidth=bandwidth)) 
}


bkde2D <- function(x,bandwidth,gridsize=c(51,51),range.x,truncate=T)

# Last changed: 25/08/95

{
   # Rename common variables

   n <- nrow(x)
   M <- gridsize
   h <- bandwidth
   tau <- 3.4      # For bivariate normal kernel.

   # Use same bandwidth in each
   # direction if only a single
   # bandwidth is given.

   if (length(h)==1) h <- c(h,h)

   # If range.x is not specified then
   # set it at its default value.
      
   if (missing(range.x))
   {
      range.x <- list(0,0)
      for (id in (1:2))
      {
         range.x[[id]] <- c(min(x[,id])-1.5*h[id],max(x[,id])+1.5*h[id])  
      }
   }

   a <- c(range.x[[1]][1],range.x[[2]][1])
   b <- c(range.x[[1]][2],range.x[[2]][2])

   # Set up grid points and bin the data
  
   gpoints1 <- seq(a[1],b[1],length=M[1])
   gpoints2 <- seq(a[2],b[2],length=M[2])

   gcounts <- linbin2D(x,gpoints1,gpoints2)

   # Compute kernel weights 

   L <- numeric()
   kapid <- list(0,0)
   for (id in (1:2))
   {
      L[id] <- min(floor(tau*h[id]*(M[id]-1)/(b[id]-a[id])),(M[id]-1))
      lvecid <- (0:L[id])
      facid <- (b[id]-a[id])/(h[id]*(M[id]-1))
      kapid[[id]] <- matrix(dnorm(lvecid*facid)/h[id])
   }
   kapp <- kapid[[1]]%*%(t(kapid[[2]]))/n  
      
   # Now combine weight and counts using the FFT
   # to obtain estimate

   P <- 2^(ceiling(log(M+L)/log(2))) # smallest powers of 2 >= M+L
   L1 <- L[1] ; L2 <- L[2] 
   M1 <- M[1] ; M2 <- M[2]               
   P1 <- P[1] ; P2 <- P[2]

   rp <- matrix(0,P1,P2)
   rp[1:(L1+1),1:(L2+1)] <- kapp
   if (L1>0) rp[(P1-L1+1):P1,1:(L2+1)] <- kapp[(L1+1):2,1:(L2+1)]
   if (L2>0) rp[,(P2-L2+1):P2] <- rp[,(L2+1):2]
                               # wrap-around version of "kapp" 

   sp <- matrix(0,P1,P2)
   sp[1:M1,1:M2] <- gcounts
                               # zero-padded version of "gcounts" 

   rp <- fft(rp)                       # Obtain FFT's of r and s
   sp <- fft(sp)
   rp <- Re(fft(rp*sp,inverse=T)/(P1*P2))[1:M1,1:M2]
                             # invert element-wise product of FFT's
                             # and truncate and normalise it

   # Ensure that rp is non-negative

   rp <- rp*matrix(as.numeric(rp>0),nrow(rp),ncol(rp))

   return(list(x1=gpoints1,x2=gpoints2,fhat=rp))
}


bkfe <- function(x,drv,bandwidth,gridsize=401,range.x,binned=F,truncate=T)

# Last changed: 18/12/95

{
   if (missing(range.x)&(binned==F)) range.x <- c(min(x),max(x))

   # Rename variables

   M <- gridsize
   a <- range.x[1]
   b <- range.x[2]
   h <- bandwidth

   # Bin the data if not already binned

   if (binned==F)
   {
      gpoints <- seq(a,b,length=gridsize)
      gcounts <- linbin(x,gpoints,truncate) 
   }
   else
   {
      gcounts <- x
      M <- length(gcounts)
      gpoints <- seq(a,b,length=M)
   }

   # Set the sample size and bin width 

   n <- sum(gcounts)
   delta <- (b-a)/(M-1)
   
   # Obtain kernel weights

   tau <- 4 + drv
   L <- min(floor(tau*h/delta),M) 
   lvec <- (0:L)
   arg <- lvec*delta/h 
 
   kappam <- dnorm(arg)/(h^(drv+1))
   hmold0 <- 1
   hmold1 <- arg
   hmnew <- 1
   if (drv >= 2) for (i in (2:drv)) 
   {
      hmnew <- arg*hmold1 - (i-1)*hmold0
      hmold0 <- hmold1   # Compute mth degree Hermite polynomial
      hmold1 <- hmnew    # by recurrence.
   }
   kappam <- hmnew*kappam

   # Now combine weights and counts to obtain estimate

   P <- 2^(ceiling(log(M+L)/log(2)))
   kappam <- c(kappam,rep(0,P-2*L-1),kappam[(L+1):2])
   Gcounts <- c(gcounts,rep(0,P-M))
   kappam <- fft(kappam)
   Gcounts <- fft(Gcounts)
   
   est <- sum(gcounts*(Re(fft(kappam*Gcounts,T))/P)[1:M])/(n^2)

   return(est)  
}
########## S-function: blkest ##########

# For obtaining preliminary estimates of
# quantities required for the "direct plug-in"
# regression bandwidth selector based on
# blocked qth degree polynomial fits.

# Last changed: 06/04/95

blkest <- function(x,y,Nval,q)
{

   n <- length(x)

   # Sort the (x,y) data with respect to
   # the x's.

   datmat <- cbind(x,y)
#   datmat <- datmat[sort.list(datmat[,1]),]
   datmat <- datmat[order(datmat[,1]),]
   x <- datmat[,1]
   y <- datmat[,2]

   # Set up arrays for FORTRAN programme "blkest"

   qq <- q + 1
   xj <- rep(0,n)
   yj <- rep(0,n)
   coef <- rep(0,qq)
   Xmat <- matrix(0,n,qq)
   wk <- rep(0,n)
   qraux <- rep(0,qq)
   sigsqe <- 0
   th22e <- 0
   th24e <- 0

   out <- .Fortran("blkest",as.double(x),as.double(y),as.integer(n),
                    as.integer(q),as.integer(qq),as.integer(Nval),as.double(xj),
                    as.double(yj),as.double(coef),as.double(Xmat),as.double(wk),
                    as.double(qraux),as.double(sigsqe),as.double(th22e),
                    as.double(th24e))
                
   return(list(sigsqe=out[[13]],th22e=out[[14]],th24e=out[[15]]))                           

}

######### End of S-function blkest ########
########## S-function: cpblock ##########

# Chooses the number of blocks for the premilinary
# step of a plug-in rule using Mallows' C_p.

# Last changed: 06/04/95

cpblock <- function(X,Y,Nmax,q)

{
   n <- length(X)

   # Sort the (X,Y) data with respect to
   # the X's.

   datmat <- cbind(X,Y)
#   datmat <- datmat[sort.list(datmat[,1]),]
   datmat <- datmat[order(datmat[,1]),]
   X <- datmat[,1]
   Y <- datmat[,2]

   # Set up arrays for FORTRAN subroutine "cp"

   qq <- q + 1
   RSS <- rep(0,Nmax)
   Xj <- rep(0,n)
   Yj <- rep(0,n)
   coef <- rep(0,qq)
   Xmat <- matrix(0,n,qq)
   Cpvals <- rep(0,Nmax)
   wk <- rep(0,n)
   qraux <- rep(0,qq)

   out <- .Fortran("cp",as.double(X),as.double(Y),as.integer(n),as.integer(q),
                   as.integer(qq),as.integer(Nmax),as.double(RSS),as.double(Xj),
                   as.double(Yj),as.double(coef),as.double(Xmat),as.double(wk),
                   as.double(qraux),as.double(Cpvals))
                
   Cpvec <- out[[14]]
   Nval <- order(Cpvec)[1]
 
   return(Nval)
}

######### End of S-function cpblock ########
dpih <- function(x,scalest="minim",level=2,gridsize=401,range.x=range(x),
                          truncate=T)

# Last changed: 16/06/95

{
   if (level>5) return("Level should be between 0 and 5")
  
   # Rename variables

   n <- length(x)
   M <- gridsize
   a <- range.x[1]
   b <- range.x[2]

   # Set up grid points and bin the data

   gpoints <- seq(a,b,length=M)
   gcounts <- linbin(x,gpoints,truncate)   
   delta <- (b-a)/(M-1)

   # Compute scale estimate

   if (scalest=="stdev") 
   {
      scalest <- sqrt(var(x))
   }
   else if (scalest=="iqr")
   {   
      scalest <- (quantile(x,3/4)-quantile(x,1/4))/1.349 
   }
   else if (scalest=="minim")
   {
      scalest <- (quantile(x,3/4)-quantile(x,1/4))/1.349
      scalest <- min(scalest,sqrt(var(x)))
   }
   
   # Perform plug-in steps 

   if (level==0)
   {
      hpi <- (24*sqrt(pi)/n)^(1/3)*scalest
   }
   else if (level==1)
   {
      alpha <- (2/(3*n))^(1/5)*sqrt(2)*scalest   # bandwidth for psi_2
      psi2hat <- bkfe(gcounts,2,alpha,range.x=c(a,b),binned=T)
      hpi <- (6/(-psi2hat*n))^(1/3)
   }
   else if (level==2)
   {
      alpha <- ((2/(5*n))^(1/7))*sqrt(2)*scalest    # bandwidth for psi_4
      psi4hat <- bkfe(gcounts,4,alpha,range.x=c(a,b),binned=T)
      alpha <- (sqrt(2/pi)/(psi4hat*n))^(1/5)       # bandwidth for psi_2
      psi2hat <- bkfe(gcounts,2,alpha,range.x=c(a,b),binned=T)
      hpi <- (6/(-psi2hat*n))^(1/3)
   }   
   else if (level==3)
   {
      alpha <- ((2/(7*n))^(1/9))*sqrt(2)*scalest    # bandwidth for psi_6
      psi6hat <- bkfe(gcounts,6,alpha,range.x=c(a,b),binned=T)
      alpha <- (-3*sqrt(2/pi)/(psi6hat*n))^(1/7)    # bandwidth for psi_4       
      psi4hat <- bkfe(gcounts,4,alpha,range.x=c(a,b),binned=T)
      alpha <- (sqrt(2/pi)/(psi4hat*n))^(1/5)       # bandwidth for psi_2
      psi2hat <- bkfe(gcounts,2,alpha,range.x=c(a,b),binned=T)
      hpi <- (6/(-psi2hat*n))^(1/3)
   }   
   else if (level==4)
   {
      alpha <- ((2/(9*n))^(1/11))*sqrt(2)*scalest    # bandwidth for psi_8
      psi8hat <- bkfe(gcounts,8,alpha,range.x=c(a,b),binned=T)
      alpha <- (15*sqrt(2/pi)/(psi8hat*n))^(1/9)     # bandwidth for psi_6
      psi6hat <- bkfe(gcounts,6,alpha,range.x=c(a,b),binned=T)
      alpha <- (-3*sqrt(2/pi)/(psi6hat*n))^(1/7)     # bandwidth for psi_4
      psi4hat <- bkfe(gcounts,4,alpha,range.x=c(a,b),binned=T)
      alpha <- (sqrt(2/pi)/(psi4hat*n))^(1/5)        # bandwidth for psi_2
      psi2hat <- bkfe(gcounts,2,alpha,range.x=c(a,b),binned=T)
      hpi <- (6/(-psi2hat*n))^(1/3)
   }   
   else if (level==5)
   {
      alpha <- ((2/(11*n))^(1/13))*sqrt(2)*scalest
      psi10hat <- bkfe(gcounts,10,alpha,range.x=c(a,b),binned=T)
      alpha <- (-105*sqrt(2/pi)/(psi10hat*n))^(1/11) # bandwidth for psi_8
      psi8hat <- bkfe(gcounts,8,alpha,range.x=c(a,b),binned=T)
      alpha <- (15*sqrt(2/pi)/(psi8hat*n))^(1/9)     # bandwidth for psi_6
      psi6hat <- bkfe(gcounts,6,alpha,range.x=c(a,b),binned=T)
      alpha <- (-3*sqrt(2/pi)/(psi6hat*n))^(1/7)     # bandwidth for psi_4
      psi4hat <- bkfe(gcounts,4,alpha,range.x=c(a,b),binned=T)
      alpha <- (sqrt(2/pi)/(psi4hat*n))^(1/5)        # bandwidth for psi_2
      psi2hat <- bkfe(gcounts,2,alpha,range.x=c(a,b),binned=T)
      hpi <- (6/(-psi2hat*n))^(1/3)
   }   

   return(hpi)  
}
dpik <- function(x,scalest="minim",level=2,kernel="normal",
                canonical=F,gridsize=401,range.x=range(x),truncate=T)

# Last changed: 18/06/97

{
   if (level>5) return("Level should be between 0 and 5")
   
   # Set kernel constants
   
   if (canonical) del0 <- 1
   else
   {
      if (kernel=="normal") del0 <- 1/((4*pi)^(1/10))
      if (kernel=="box") del0 <- (9/2)^(1/5)
      if (kernel=="epanech") del0 <- 15^(1/5)
      if (kernel=="biweight") del0 <- 35^(1/5)
      if (kernel=="triweight") del0 <- (9450/143)^(1/5)
   }

   # Rename variables

   n <- length(x)
   M <- gridsize
   a <- range.x[1]
   b <- range.x[2]

   # Set up grid points and bin the data

   gpoints <- seq(a,b,length=M)
   gcounts <- linbin(x,gpoints,truncate)   
   delta <- (b-a)/(M-1)
 
   # Compute scale estimate

   if (scalest=="stdev") 
   {
      scalest <- sqrt(var(x))
   }
   else if (scalest=="iqr")
   {   
      scalest <- (quantile(x,3/4)-quantile(x,1/4))/1.349 
   }
   else if (scalest=="minim")
   {
      scalest <- (quantile(x,3/4)-quantile(x,1/4))/1.349
      scalest <- min(scalest,sqrt(var(x)))
   }
   
   # Perform plug-in steps 

   if (level==0)
   {
      psi4hat <- 3/(8*sqrt(pi)*scalest^5)
   }
   else if (level==1)
   {
      alpha <- (2*(sqrt(2)*scalest)^7/(5*n))^(1/7)   # bandwidth for psi_4
      psi4hat <- bkfe(gcounts,4,alpha,range.x=c(a,b),binned=T)
   }
   else if (level==2)
   {
      alpha <- (2*(sqrt(2)*scalest)^9/(7*n))^(1/9)   # bandwidth for psi_6
      psi6hat <- bkfe(gcounts,6,alpha,range.x=c(a,b),binned=T)
      alpha <- (-3*sqrt(2/pi)/(psi6hat*n))^(1/7)   # bandwidth for psi_4 
      psi4hat <- bkfe(gcounts,4,alpha,range.x=c(a,b),binned=T)
   }   
   else if (level==3)   
   {
      alpha <- (2*(sqrt(2)*scalest)^11/(9*n))^(1/11) # bandwidth for psi_8
      psi8hat <- bkfe(gcounts,8,alpha,range.x=c(a,b),binned=T)
      alpha <- (15*sqrt(2/pi)/(psi8hat*n))^(1/9)   # bandwidth for psi_6 
      psi6hat <- bkfe(gcounts,6,alpha,range.x=c(a,b),binned=T)
      alpha <- (-3*sqrt(2/pi)/(psi6hat*n))^(1/7)   # bandwidth for psi_4 
      psi4hat <- bkfe(gcounts,4,alpha,range.x=c(a,b),binned=T)
   }
   else if (level==4)
   {
      alpha <- (2*(sqrt(2)*scalest)^13/(11*n))^(1/13) # bandwidth for psi_10
      psi10hat <- bkfe(gcounts,10,alpha,range.x=c(a,b),binned=T)
      alpha <- (-105*sqrt(2/pi)/(psi10hat*n))^(1/11)# bandwidth for psi_8 
      psi8hat <- bkfe(gcounts,8,alpha,range.x=c(a,b),binned=T)
      alpha <- (15*sqrt(2/pi)/(psi8hat*n))^(1/9)    # bandwidth for psi_6 
      psi6hat <- bkfe(gcounts,6,alpha,range.x=c(a,b),binned=T)
      alpha <- (-3*sqrt(2/pi)/(psi6hat*n))^(1/7)    # bandwidth for psi_4 
      psi4hat <- bkfe(gcounts,4,alpha,range.x=c(a,b),binned=T)
   }      
   else if (level==5)
   {
      alpha <- (2*(sqrt(2)*scalest)^15/(13*n))^(1/15) # bandwidth for psi_12
      psi12hat <- bkfe(gcounts,12,alpha,range.x=c(a,b),binned=T)
      alpha <- (945*sqrt(2/pi)/(psi12hat*n))^(1/13) # bandwidth for psi_10 
      psi10hat <- bkfe(gcounts,10,alpha,range.x=c(a,b),binned=T)
      alpha <- (-105*sqrt(2/pi)/(psi10hat*n))^(1/11)# bandwidth for psi_8 
      psi8hat <- bkfe(gcounts,8,alpha,range.x=c(a,b),binned=T)
      alpha <- (15*sqrt(2/pi)/(psi8hat*n))^(1/9)    # bandwidth for psi_6 
      psi6hat <- bkfe(gcounts,6,alpha,range.x=c(a,b),binned=T)
      alpha <- (-3*sqrt(2/pi)/(psi6hat*n))^(1/7)    # bandwidth for psi_4 
      psi4hat <- bkfe(gcounts,4,alpha,range.x=c(a,b),binned=T)
   }

   hpi <- del0*(1/(psi4hat*n))^(1/5)
   
   return(hpi)  
}
########## S-function: dpill ##########

# Computes a direct plug-in selector of the
# bandwidth for local linear regression as
# described in the 1996 J. Amer. Statist. Assoc. 
# paper by Ruppert, Sheather and Wand.

# Last changed: 16/06/95 

dpill <- function(x,y,blockmax=5,divisor=20,trim=0.01,proptrun=0.05,gridsize=401,
                  range.x=range(x),truncate=T)

{
   # Trim the 100(trim)% of the data from each end (in the x-direction).

   xy <- cbind(x,y)
#   xy <- xy[sort.list(xy[,1]),]
   xy <- xy[order(xy[,1]),]
   x <- xy[,1]
   y <- xy[,2]
   indlow <- floor(trim*length(x)) + 1
   indupp <- length(x) - floor(trim*length(x))  

   x <- x[indlow:indupp]
   y <- y[indlow:indupp]

   # Rename common parameters
   
   n <- length(x)
   M <- gridsize
   a <- range.x[1]
   b <- range.x[2]
   delta <- (b-a)/(M-1)

   # Bin the data 

   gpoints <- seq(a,b,length=M)
   out <- rlbin(x,y,gpoints,truncate)
   xcounts <- out$xcounts
   ycounts <- out$ycounts

   # Choose the value of N using Mallow's C_p

   Nmax <- max(min(floor(n/divisor),blockmax),1)

   Nval <- cpblock(x,y,Nmax,4)

   # Estimate sig^2, theta_22 and theta_24 using quartic fits 
   # on "Nval" blocks.

   out <- blkest(x,y,Nval,4)
   sigsqQ <- out$sigsqe
   th22Q <- out$th22e
   th24Q <- out$th24e

   # Estimate theta_22 using a local cubic fit
   # with a "rule-of-thumb" bandwidth: "gamseh"

   gamseh <- (sigsqQ*(b-a)/(abs(th24Q)*n))
   if (th24Q<0) gamseh <- (3*gamseh/(8*sqrt(pi)))^(1/7)
   if (th24Q>0) gamseh <- (15*gamseh/(16*sqrt(pi)))^(1/7)

   mddest <- locpoly(xcounts,ycounts,drv=2,bandwidth=gamseh,
                     range.x=range.x,binned=T)$y

   llow <- floor(proptrun*M) + 1
   lupp <- M - floor(proptrun*M)
   th22kn <- sum((mddest[llow:lupp]^2)*xcounts[llow:lupp])/n

   # Estimate sigma^2 using a local linear fit
   # with a "direct plug-in" bandwidth: "lamseh"

   C3K <- (1/2) + 2*sqrt(2) - (4/3)*sqrt(3)
   C3K <- (4*C3K/(sqrt(2*pi)))^(1/9)
   lamseh <- C3K*(((sigsqQ^2)*(b-a)/((th22kn*n)^2))^(1/9))

   # Now compute a local linear kernel estimate of
   # the variance.

   mest <- locpoly(xcounts,ycounts,bandwidth=lamseh,
                   range.x=range.x,binned=T)$y

   Sdg <- sdiag(xcounts,bandwidth=lamseh,
                range.x=range.x,binned=T)$y

   SSTdg <- sstdiag(xcounts,bandwidth=lamseh,
                    range.x=range.x,binned=T)$y

   sigsqn <- sum(y^2) - 2*sum(mest*ycounts) + sum((mest^2)*xcounts)

   sigsqd <- n - 2*sum(Sdg*xcounts) + sum(SSTdg*xcounts)

   sigsqkn <- sigsqn/sigsqd
  
   # Combine to obtain final answer.

   hhat <- (sigsqkn*(b-a)/(2*sqrt(pi)*th22kn*n))^(1/5)
 
   return(hhat)  
}

######### End of S-function dpill.S ########
######## S-function: linbin ########

# For application of linear binning to a univariate 
# data set.

# Last changed: 16/06/95

linbin <- function(X,gpoints,truncate=T)

{
   n <- length(X)
   M <- length(gpoints)  
   trun <- 0
   if (truncate) trun <- 1
   a <- gpoints[1]
   b <- gpoints[M]
   out <- .Fortran("linbin",as.double(X),as.integer(n),
           as.double(a),as.double(b),as.integer(M),
           as.integer(trun),double(M))
   return(out[[7]])
}

########## End of S-function linbin ##########
######### S-function: linbin2D #########
 
# Creates the grid counts from a bivariate data set X 
# over an equally-spaced set of grid points
# contained in "gpoints" using the linear 
# binning strategy. Note that the FORTRAN subroutine
# "lbtwod" is called. 

# Last changed: 25/08/95

linbin2D <- function(X,gpoints1,gpoints2)
{
   n <- nrow(X)
   X <- c(X[,1],X[,2]) 
   M1 <- length(gpoints1)
   M2 <- length(gpoints2)
   a1 <- gpoints1[1]
   a2 <- gpoints2[1]
   b1 <- gpoints1[M1]
   b2 <- gpoints2[M2]
   out <- .Fortran("lbtwod",as.double(X),as.integer(n),
           as.double(a1),as.double(a2),as.double(b1),as.double(b2),
           as.integer(M1),as.integer(M2),double(M1*M2))
   return(matrix(out[[9]],M1,M2))
}

######## End of S-function linbin2D ########
######## S-function: locpoly ########

# For computing a binned local polynomial 
# regression estimator of a univariate regression
# function or its derivative.
# The data are discretised on an equally
# spaced grid. The bandwidths are discretised on a 
# logarithmically spaced grid.

# Last changed: 07/07/95 
 
locpoly <- function(x,y,drv=0,degree,kernel="normal",
                    bandwidth,gridsize=401,bwdisc=25,range.x,
                    binned=F,truncate=T)

{  
   if (missing(degree)) degree <- drv + 1

   if (missing(range.x)&(binned==F)) 
   if (missing(y))
   {
      extra <- 0.05*(max(x) - min(x))
      range.x <- c(min(x)-extra,max(x)+extra)
   }
   else
   {
      range.x <- c(min(x),max(x))
   }

   # Rename common variables

   M <- gridsize
   Q <- bwdisc
   a <- range.x[1]
   b <- range.x[2]
   pp <- degree + 1
   ppp <- 2*degree + 1
   tau <- 4

   # Decide whether a density estimate or regression
   # estimate is required.

   if (missing(y))     # obtain density estimate
   {
      y <- NULL
      n <- length(x)
      gpoints <- seq(a,b,length=M)    
      xcounts <- linbin(x,gpoints,truncate)
      ycounts <- (M-1)*xcounts/(n*(b-a))
      xcounts <- rep(1,M)
   }
   else                # obtain regression estimate
   {
      
      # Bin the data if not already binned

      if (binned==F)
      {
         gpoints <- seq(a,b,length=M) 
         out <- rlbin(x,y,gpoints,truncate)
         xcounts <- out$xcounts        
         ycounts <- out$ycounts
      }
      else
      {
         xcounts <- x
         ycounts <- y
         M <- length(xcounts) 
         gpoints <- seq(a,b,length=M)
      }
   }

   # Set the bin width

   delta <- (b-a)/(M-1)

   # Discretise the bandwidths

   if (length(bandwidth)==M)
   {

      hlow <- sort(bandwidth)[1]
      hupp <- sort(bandwidth)[M]

      hdisc <- exp(seq(log(hlow),log(hupp),length=Q))
  
      # Determine value of L for each member of "hdisc"
      
      Lvec <- floor(tau*hdisc/delta)

      # Determine index of closest entry of "hdisc" 
      # to each member of "bandwidth"

      if (Q > 1)
      {
         lhdisc <- log(hdisc)
         gap <- (lhdisc[Q]-lhdisc[1])/(Q-1)
         if (gap==0)
            indic <- rep(1,M)
         else
         {
            tlhvec <- ((log(bandwidth) - log(sort(bandwidth)[1]))/gap) + 1
            indic <- round(tlhvec)
         }
      }
      else indic <- rep(1,M)
   }
   else if (length(bandwidth)==1)
   {
      indic <- rep(1,M)
      Q <- 1 
      Lvec <- rep(floor(tau*bandwidth/delta),Q)
      hdisc <- rep(bandwidth,Q)     
   }
   else 
   {
      print("Bandwidth must be a scalar or an array of length gridsize")
      return()
   }

   # Allocate space for the kernel vector and final estimate

   dimfkap <- 2*sum(Lvec) + Q  
   fkap <- rep(0,dimfkap)
   curvest <- rep(0,M)
   midpts <- rep(0,Q)
   ss <- matrix(0,M,ppp)
   tt <- matrix(0,M,pp)
   Smat <- matrix(0,pp,pp)
   Tvec <- rep(0,pp)
   ipvt <- rep(0,pp)

   # Call FORTRAN routine "locpol"

   out <- .Fortran("locpol",as.double(xcounts),as.double(ycounts),
                  as.integer(drv),as.double(delta),as.double(hdisc),
                  as.integer(Lvec),as.integer(indic),as.integer(midpts),
                  as.integer(M),as.integer(Q),as.double(fkap),as.integer(pp),
                  as.integer(ppp),as.double(ss),as.double(tt),
                  as.double(Smat),as.double(Tvec),as.integer(ipvt),
                  as.double(curvest))   
   
   curvest <- gamma(drv+1)*out[[19]]

   return(list(x=gpoints,y=curvest))
}

########## End of S-function locpoly ##########
######## S-function: rlbin ########

# For application of linear binning to a regression
# data set.

# Last changed: 16/06/95 

rlbin <- function(X,Y,gpoints,truncate=T)

{
   n <- length(X)
   M <- length(gpoints)  
   if (truncate) trun <- 1
   else trun <- 0
   a <- gpoints[1]
   b <- gpoints[M]
   out <- .Fortran("rlbin",as.double(X),as.double(Y),as.integer(n),
           as.double(a),as.double(b),as.integer(M),as.integer(trun),
           double(M),double(M))
   return(list(xcounts=out[[8]],ycounts=out[[9]]))
}

########## End of S-function rlbin ##########
########## S-function: sdiag ##########

# For computing the binned diagonal entries of a smoother
# matrix for local polynomial kernel regression.

# Last changed: 16/06/95


sdiag <- function(x,drv=0,degree=1,kernel="normal",
                    bandwidth,gridsize=401,bwdisc=25,range.x,
                    binned=F,truncate=T)
 
{

   if (missing(range.x)&(binned==F)) range.x <- c(min(x),max(x))

   # Rename common variables
   
   M <- gridsize
   Q <- bwdisc
   a <- range.x[1]  
   b <- range.x[2]
   pp <- degree + 1
   ppp <- 2*degree + 1
   tau <- 4

   # Bin the data if not already binned

   if (binned==F)
   {
      gpoints <- seq(a,b,length=M) 
      xcounts <- linbin(x,gpoints,truncate)
   }
   else
   {
      xcounts <- x
      M <- length(xcounts) 
      gpoints <- seq(a,b,length=M)
   }

   # Set the bin width

   delta <- (b-a)/(M-1)

   # Discretise the bandwidths

   if (length(bandwidth)==M)
   {

      hlow <- sort(bandwidth)[1]
      hupp <- sort(bandwidth)[M]

      hdisc <- exp(seq(log(hlow),log(hupp),length=Q))
  
      # Determine value of L for each member of "hdisc"
      
      Lvec <- floor(tau*hdisc/delta)

      # Determine index of closest entry of "hdisc" 
      # to each member of "bandwidth"

      if (Q > 1)
      {
         lhdisc <- log(hdisc)
         gap <- (lhdisc[Q]-lhdisc[1])/(Q-1)
         if (gap==0)
            indic <- rep(1,M)
         else
         {
            tlhvec <- ((log(bandwidth) - log(sort(bandwidth)[1]))/gap) + 1
            indic <- round(tlhvec)
         }
      }
      else indic <- rep(1,M)
   }
   else if (length(bandwidth)==1)
   {
      indic <- rep(1,M)
      Q <- 1 
      Lvec <- rep(floor(tau*bandwidth/delta),Q)
      hdisc <- rep(bandwidth,Q)     
   }
   else 
   {
      print("Bandwidth must be a scalar or an array of length gridsize")
      return()
   }
 
   dimfkap <- 2*sum(Lvec) + Q
   fkap <- rep(0,dimfkap)
   midpts <- rep(0,Q)
   ss <- matrix(0,M,ppp)
   Smat <- matrix(0,pp,pp)
   work <- rep(0,pp)
   det <- rep(0,2)
   ipvt <- rep(0,pp)
   Sdg <- rep(0,M)

   out <- .Fortran("sdiag",as.double(xcounts),as.double(delta),
                    as.double(hdisc),as.integer(Lvec),as.integer(indic),
                    as.integer(midpts),as.integer(M),as.integer(Q),
                    as.double(fkap),as.integer(pp),as.integer(ppp),
                    as.double(ss),as.double(Smat),as.double(work),
                    as.double(det),as.integer(ipvt),as.double(Sdg))

   Sdg <- out[[17]]

   return(list(x=gpoints,y=Sdg))
}

######## End of S-function sdiag ########

########## S-function: sstdiag ##########

# For computing the binned diagonal entries of SS^T
# where S is a smoother matrix for local polynomial 
# kernel regression.

# Last changed: 16/06/95

sstdiag <- function(x,drv=0,degree=1,kernel="normal",
                    bandwidth,gridsize=401,bwdisc=25,range.x,
                    binned=F,truncate=T)
 
{

   if (missing(range.x)&(binned==F)) range.x <- c(min(x),max(x))

   # Rename common variables
   
   M <- gridsize
   Q <- bwdisc
   a <- range.x[1]  
   b <- range.x[2]
   pp <- degree + 1
   ppp <- 2*degree + 1
   tau <- 4

   # Bin the data if not already binned

   if (binned==F)
   {
      gpoints <- seq(a,b,length=M) 
      xcounts <- linbin(x,gpoints,truncate)
   }
   else
   {
      xcounts <- x
      M <- length(xcounts) 
      gpoints <- seq(a,b,length=M)
   }

   # Set the bin width

   delta <- (b-a)/(M-1)

   # Discretise the bandwidths

   if (length(bandwidth)==M)
   {

      hlow <- sort(bandwidth)[1]
      hupp <- sort(bandwidth)[M]

      hdisc <- exp(seq(log(hlow),log(hupp),length=Q))
  
      # Determine value of L for each member of "hdisc"
      
      Lvec <- floor(tau*hdisc/delta)

      # Determine index of closest entry of "hdisc" 
      # to each member of "bandwidth"

      if (Q > 1)
      {
         lhdisc <- log(hdisc)
         gap <- (lhdisc[Q]-lhdisc[1])/(Q-1)
         if (gap==0)
            indic <- rep(1,M)
         else
         {
            tlhvec <- ((log(bandwidth) - log(sort(bandwidth)[1]))/gap) + 1
            indic <- round(tlhvec)
         }
      }
      else indic <- rep(1,M)
   }
   else if (length(bandwidth)==1)
   {
      indic <- rep(1,M)
      Q <- 1 
      Lvec <- rep(floor(tau*bandwidth/delta),Q)
      hdisc <- rep(bandwidth,Q)     
   }
   else 
   {
      print("Bandwidth must be a scalar or an array of length gridsize")
      return()
   }
 
   dimfkap <- 2*sum(Lvec) + Q
   fkap <- rep(0,dimfkap)
   midpts <- rep(0,Q)
   ss <- matrix(0,M,ppp)
   uu <- matrix(0,M,ppp)
   Smat <- matrix(0,pp,pp)
   Umat <- matrix(0,pp,pp)
   work <- rep(0,pp)
   det <- rep(0,2)
   ipvt <- rep(0,pp)
   SSTd <- rep(0,M)

   out <- .Fortran("sstdg",as.double(xcounts),as.double(delta),
                    as.double(hdisc),as.integer(Lvec),as.integer(indic),
                    as.integer(midpts),as.integer(M),as.integer(Q),
                    as.double(fkap),as.integer(pp),as.integer(ppp),
                    as.double(ss),as.double(uu),as.double(Smat),
                    as.double(Umat),as.double(work),as.double(det),
                    as.integer(ipvt),as.double(SSTd))

   SSTd <- out[[19]]

   return(list(x=gpoints,y=SSTd))
}

######## End of S-function sstdiag ########

######### S-PLUS function: .First.lib ##########

# S-PLUS .First.lib function for ModuleX

.First.lib <- function(lib, pkg)
{
   library.dynam("KernSmooth", pkg, lib)
   cat("KernSmooth 2.22 installed\n")
   cat("Copyright M. P. Wand 1997\n")
}

############# End of .First.lib ###############
