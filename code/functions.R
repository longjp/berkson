###
###
### FUNCTIONS FOR GENERATING RESULTS PRESENTED IN
### "FINITE SAMPLE RESULTS" SECTION OF
### "Kernel Density Estimation with Berkson Error"
### by Long, El Karoui, and Rice
###
### by: James Long
### date: July 12, 2013
###

## IMPLEMENTATION NOTE: ComputeOmega, ComputeMISE, and
## OptimizeMISE work with a scalar bandwidth
## and in arbitrary dimension p. ComputeIntF2, ComputeW2dv,
## ComputeW4du, and ComputeOptimalHfy work for a
## scalar bandwidth but in only p=1 dimension.

library('ks')
library(mvtnorm) ## added December 8,2020 because new version of ks does not import this

ComputeHy <- function(n,v.e,v.x){
    ## compute asymptotically optimal smoothing for
    ## a normal density and normal error
    return(sqrt((4 / (3*n))*((v.e + v.x)^(5/2) / v.e^(3/2) - (v.e + v.x))))
}

ComputeOmega <- function(a,S,evar,means,mvar){
  ## computes \Omega_2, \Omega_1, and \Omega_0
  ## see "Finite Sample Results" section in "Kernel Density with Berkson
  ## Error" by Long, El Karoui, and Rice for precise definitions
  ##
  ## used by ComputeMISE
  ##
  ## arguments:
  ##     a : 0,1,2 corresponding to \Omega_0, \Omega_1, \Omega_2
  ##     S : the squared bandwidth times variance of kernel.
  ##         must be a non-negative scalar in current implementation
  ##  evar : p x p covariance matrix for error (covariance is
  ##         assumed normal, mean 0)
  ## means : p x k vector specifying means of components of multivariate
  ##         gaussian density (p = dimension, k = # of comp.)
  ##  mvar : p x p x k array of covariances of components
  ##         (p = dimension, k = # of comp.)
  ##
  ##     value
  ##  \Omega_a
  p <- dim(means)[1]
  k <- dim(means)[2]
  omega <- matrix(0,nrow=k,ncol=k)
  ## omega is symmetric so only compute upper triangle
  for(ii in 1:k){
    for(jj in ii:k){
      omega[ii,jj] <- dmvnorm(rep(0,p),
                              means[,ii] - means[,jj],
                              (a*diag(S,p) + 2*evar + mvar[,,ii] + mvar[,,jj]))
    }
  }
  ## replace lower triangle 0s with upper triangle values
  omega <- omega + t(omega)
  diag(omega) <- diag(omega)/2
  return(omega)
}


ComputeMISE <- function(S,n,evar,den){
  ## computes mise when f_X is mixture of gaussians
  ## and kernel and error are normal
  ##
  ## arguments:
  ##      S : covariance of kernel
  ##      n : sample size
  ##   evar : p x p matrix specifying variance of error
  ##    den : list specifying f_X (see ConstructX functions
  ##          for details)
  ##
  ##   value
  ##   the exact MISE
  p <- dim(evar)[1]
  term1 <- (1/n)*dmvnorm(rep(0,p),rep(0,p),2*diag(S,p)+2*evar)
  omega2 <- ComputeOmega(2,S,evar,den$means,den$mvar)
  omega1 <- ComputeOmega(1,S,evar,den$means,den$mvar)
  omega0 <- ComputeOmega(0,S,evar,den$means,den$mvar)
  return((term1 + t(den$ws) %*% ((1 - 1/n)*omega2 -
                                 2*omega1 + omega0)%*%den$ws))
}



OptimizeMISE <- function(n,evar,den){
  ## returns the squared bandwidth (S) that minimizes the MISE for
  ## a particular sample size n, error variance evar, and density den
  ## NOTE: S (i.e. H) must be a scalar
  ## 
  ## arguments:
  ##    n : sample size
  ## evar : covariance of error
  ##  den : list specifying f_X density (see ConstructX functions for
  ##        details and formatting of den)
  ##
  ## value
  ##  opt : list where MISE(opt$minimum) = opt$objective
  FrozenComputeMISE <- function(S){   # construct a 1-d computeMISE function
    return(ComputeMISE(S,n,evar,den))
  }
  opt <- optimize(FrozenComputeMISE,lower=0,upper=10,
                  tol = .Machine$double.eps)
  return(opt)
}




ComputeIntF2 <- function(den,deriv=0){
  ## \int (f^{(s)}(x))^2 dx for s = 1,2
  ## where f is a 1-d gaussian mixture
  ##
  ## see Theorem 4.1 in "Exact Mean Integrated Squared Error"
  ## by Marron and Wand for formula implemented here
  ##
  ##   arguments:
  ##      den : list specifying the gaussian mixture (see
  ##               ConstructX functions for formatting of list)
  ##    deriv : 1,2 (for value of s)
  ##
  ## value
  ##   \int (f^{(s)}(x))^2 dx for s = 1,2 where s = deriv
  k <- length(den$ws)
  ws <- matrix(den$ws,nrow=k,ncol=1)
  mvar <- den$mvar
  means <- den$means
  M <- matrix(0,nrow=k,ncol=k)
  ## M is symmetric, so only compute upper triangle
  if (deriv==1){
    for(ii in 1:k){
      for(jj in ii:k){
        vr <- mvar[1,1,jj] + mvar[1,1,ii]
        mn <- means[1,jj] - means[1,ii]
        ## second derivative of mean=mn,var=vr normal evaluated
        M[ii,jj] <- (vr^(-1) - vr^(-2)*mn^2)*dnorm(0,mn,sqrt(vr))
      }
    }
  }
  if (deriv==2){
    for(ii in 1:k){
      for(jj in ii:k){
        vr <- mvar[1,1,jj] + mvar[1,1,ii]
        mn <- means[1,jj] - means[1,ii]
        ## fourth derivative of mean=mn,var=vr normal evaluated
        M[ii,jj] <- (((-vr^(-1) + vr^(-2)*mn^2)^2 + 2/vr^2
                      - 4*mn^2/vr^3)*dnorm(0,mn,sqrt(vr)))
      }
    }
  }
  ## replace lower triangle 0s with upper triangle values
  M <- M + t(M)
  diag(M) <- diag(M)/2
  return(t(ws)%*%M%*%ws)
}


ComputeW2dv <- function(den,evar){
  ## returns \int \omega^2 d\nu(\omega)
  ## when f_X is 1-d gaussian mixture and
  ## f_\epsilon is gaussian
  ##
  ## arguments:
  ##     den : specification of the gaussian mixture
  ##           see ConstructX functions for precise
  ##           specifications
  ##    evar : variance of error
  ##
  den1 <- list()
  den1$means <- matrix(c(0),nrow=1,ncol=1)
  den1$mvar <- array(evar,c(dim(evar),1))
  den1$ws <- 1
  w2fe <- 2*pi*ComputeIntF2(den1,deriv=1)
  ## now add variance of error to variances of
  ## each component in Gaussian mixture
  evar <- array(evar,c(dim(evar),dim(den$mvar)[3]))
  den$mvar <- den$mvar + evar    
  w2fefx <- 2*pi*ComputeIntF2(den,deriv=1)
  return(w2fe - w2fefx)
}

ComputeW4du <- function(den,evar){
  ## returns \int \omega^4 d\mu(\omega)
  ## when f_X is 1-d gaussian mixture and
  ## f_\epsilon is gaussian
  ##
  ## arguments:
  ##     den : specification of the gaussian mixture
  ##    evar : variance of error
  ##
  evar <- array(evar,c(dim(evar),dim(den$mvar)[3]))
  den$mvar <- den$mvar + evar    
  return(2*pi*ComputeIntF2(den,deriv=2))
}
    

ComputeOptimalHfy <- function(den,n,evar){
  ## computes asymptotically optimal h, h^* from
  ## Section 2.3.2
  ##
  ## only does 1-d case (limited by ComputeIntF2)
  ##
  ## arguments:
  ##    den : list specifying density (see constructX functions for details)
  ##      n : integer specifying sample size
  ##   evar : variance of error term
  w2dv <- ComputeW2dv(den,evar)
  w4du <- ComputeW4du(den,evar)
  return(sqrt((2*w2dv) / (n * w4du)))
}






####
#### THE FUNCTIONS BELOW ARE PRIMARILY FOR
#### GENERATING PLOTS AND TABLES FOR RESULTS
####




PlotDen <- function(den,main=TRUE,add=FALSE,lwd=3.5){
  ## plot normal mixture density (assumes den is 1-d)
  ##
  ## argument
  ##    den : list specifying density (see constructX functions for details)
  ##   main : title the plot den$main
  ##    add : should curve be added

  ## find 3 s.d. less than min component mean
  from <- min(den$means[1,] - 3*sqrt(den$mvar[1,1,]))
  ## find 3 s.d. greater than max component mean
  to <- max(den$means[1,] + 3*sqrt(den$mvar[1,1,]))
  ## vectorized evaluation of mixture
  Func <- function(x){
    return(sapply(x,function(x){
      sum(den$ws*dnorm(x,mean=den$means[1,],
                       sd=sqrt(den$mvar[1,1,])))}))
  }
  ## plot
  if (main){
    main <- den$name
  }
  else {
    main <- NULL
  }
  curve(Func,from=from,to=to,lwd=lwd,
        xlab="X",ylab="Density",main=main,
        add=add)
}


GenerateData <- function(den,n){
  ## generates data from a normal mixture
  ##
  ##  argument
  ##   den : density (see ConstructX functions for details)
  ##     n : number of samples
  ##
  ##  value
  ## draws : in dim. 1, vector of length n, in dim p,
  ##         matrix with dim(draws) = c(n,p)
  ##
  ## because rmvnorm.mixt only does p > 1 case and
  ##         rnorm.mixt only does p = 1 case,
  ##         have to treat separately
  if (dim(den$means)[1] > 1){
    p <- dim(den$means)[1]
    k <- dim(den$means)[2]
    mus <- t(den$means)
    ## here Sigmas is stacked variances of components
    Sigmas <- matrix(den$mvar,nrow=p*k,ncol=p,byrow=TRUE)
    draws <- rmvnorm.mixt(n=n,mus=mus,Sigmas=Sigmas,props=den$ws)
  }
  else {
    mus <- den$means[1,]
    ## here sigmas is s.d.
    sigmas <- sqrt(den$mvar[1,1,])
    draws <- rnorm.mixt(n=n,mus=mus,sigmas=sigmas,props=den$ws)
  }
  return(draws)
}


MakeRelativeErrorTable <- function(dens,n,evars){
  ## Compare MISEs using different amounts of smoothing for several
  ## densities and several error variances
  ##
  ## arguments:
  ##    dens : list of densities
  ##       n : sample size
  ##   evars : list of variances for error term
  ##
  ## value: returns character matrix with relative errors of different
  ##        methods
  results <- matrix(" ",nrow=length(evars),ncol=length(dens)+1)
  results[,1] <- sapply(evars,function(x){x[1,1]})
  for(ii in 1:length(dens)){
    den <- dens[[ii]]
    means <- den$means
    mvar <- den$mvar
    ws <- den$ws
    ## find mise with h=0
    mise0 <- sapply(evars,
                    function(evar){
                      ComputeMISE(0,n,evar,den)})
    ## find mise with h=hx
    p <- dim(means)[1]
    evar <- matrix(0,nrow=p,ncol=p) # set error variance to 0
    sx <- OptimizeMISE(n,evar,den)$minimum
    misex <- sapply(evars, 
                    function(evar){
                      ComputeMISE(sx,n,evar,den)})
    ## find mise with h=hy (best procedure)
    misey <- sapply(evars,
                    function(evar){
                      OptimizeMISE(n,evar,den)$objective})
    ## store ratios of errors in results
    results[,ii+1] <- paste("(",format(mise0 / misey,
                                       digits=2,nsmall=2),",",
                            format(misex / misey,
                                   digits=2,nsmall=2),
                            ")",sep="")
  }
  colnames(results) <- c("$\\sigma_{\\epsilon}^2$",
                         lapply(dens,function(x){x$name}))
  return(results)
}

MakeRelativeErrorTableEst <- function(dens,n,evars,N=500){
  ## Compare MISEs using different amounts of smoothing for several
  ## densities and several error variances. this function
  ## uses estimated smoothing parameters (sheather jones for sx)
  ## and min(rule of thumb,sheather jones) for sy
  ##
  ## arguments:
  ##    dens : list of densities
  ##       n : sample size
  ##   evars : list of variances for error term
  ##
  ## value: returns character matrix with relative errors of different
  ##        methods
  results <- matrix(" ",nrow=length(evars),ncol=length(dens)+1)
  results[,1] <- sapply(evars,function(x){x[1,1]})
  for(ii in 1:length(dens)){
    den <- dens[[ii]]
    means <- den$means
    mvar <- den$mvar
    ws <- den$ws
    ## find mise with h=0
    mise0 <- sapply(evars,
                    function(evar){
                      ComputeMISE(0,n,evar,den)})
    ## find mise with h=hx
    draws <- matrix(GenerateData(den,n*N),nrow=N)
    sx <- apply(draws,1,bw.SJ)^2
    out <- vapply(sx,function(s){vapply(evars,function(evar){
        ComputeMISE(s,n,evar,den)},c(0))},rep(0,length(evars)))
    misex <- rowMeans(out)
    ## mise with min of rule of thumb and sx
    out <- matrix(0,nrow=length(evars),ncol=N)
    draws.var <- apply(draws,1,var)
    for(jj in 1:length(evars)){
        sy <- vapply(draws.var,function(x){ComputeHy(n,evars[[jj]][1,1],x)},c(0))^2
        sy <- apply(cbind(sy,sx),1,min)
        out[jj,] <- vapply(sy,function(s){ComputeMISE(s,n,evars[[jj]],den)},c(0))
    }
    misey_smooth <- rowMeans(out)
    ## find mise with h=hy (best procedure)
    misey <- sapply(evars,
                    function(evar){
                      OptimizeMISE(n,evar,den)$objective})
    ## store ratios of errors in results
    results[,ii+1] <- paste("(",format(mise0 / misey,
                                       digits=2,nsmall=2),",",
                            format(misex / misey,
                                   digits=2,nsmall=2),",",
                            format(misey_smooth / misey,
                                   digits=2,nsmall=2),
                            ")",sep="")
  }
  colnames(results) <- c("$\\sigma_{\\epsilon}^2$",
                         lapply(dens,function(x){x$name}))
  return(results)
}

PlotQuantiles <- function(den,
                          n,
                          evar,
                          sx=0,
                          sy=0,
                          q1=.1,
                          q2=.9,
                          n.sim=100,
                          plot0=TRUE,
                          plotx=TRUE,
                          legend.loc='topright'){
  ## plots quantiles for n.sim density estimates of f_Y (= den convolved
  ## with evar) using optimal smoothing and either no smoothing and/or
  ## f_X optimal smoothing (depending on value of plot0 and plotx)
  ##
  ## arguments:
  ##        den : list specifying f_X density (see ConstructX functions for
  ##              details and formatting of den)
  ##          n : sample size
  ##       evar : the variance of the error, a 1 x 1 matrix
  ##         sx : the squared bandwidth for f_X smoothing
  ##         sy : the squared bandwidth for f_Y smoothing
  ##         q1 : the lower quantile to plot (should = 1 - q2)
  ##         q2 : the upper quantile to plot (should = 1 - q1)
  ##      n.sim : the number of density estimates to base quantiles off of
  ##      plot0 : plot / no not plot 0 smoothing quantiles
  ##      plotx : plot / no not plot f_X optimal smoothing quantiles
  ## legend.loc : location of legend   
  ##
  ## value:
  ##    produces a plot
  line.names <- c(expression(f[Y]),expression(h[Y]))
  line.colors <- c("black","orange")
  line.type <- c(1,2)
  ## find x range for plotting
  from <- min(den$means[1,] - 3*sqrt(den$mvar[1,1,] + evar))
  to <- max(den$means[1,] + 3*sqrt(den$mvar[1,1,] + evar))
  ## generate n observations n.sim number of times
  data1 <- matrix(GenerateData(den,n*n.sim),ncol=n.sim)
  ## figure x values chosen by R for constructing density est.
  x <- density(data1[,1],from=from,to=to)$x

  ## estimate fy density using 3 bandwidths
  dx <- apply(data1,2,function(x){density(x,bw=sqrt(sx+evar),
                                          from=from,to=to)$y})
  dx.quants <- apply(dx,1,function(x){quantile(x,probs=c(q1,q2))})
  dy <- apply(data1,2,function(x){density(x,bw=sqrt(sy+evar),
                                          from=from,to=to)$y})
  dy.quants <- apply(dy,1,function(x){quantile(x,probs=c(q1,q2))})
  d0 <- apply(data1,2,function(x){density(x,bw=sqrt(evar),
                                          from=from,to=to)$y})
  d0.quants <- apply(d0,1,function(x){quantile(x,probs=c(q1,q2))})

  ## plot everything
  ymax <- max(dx.quants,dy.quants,d0.quants)
  plot(c(min(x),max(x)),c(0,ymax),col=0,xlab="Y",ylab="Density")
  if (plotx){
    points(x,dx.quants[1,],type='l',col='blue',lwd=3.5,lty=3)
    points(x,dx.quants[2,],type='l',col='blue',lwd=3.5,lty=3)
    line.names[length(line.names)+1] <- expression(h[X])
    line.colors[length(line.colors)+1] <- "blue"
    line.type[length(line.type)+1] <- 3
  }
  points(x,dy.quants[1,],type='l',col='orange',lwd=3,lty=2)
  points(x,dy.quants[2,],type='l',col='orange',lwd=3,lty=2)
  if (plot0){
    points(x,d0.quants[1,],type='l',col='blue',lwd=3.5,lty=3)
    points(x,d0.quants[2,],type='l',col='blue',lwd=3.5,lty=3)
    line.names[length(line.names)+1] <- "h=0"
    line.colors[length(line.colors)+1] <- "blue"
    line.type[length(line.type)+1] <- 3
  }
  ## add error (evar) to density variances, then plot
  evar <- array(evar,c(dim(evar),dim(den$mvar)[3]))
  den$mvar <- den$mvar + evar        
  PlotDen(den,main=FALSE,add=TRUE)
  legend(legend.loc,line.names,col=line.colors,
         lty=line.type,lwd=3,cex=1.5)
}                            


PlotSamples <- function(den,
                        n,
                        evar,
                        s=0,
                        n.sim=10){
  ## plots n.sim density estimates for f_Y
  ## along with true density f_Y
  ##
  ## arguments:
  ##     den : list specifying f_X density (see ConstructX functions for
  ##           details and formatting of den)
  ##       n : sample size
  ##    evar : the variance of the error, a 1 x 1 matrix
  ##       s : the squared bandwidth
  ##   n.sim : the number of density estimates to plot (too many crowds
  ##           the plot, too few gives poor idea of sampling distribution)
  ##
  ## value:
  ##    produces a plot
  ## find x range for plotting
  from <- min(den$means[1,] - 3*sqrt(den$mvar[1,1,] + evar))
  to <- max(den$means[1,] + 3*sqrt(den$mvar[1,1,] + evar))
  ## generate n observations n.sim number of times
  data1 <- matrix(GenerateData(den,n*n.sim),ncol=n.sim)

  ## find the x grid of points
  x <- density(data1[,1],from=from,to=to)$x
  ## find the function values at x grid
  d <- apply(data1,2,function(x){density(x,bw=sqrt(s+evar),
                                         from=from,to=to)$y})

  ## plot lines
  col1 <- "#00000060"
  ymax <- max(d)
  plot(c(from,to),c(0,ymax),col=0,xlab="Y",ylab="Density")
  apply(d,2,function(y){lines(x,y,type='l',col=col1)})
  evar <- array(evar,c(dim(evar),dim(den$mvar)[3]))
  den$mvar <- den$mvar + evar    
  PlotDen(den,add=TRUE)    
}


PlotConvergenceRate <- function(den,evars){
  ## plots ratio of optimal to asymptotic h
  ## as a function of log n,
  ## different lines for diff error distributions
  ##
  ## arguments:
  ##     den : list specifying f_X density (see ConstructX functions for
  ##           details and formatting of den)
  ##   evars : the error variances used for producing this lines
  ##
  ## value:
  ##    produces a plot

  n <- sapply((20:80)/20,function(x){10^x})
  n.evars <- length(evars)
  ## compute h that minimizes mise and asymptotically optimal h
  ## h[,,3] stores ratio of these quantities, which we plot
  hopts <- array(0,c(length(n),n.evars,3))
  for(ii in 1:n.evars){
    hopts[,ii,1] <- ComputeOptimalHfy(den,n,evars[[ii]])
    hopts[,ii,2] <- sqrt(sapply(n,function(x){
      OptimizeMISE(x,evars[[ii]],den)$minimum}))
    hopts[,ii,3] <- hopts[,ii,2] / hopts[,ii,1]
  }
  ## plot h[,,3]
  ymin <- min(hopts[,,3])
  ymax <- max(hopts[,,3])
  par(mar=(c(5,6,4,1)+0.1))
  plot(c(min(n),max(n)),c(ymin,ymax),col=0,log='xy',
       xlab="n",ylab=expression("h"[Y]*"/h"[Y]*"*"),
       main=den$name,
       cex.lab=1.4)
  cols <- rainbow(n.evars)
  abline(h=1,col='grey')
  for(ii in 1:n.evars){
    lines(n,hopts[,ii,3],type='l',
          col=cols[ii],lwd=4,lty=ii)
  }
  legend("bottomright",
         as.character(sapply(evars,function(x){x[1,1]})),
         title="Error Variance",
         col=rainbow(length(evars)),
         lty=1:length(evars),
         lwd=3)
}


## ConstructX functions return a list den
## that specifies a multivariate normal mixture density.
##
## den$means :  p x k matrix where den$means[,ii] is the
##              mean of the ii mixture component
##  den$mvar :  p x p x k array where den$mvar[,,ii] is the
##              covariance of the ii component
##    den$ws :  length k vector where den$ws[ii] is weight of
##              ii mixture component
##  den$name :  name of density, used for plotting and filenames
##
ConstructNormal <- function(){
  den <- list()
  den$means <- matrix(c(0),nrow=1,ncol=1)
  den$mvar <- array(1,c(1,1,1))
  den$ws <- 1
  den$name <- "Normal"
  return(den)
}

ConstructBimodal1 <- function(){
  den <- list()
  den$means <- matrix(c(0,3),nrow=1,ncol=2)
  den$mvar <- array(1,c(1,1,2))
  den$ws <- c(.7,.3)
  den$name <- "Bimodal 1"
  return(den)
}

ConstructBimodal2 <- function(){
  den <- list()
  den$means <- matrix(c(-6,6),nrow=1,ncol=2)
  den$mvar <- array(1,c(1,1,2))
  den$ws <- c(.5,.5)
  den$name <- "Bimodal 2"
  return(den)
}

ConstructTrimodal <- function(){
  den <- list()
  den$means <- matrix(c(-4,0,3),nrow=1,ncol=3)
  den$mvar <- array(c(2,.3,1),c(1,1,3))
  den$ws <- c(.4,.2,.4)
  den$name <- "Trimodal"
  return(den)
}



## See specification of ConstructX functions above
ConstructMultiNormal <- function(){
  den <- list()
  den$means <- matrix(0,nrow=3,ncol=1)
  den$mvar <- array(c(1,0,0,0,1,0,0,0,1),c(3,3,1))
  den$ws <- 1
  den$name <- "Multi. Normal"
  return(den)
}

ConstructMulti2Comp1 <- function(){
  den <- list()
  loc <- .64
  den$means <- matrix(c(0,0,0,1,1,1),nrow=3,ncol=2)
  den$mvar <- array(c(1,loc,0,loc,1,loc,0,loc,1,
                      1,-loc,0,-loc,1,-loc,0,-loc,1),
                    c(3,3,2))
  den$ws <- c(.7,.3)
  den$name <- "Multi. 2-Comp 1"
  return(den)
}  

ConstructMulti2Comp2 <- function(){
  den <- list()
  den$means <- matrix(c(-6,0,0,6,0,0),nrow=3,ncol=2)
  den$mvar <- array(c(1,0,0,0,1,0,0,0,1),c(3,3,2))
  den$ws <- c(.5,.5)
  den$name <- "Multi. 2-Comp 2"
  return(den)
}  

ConstructMulti3Comp <- function(){
  den <- list()
  loc <- .64
  den$means <- matrix(c(0,0,0,1,1,1,0,0,0),nrow=3,ncol=3)
  den$mvar <- array(c(1,loc,0,loc,1,loc,0,loc,1,
                      1,-loc,0,-loc,1,-loc,0,-loc,1,
                      1,-loc,0,-loc,1,-loc,0,-loc,1),
                    c(3,3,3))
  den$ws <- c(.4,.2,.4)
  den$name <- "Multi. 3-Comp"
  return(den)
}  

