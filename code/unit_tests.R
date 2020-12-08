### 
### 
### SET OF UNIT-LIKE TESTS FOR functions.R 
###  
###
### by James Long 
### date: NOVEMBER 5, 2013 
### 

source('functions.R')


## 
## TEST 1: Compare ComputeIntF2 to numerical integration
##

## numerical integration function
Func <- function(den,deriv=1,points=5000){
    ## lower limit for integration
    lower <- min(den$means[1,] - 4*sqrt(den$mvar[1,1,]))
    ## upper limit for integration
    upper <- max(den$means[1,] + 4*sqrt(den$mvar[1,1,]))
    x <- ((0:points)/points)*(upper - lower) + lower
    vals <- sapply(x,function(d){
      dnorm.mixt(d,
                 mus=den$means[1,],
                 sigmas=sqrt(den$mvar[1,1,]),
                 props=den$ws)})
    deriv1 <- diff(vals) / (x[2]-x[1])
    if (deriv==1){
        return(sum(deriv1^2)*(x[2]-x[1]))
    }
    if (deriv==2){
        deriv2 <- diff(deriv1) / (x[2]-x[1])
        return(sum(deriv2^2)*(x[2]-x[1]))
    }
}

dens <- list(ConstructNormal(),
             ConstructBimodal1(),
             ConstructBimodal2(),
             ConstructTrimodal())


## do Func and ComputeIntF2 produce same output?
for(ii in dens){
    ## first derivative
    print("first deriv test")
    print(Func(ii,deriv=1))
    print(ComputeIntF2(ii,deriv=1))
    ## second derivative
    print("second deriv test")
    print(Func(ii,deriv=2))
    print(ComputeIntF2(ii,deriv=2))
}








#####
##### TEST 2: A few tests for ComputeOmega
##### 

## should throw error b/c variance is 0
a <- 0
S <- 0
evar <- matrix(0,nrow=1,ncol=1)
means <- matrix(0,nrow=1,ncol=1)
mvar <- array(0,c(1,1,1))
ComputeOmega(a,S,evar,means,mvar)

## output should be equal
a <- 0
S <- 0
evar <- matrix(3,nrow=1,ncol=1)
means <- matrix(0,nrow=1,ncol=1)
mvar <- array(0,c(1,1,1))
ComputeOmega(a,S,evar,means,mvar)
dnorm(0,0,sqrt(6))

## 2 dimensional, independent, output equal
a <- 0
S <- 0
evar <- matrix(c(3,0,0,3),nrow=2,ncol=2)
means <- matrix(0,nrow=2,ncol=1)
mvar <- array(0,c(2,2,1))
ComputeOmega(a,S,evar,means,mvar)
dnorm(0,0,sqrt(6))^2





#####
##### TEST 3: Visual tests of OptimizeMISE, first
##### run function and then plot to make sure we are
##### at the minimum
#####
dens <- list(ConstructNormal(),
             ConstructBimodal1(),
             ConstructBimodal2(),
             ConstructTrimodal())


evars <- c(2,1,.5,.25,.125)
evars <- lapply(evars,function(x){matrix(x,nrow=1,ncol=1)})

## plot MISE and minimum found by OptimizeMISE
n <- 50
par(mfcol=c(5,4))
for(den in dens){
  for(evar in evars){
    opt <- OptimizeMISE(n,evar,den)
    Ss <- ((0:100) / 100 + .5) * opt$minimum
    responses <- vapply(Ss,function(S){ComputeMISE(S,n,evar,den)},0)
    plot(Ss,responses,main=paste("evar:",evar[1,1],"  den:", den$name,sep=""))
    abline(v=opt$minimum)
    abline(h=opt$objective)
  }
}


n <- 100
par(mfcol=c(5,4))
evar0 <- matrix(0,nrow=1,ncol=1)
for(den in dens){
  for(evar in evars){
    opt <- OptimizeMISE(n,evar,den)
    Ss <- ((0:100) / 100 + .5) * opt$minimum
    responsesy <- vapply(Ss,function(S){ComputeMISE(S,n,evar,den)},0)
    plot(Ss,responses,main=paste("evar:",evar[1,1],"  den:", den$name,sep=""))
    abline(v=opt$minimum)
    abline(h=opt$objective)
  }
}

dens <- list(ConstructMultiNormal(),
             ConstructMulti2Comp1(),
             ConstructMulti2Comp2(),
             ConstructMulti3Comp())

evars <- c(2,1,.5,.25,.125)
evars <- lapply(evars,function(x){matrix(c(x,0,0,0,x,0,0,0,x),nrow=3,ncol=3)})


n <- 100
par(mfcol=c(5,4))
for(den in dens){
  for(evar in evars){
    opt <- OptimizeMISE(n,evar,den)
    Ss <- ((0:100) / 100 + .5) * opt$minimum
    responsesy <- vapply(Ss,function(S){ComputeMISE(S,n,evar,den)},0)
    plot(Ss,responses,main=paste("evar:",evar[1,1],"  den:", den$name,sep=""))
    abline(v=opt$minimum)
    abline(h=opt$objective)
  }
}


n <- 500
par(mfcol=c(5,4))
evar0 <- matrix(0,nrow=1,ncol=1)
for(den in dens){
  for(evar in evars){
    opt <- OptimizeMISE(n,evar,den)
    Ss <- ((0:100) / 100 + .5) * opt$minimum
    responsesy <- vapply(Ss,function(S){ComputeMISE(S,n,evar,den)},0)
    plot(Ss,responses,main=paste("evar:",evar[1,1],"  den:", den$name,sep=""))
    abline(v=opt$minimum)
    abline(h=opt$objective)
  }
}
