### 
### 
### DENSITY ESTIMATES FOR NO2 EXPOSURE IN WATERTOWN 
### (data from article by Aurore Delaigle)
### 
### by James Long 
### date:  October 13, 2013
### 

##
## reference notes:
## DELAIGLE = "Nonparametric density estimation from data
## with a mixture of Berkson and classical errors" by
## Aurore Delaigle in Canadian Journal of Statistics
##

rm(list=ls(all=TRUE))
set.seed(20130712)

## functions.R is loaded so that
## we can check that optimal h is being used
## correctly by comparing ComputeHy (see below)
## to the more general function ComputeOptimalHfy
## in functions.R, basically we use functions.R as
## a unit test
source("functions.R")

## read in data
conc <- read.table("data/HEALTH.txt")

## column names come from Readme_plus_PORTAGE.txt
names(conc) <- c("id","wheeze","wheeze.pers","kitchen","bedroom")

## using parameter estimates in DELAIGLE section 3.2
## compute NO2 exposure = a_0 + a_1*log(kitchen) + a_2*log(bedroom)
conc$kit.bed <- 1.22 + .3*log(conc$kitchen) + .33*log(conc$bedroom)


## Note: NO2 exposure in children is modeled as
## conc$kit.bed + epsilon where epsilon is standard normal.
## here we estimate the density for children's exposure
## trying several variances for epsilon, including
## the estimate of .06 used in DELAIGLE


## plot of density without any error
no.error <- density(conc$kit.bed)

pdf("figs/watertown_no_error.pdf")
plot(no.error,lwd=2,main="Density Estimate Ignoring Error")
dev.off()

n <- nrow(conc)
ConstructNormal <- function(){
  den <- list()
  den$means <- matrix(c(0),nrow=1,ncol=1)
  den$mvar <- array(1,c(1,1,1))
  den$ws <- 1
  den$name <- "Normal"
  return(den)
}
den <- ConstructNormal()
evar <- matrix(.06,nrow=1,ncol=1)

## check that ComputeHy is working by
## comparing to function ComputeOptimalHfy
## which computes optimal Hy when fx is
## gaussian mixture
ComputeOptimalHfy(den,n,evar)
ComputeHy(n,evar[1,1],1)


##  Creates a plot with no smoothing, hx smoothing, and hy smoothing lines
##     obs : vector with the observations (i.i.d. draws from f_X)
##     v.e : variance of the error epsilon (assumed normal, mean 0)
MakePlot <- function(obs,v.e){
    ## compute the three densities (no smoothing, hx, and hy)
    no.smooth <- density(obs,bw=sqrt(v.e))
    hx <- bw.SJ(obs)
    hx.smooth <- density(obs,bw=sqrt(v.e + hx^2))
    hy <- min(ComputeHy(length(obs),v.e,var(obs)),hx)
    hy.smooth <- density(obs,bw=sqrt(v.e + hy^2))
    ## limits of plot 
    y.max <- max(no.smooth$y,hx.smooth$y,hy.smooth$y)
    x.min <- min(no.smooth$x,hx.smooth$x,hy.smooth$x)
    x.max <- max(no.smooth$x,hx.smooth$x,hy.smooth$x)
    ## make plot
    par(mar=c(5,5,.5,.5))
    plot(c(x.min,x.max),c(0,y.max),
         ylab="Density",xlab="Y",col=0,cex.lab=2)
    lines(no.smooth,col='red',lwd=3,lty=1)
    lines(hx.smooth,col='blue',lwd=3,lty=2)
    lines(hy.smooth,col='black',lwd=3,lty=3)
    hy.r <- round(hy,2)
    hx.r <- round(hx,2)
    legend.names <- c(expression("h=0"),
                      substitute(paste(tilde(h)[X],"=",hx.r),list(hx.r=hx.r)),
                      substitute(paste(tilde(h)[Y],"=",hy.r),list(hy.r=hy.r)))
    legend("topleft",legend.names,
           lty=c(1,2,3),col=c("red","blue","black"),
           lwd=3,cex=2)
}

## make 3 plots with error variance (.06)
## used in DELAIGLE as well as 10 times and
## 1/10 this value i.e. (.6 and .006)
v.e <- .06 * 10
pdf("figs/error_large.pdf")
MakePlot(conc$kit.bed,v.e)
dev.off()

v.e <- .06
pdf("figs/error_true.pdf")
MakePlot(conc$kit.bed,v.e)
dev.off()

v.e <- .06 / 10
pdf("figs/error_small.pdf")
MakePlot(conc$kit.bed,v.e)
dev.off()
