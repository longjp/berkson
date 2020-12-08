### 
### 
### CONSTRUCT TABLES AND PLOTS FOR
### ''KERNEL DENSITY ESTIMATION WITH BERKSON ERROR''
### by Long, El Karoui, and Rice
###
### by: James Long 
### date: July 12, 2013
### 


rm(list=ls(all=TRUE))
## July 12, 2013 - first date of work on this file
set.seed(20130712)

## computing / optimizing MISE and bandwidth 
## along with generating plots performed
## by functions in functions.R
source("functions.R")
## for outputting tables in latex
library('xtable')
## for 3-d plotting of draws from the multivariate densities
library('rgl')

## load all the densities into list dens
dens <- list(ConstructNormal(),
             ConstructBimodal1(),
             ConstructBimodal2(),
             ConstructTrimodal())


## plot all densities in one figure
pdf("figs/densities.pdf")
par(mfrow=c(2,2))
for(den in dens){
    PlotDen(den,lwd=2)
}
dev.off()



### tables of MISE ratios, specifically
### (MISE(0)/MISE(h_Y), MISE(h_X)/MISE(h_Y))
### for densities in dens with
### several error variances (evars)
### and 50 and 100 observations

evars <- c(2,1,.5,.25,.125)
## each error variance must be a matrix for
## compatibility with multi-d case, so turn
## evars into a list of 1x1 matrices
evars <- lapply(evars,function(x){matrix(x,nrow=1,ncol=1)})

## 50 observations, tab contains the results
n <- 50
tab <- MakeRelativeErrorTable(dens,n,evars)
cap <- paste("n=",n,".",sep="")
cap <- "Each entry is $\\left(\\frac{MISE(0)}{MISE(h_Y)},\\frac{MISE(h_X)}{MISE(h_Y)}\\right)$ for $n=50$. These ratios are always greater than $1$ because $h_Y$ is the minimizer of the $\\MISE$. As expected, $\\MISE(0)$ performs well when $\\sigma_{\\epsilon}^2$ (the error variance) is large but poorly when $\\sigma_{\\epsilon}^2$ is small. $\\MISE(h_X)$ performs well when $\\sigma_{\\epsilon}^2$ is small but poorly when $\\sigma_{\\epsilon}^2$ is large."
x.tab <- xtable(tab,
                align="rr|cccc",
                caption=cap,
                label=paste("tab:",n,sep=""))
print(x.tab,
      file=paste("figs/",n,".tex",sep=""),
      type='latex',
      include.rownames=FALSE,
      sanitize.text.function = function(x){x})


## 100 observations, tab contains the results
n <- 100
tab <- MakeRelativeErrorTable(dens,n,evars)
cap <- "The entries here are the same as Table \\ref{tab:50} but for $n=100$. This larger $n$ generally improves performance for $\\MISE(0)$ and worsens the performance of $\\MISE(h_X)$ (relative to $\\MISE(h_Y)$). This is predicted by our asymptotic theory, since as $n \\rightarrow \\infty$, $\\frac{MISE(0)}{MISE(h_Y)} \\rightarrow 1$ while $\\frac{MISE(h_X)}{MISE(h_Y)} \\rightarrow \\infty$. However at $n=100$, using $h_X$ still generally outperforms no smoothing."
x.tab <- xtable(tab,
                align="rr|cccc",
                caption=cap,
                label=paste("tab:",n,sep=""))
print(x.tab,
      file=paste("figs/",n,".tex",sep=""),
      type='latex',
      include.rownames=FALSE,
      sanitize.text.function = function(x){x})

## case study: normal with error var 2
## here hx is much worse than hy
## plot quantiles and n.sim density estimates
n <- 50
evar <- matrix(2,nrow=1,ncol=1)
den <- ConstructNormal()
sx <- OptimizeMISE(n,matrix(0,nrow=1,ncol=1),den)$minimum
sy <- OptimizeMISE(n,evar,den)$minimum
print(sqrt(sx))
print(sqrt(sy))
n.sim <- 10

pdf("figs/normal_conf.pdf",width=6,height=6)
PlotQuantiles(den,n,evar,sx,sy,plot0=FALSE)
dev.off()

pdf("figs/normal_samples_hx.pdf",height=6,width=6)
PlotSamples(den,n,evar,s=sx,n.sim=n.sim)
dev.off()

pdf("figs/normal_samples_hy.pdf",height=6,width=6)
PlotSamples(den,n,evar,s=sy,n.sim=n.sim)
dev.off()


#### case study: bimodal1 with error var .125
#### here h0 is much worse than hy
#### plot quantiles and n.sim density estimates 
n <- 50
evar <- matrix(.125,nrow=1,ncol=1)
den <- ConstructBimodal1()
sx <- OptimizeMISE(n,matrix(0,nrow=1,ncol=1),den)$minimum
sy <- OptimizeMISE(n,evar,den)$minimum
print(sqrt(sx))
print(sqrt(sy))
n.sim <- 10

pdf("figs/bimodal1_conf.pdf",width=6,height=6)
PlotQuantiles(den,n,evar,sx=sx,sy=sy,plotx=FALSE)
dev.off()

pdf("figs/bimodal1_samples_h0.pdf",height=6,width=6)
PlotSamples(den,n,evar,s=0,n.sim=n.sim)
dev.off()

pdf("figs/bimodal1_samples_hy.pdf",height=6,width=6)
PlotSamples(den,n,evar,s=sy,n.sim=n.sim)
dev.off()



## plot (exact optimal h) / (asymptotically optimal h)
## as a function of n for different error distributions

## choose error variances
evars <- c(2,1,.5,.25,.125)
evars <- lapply(evars,function(x){matrix(x,nrow=1,ncol=1)})
## for each density create plot
for(den in dens){
    fname <- gsub(" ","",paste("figs/",den$name,"_h_conv.pdf",sep=""))
    pdf(fname,height=6,width=6)
    PlotConvergenceRate(den,evars)
    dev.off()
}



##
## MULTI (3) DIMENSIONAL ANALYSIS
##

dens <- list(ConstructMultiNormal(),
             ConstructMulti2Comp1(),
             ConstructMulti2Comp2(),
             ConstructMulti3Comp())

## make 3-d plots of draws from the four densities
## can be viewed in java compliant web browswer by
## opening index.html in appropriate folder in figs/
for(den in dens){
  draws <- GenerateData(den,500)
  plot3d(draws[,1],draws[,2],draws[,3],
         xlab="x1",ylab="x2",zlab="x3",col="red", size=3)
  writeWebGL(paste("figs/",gsub(" ","",den$name),sep=""),
             width=1000,height=1000)
  rgl.close()
}

## construct error table with n=100
n <- 100
evars <- c(2,1,.5,.25,.125)
evars <- lapply(evars,function(x){matrix(c(x,0,0,0,x,0,0,0,x),nrow=3,ncol=3)})
tab <- MakeRelativeErrorTable(dens,n,evars)
cap <- "Three dimensional finite sample results for $n=100$. Generally, $h_X$ and no smoothing perform worse relative to $h_Y$ here than for $n=100$ in one dimension (see Table \\ref{tab:100})."
x.tab <- xtable(tab,align="rr|cccc",
                caption=cap,
                label=paste("tab:multi",n,sep=""))
print(x.tab,file=paste("figs/multi",n,".tex",sep=""),
      type='latex',include.rownames=FALSE,
      sanitize.text.function = function(x){x})


## construct error table with n=500
n <- 500
evars <- c(2,1,.5,.25,.125)
evars <- lapply(evars,function(x){matrix(c(x,0,0,0,x,0,0,0,x),nrow=3,ncol=3)})
tab <- MakeRelativeErrorTable(dens,n,evars)
cap <- "Three dimensional finite sample results for $n=500$. $\\MISE(h_X) / MISE(h_Y)$ is larger and $\\MISE(0)/MISE(h_Y)$ is smaller here relative to Table \\ref{tab:multi100} where the sample size was $100$."
x.tab <- xtable(tab,align="rr|cccc",
                caption=cap,
                label=paste("tab:multi",n,sep=""))
print(x.tab,file=paste("figs/multi",n,".tex",sep=""),
      type='latex',include.rownames=FALSE,
      sanitize.text.function = function(x){x})
