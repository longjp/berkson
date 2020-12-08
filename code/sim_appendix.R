### 
### 
### CONSTRUCT TABLES AND PLOTS FOR
### ''KERNEL DENSITY ESTIMATION WITH BERKSON ERROR''
### by Long, El Karoui, and Rice
###
### MAKES TABLE IN APPENDIX WHERE BANDWIDTHS ARE
### ESTIMATED
###
### by: James Long 
### date: October 18, 2015
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


##
## ONE DIMENSIONAL ANALYSIS
##

## load all the densities into list dens
dens <- list(ConstructNormal(),
             ConstructBimodal1(),
             ConstructBimodal2(),
             ConstructTrimodal())

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

## make results table with n obs
n <- 50
tab <- MakeRelativeErrorTableEst(dens,n,evars)
cap <- "This table repeats the simulations of Table \\ref{tab:50} but with estimated bandwidths. Each entry is $\\left(\\frac{MISE(0)}{MISE(h_Y)},\\frac{MISE(\\widetilde{h}_X)}{MISE(h_Y)},\\frac{MISE(\\widetilde{h}_Y)}{MISE(h_Y)}\\right)$ for $n=50$. These ratios are always greater than $1$ because $h_Y$ is the minimizer of the $\\MISE$."
x.tab <- xtable(tab,
                align="rr|cccc",
                caption=cap,
                label=paste("tab:est",n,sep=""))
print(x.tab,
      file=paste("figs/est",n,".tex",sep=""),
      type='latex',
      include.rownames=FALSE,
      sanitize.text.function = function(x){x})


## make results table with n obs
n <- 100
tab <- MakeRelativeErrorTableEst(dens,n,evars)
cap <- "The simulations here are identical to Table \\ref{tab:est50} but for n=100."
x.tab <- xtable(tab,
                align="rr|cccc",
                caption=cap,
                label=paste("tab:est",n,sep=""))
print(x.tab,
      file=paste("figs/est",n,".tex",sep=""),
      type='latex',
      include.rownames=FALSE,
      sanitize.text.function = function(x){x})
