## get rid of old plots/tables
rm(list=ls())
unlink("figs",recursive=TRUE)
unlink("../ms/figs",recursive=TRUE)
dir.create("figs")
dir.create("../ms/figs")

## run simulations and real data
rm(list=ls())
source('sim.R')
rm(list=ls())
source('sim_appendix.R')
rm(list=ls())
source('watertown.R')

## copy output to ../ms/figs/
rm(list=ls())
fs <- list.files("figs",full.names=TRUE)
file.copy(fs,to="../ms/figs/")
