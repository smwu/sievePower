library("rslurm")

# Calculation of power for 60 scenarios: 
# (AMP-B, AMP-C, pooled) x (2:1, 1:1) x (IC50, IC80) x (0, 0.3, 0.5, 0.7, 0.9)

rm(list=ls(all=TRUE))

dataDir <- "../data"
# Ted, you obviously need to change 'codeDir' to wherever you have your local repo
codeDir <- "."

library(survival)
library(np)

source(file.path(codeDir, "../markHR.R"))
source(file.path(codeDir, "../covEst.R"))
source(file.path(codeDir, "markDensSim_tholzman_20181029.R"))

dataB <- read.csv(file.path(dataDir, "catnap_vrc01_neut_b.csv"))
dataC <- read.csv(file.path(dataDir, "catnap_vrc01_neut_c.csv"))
dataBandC <- read.csv(file.path(dataDir, "catnap_vrc01_neut_all.csv"))

#nMC <- 1
nMC <- 1000
pars <- data.frame(
    seed=seq(1,nMC)
    )
pars$trial <- rep(I(list(c("704", "703", "704and703"))),times=nMC)
pars$randRatio <- rep(I(list(c("2:1", "1:1"))),times=nMC)
pars$IC <- rep(I(list(c("IC50", "IC80"))),times=nMC)
pars$VEcoord <- rep(I(list(c(0, 0.3, 0.5, 0.7, 0.9))),times=nMC)
pars$taumax <- rep(80,times=nMC)
pars$alpha1sided <- rep(0.05,times=nMC)

results <- slurm_apply(
    getInferenceOneMC,
    pars,
    add_objects=c("dataB","dataC","dataBandC"),
    jobname="markPower",
    nodes=1000,
    cpus_per_node=1,
    submit=FALSE
    )

print(results)
save(results, file="results.RData")

