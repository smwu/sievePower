.tmplib <- lapply(c('base','methods','datasets','utils','grDevices','graphics','stats','rslurm','survival','np'), 
                  library, character.only = TRUE, quietly = TRUE)

source("../../markHR.R")
source("../../covEst.R")
source("../markDensSim_tholzman_20181029.R")

load('add_objects.RData')
.rslurm_func <- function(seed, trial, randRatio, IC, VEcoord, taumax, alpha1sided){
  scenario <- c(outer(c(outer(c(outer(trial, randRatio, paste1)), IC, paste1)), VEcoord, paste1))
  out <- vector("list", length=length(scenario))
  names(out) <- scenario
  
  for (i in 1:length(trial)){
    if (trial[i]=="704"){ 
      data <- dataB
      # 'Np' is number of subjects in placebo group
      Np <- 900
      # 'np' is the expected number of cases in placebo group
      np <- 34
    } else if (trial[i]=="703"){ 
      data <- dataC
      Np <- 634
      np <- 34
    } else {
      data <- dataBandC
      Np <- 900 + 634
      np <- 68
    }
    
    # rate parameters for failure time T and censoring time C
    lambdaT <- (log(1 - (1 + 0.1*Np/np)*(np/Np)))/(-taumax*(1 + 0.1*Np/np))
    lambdaC <- 0.1*Np/np*lambdaT
    
    for (k in 1:length(IC)){
      mark <- data[,paste0(tolower(IC[k]),".geometric.mean.imputed.log10")]
      
      # nonparametric kernel density estimate of the mark variable
      densbw <- npudensbw(~ mark, ckertype="epanechnikov")
      dens <- npudens(densbw)
      
      for (l in 1:length(VEcoord)){
        coeff <- getAlphaBetaGamma(VEcoord[l], dens, "mark")
        
        for (j in 1:length(randRatio)){
          scenarioLabel <- paste0(trial[i],"_",randRatio[j],"_",IC[k],"_",VEcoord[l])
          out[[scenarioLabel]] <- simulOne(Np=Np, np=np, lambdaT=lambdaT, lambdaC=lambdaC, alpha=coeff$alpha, beta=coeff$beta, gamma=coeff$gamma, 
                                           taumax=taumax, dens=dens, varName="mark", randomRatio=ifelse(randRatio[j]=="2:1",2,1), seed=seed)
        }
      }
    }
  }
  
  return(out)
}

.rslurm_params <- readRDS('params.RDS')
.rslurm_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
.rslurm_istart <- .rslurm_id * 1 + 1
.rslurm_iend <- min((.rslurm_id + 1) * 1, nrow(.rslurm_params))
.rslurm_result <- do.call(parallel::mcMap, c(.rslurm_func,
    .rslurm_params[.rslurm_istart:.rslurm_iend, , drop = FALSE],
    mc.cores = 1))
               
saveRDS(.rslurm_result, file = paste0('results_', .rslurm_id, '.RDS'))
