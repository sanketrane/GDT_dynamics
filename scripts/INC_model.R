## Simple linear ODE model -- Homogeneous model 
## clearing the environment
rm(list = ls())  
gc()    
setwd("/opt/mesh/eigg/sanket/GDT_dynamics")

####################################################################################
## Installing r-stan pachage on the go:
if(!("rstan" %in% rownames(installed.packages())) ){
  install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies=TRUE)
}

## Installing loo:
if(!("loo" %in% rownames(installed.packages())) ){
  install.packages("loo", repos = "https://cloud.r-project.org/", dependencies=TRUE)
}

## Installing tidyverse:
if(!("tidyverse" %in% rownames(installed.packages())) ){
  install.packages("tidyverse", repos = "https://cloud.r-project.org/", dependencies=TRUE)
}
####################################################################################

require(rstan)
require(parallel)
require(loo)
require(tidyverse)

modelName <- "INC_model"
data_derived1 <- "source_gdt.csv"    # name of the file for precursor pop
data_derived2 <- "counts_gdt.csv"
data_derived3 <- "Nfd_gdt.csv"

## Setting all the directories for opeartions
projectDir <- getwd()
scriptDir <- file.path(projectDir, "scripts")
modelDir <- file.path(projectDir, "models")
dataDir <- file.path(projectDir, "data")
toolsDir <- file.path(scriptDir, "tools")
outputDir <- file.path(projectDir, "output")
saveDir <- file.path(outputDir, paste(modelName, "_", substr(data_derived1, 1,2), sep=""))

# loadiong the scr# loadiong the script that contains functions for plotting stan parameters
source(file.path(toolsDir, "stanTools.R"))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# set different seeds for different tasks
i <- as.numeric(Sys.getenv("SLURM_PROCID"))
seed <- 4010 + i
################################################################################################
################################################################################################
source_sorted <- read_csv(file.path(dataDir, data_derived1))%>% arrange(time.post.BMT)
counts_sorted <- read_csv(file.path(dataDir, data_derived2))%>% arrange(time.post.BMT)
Nfd_sorted <- read_csv(file.path(dataDir, data_derived3))%>% arrange(time.post.BMT)

# time points in data
data_time <- counts_sorted$time.post.BMT    
solve_time <- c(0, unique(data_time))   # unique time points to solve ode

#keep track of index of time point in relation to solve_time
time_index <- purrr::map_dbl(data_time, function(x) which(x == solve_time))

# time sequence for prediction
ts_pred = seq(from = 0, to = 600, by = 1)

source_list = data.frame("nu" = 0.001325649, "theta0" = 10.56, "chiEst" = 0.81, "qEst" = 0.066)

## create data set
data <- list(
  numObs = nrow(counts_sorted),
  Nd_0 = 0,
  solve_time = solve_time,
  num_index = length(solve_time),
  time_index = time_index,
  dpBMT = data_time,
  counts = counts_sorted$total_counts,
  Nfd = Nfd_sorted$Nfd,
  numPred = length(ts_pred),
  ts_pred = ts_pred,
  theta0 = exp(source_list$theta0),
  nu = source_list$nu,
  chiEst = source_list$chiEst,
  qEst = source_list$qEst
)


## create initial estimates
init <- function() list(
  psi = exp(rnorm(1, log(0.5), 0.2)),
  lambda = exp(rnorm(1,log(0.02), 1)),
  lambda_inc = exp(rnorm(1,log(0.02), 1)),
  
  y0_Log = rnorm(1, 17 , 0.1),
  alpha = exp(rnorm(1, log(0.5), 0.2)),
  
  sigma1 = exp(rnorm(1,log(1.5), 1)),
  sigma2 = exp(rnorm(1,log(1.5), 1)))

## Specify the variables for which you want history and density plots
parametersToPlot <- c("psi", "lambda", "lambda_inc", "y0_Log", "alpha", "sigma1", "sigma2")

## Additional variables to monitor
otherRVs <- c("y1_mean_pred", "countspred", "y2_mean_pred", "fdpred", "log_lik")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
## run Stan

# parameters for running fits
nChains <- 3
nPost <- 500 ## Number of post-burn-in samples per chain after thinning
nBurn <- 500 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

fit <- stan(file = file.path(modelDir, paste(modelName, ".stan", sep = "")),
            data = data,
            pars = parameters,
            iter = nIter,
            warmup = nBurnin,
            thin = nThin, 
            init = init,
            control = list(adapt_delta = 0.9),
            chains = nChains)

################################################################################################
# save results in output directory
if (!file.exists(saveDir)){
  dir.create(saveDir)
}

# saving output file as 
output_filename=paste(Sys.getenv("SLURM_JOB_ID"), "_P", Sys.getenv("SLURM_PROCID"), "_",  modelName, ".rds",  sep="")

# saving the stan fit object for individual runs as 
saveRDS(fit, file = file.path(saveDir, output_filename))

## calculating an output from individual runs for the validation of a sucessful run
loo_loglik <- extract_log_lik(fit, parameter_name = "log_lik", merge_chains = TRUE)
loo_ic <- loo(loo_loglik)

## this value is added to the job specific run_info file which is gets written in the save directory after each job copmpletion

# information about the cluster run
cluster_run_info <- paste("Time = ", timestamp(),
                          "|| Job ID = ", Sys.getenv("SLURM_JOB_ID"),
                          "|| Modelname = ", modelName,
                          "|| crude_loo_ic = ", loo_ic$estimates[3],
                          "|| Host = ", Sys.getenv("SLURM_SUBMIT_HOST"),
                          "|| Nodes = ",  Sys.info()["nodename"],
                          "|| Process_ID = P", Sys.getenv("SLURM_PROCID"))

# writing the cluster run info as a data file in the out dir 
write.table(cluster_run_info, file = file.path(saveDir, paste("job_", Sys.getenv("SLURM_JOB_ID"), "_run_info.txt")), append = TRUE)
