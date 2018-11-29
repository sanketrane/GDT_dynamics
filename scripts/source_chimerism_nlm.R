## Simple linear model example

rm(list = ls())
gc()

library(rstan)
library(bayesplot)
library(tidyverse)
library(parallel)

modelName <- "source_chimerism_nlm"
data_derived <- "source_gdt.csv"

## Relative paths assuming the working directory is the script directory
## containing this script
projectDir <- getwd()
scriptDir <- file.path(projectDir, "scripts")
figDir <- file.path(projectDir, "deliv", "figures", modelName)
tabDir <- file.path(projectDir, "deliv", "tables", modelName)
dataDir <- file.path(projectDir, "data", data_derived)
modelDir <- file.path(projectDir, "models")
outDir <- file.path(modelDir, modelName)
toolsDir <- file.path(scriptDir, "tools")

#source(file.path(scriptDir, "pkgSetup.R"))
source(file.path(toolsDir, "stanTools.R"))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#set.seed(1027) ## not required but assures repeatable results

################################################################################################

## import the data set
data_imp <- read_csv(dataDir)

ts_pred <- seq(from= 0, to=max(data_imp$time.post.BMT), by=1)

## create data set
data <- with(data_imp,
             list(
               numObs = nrow(data_imp),
               Time = data_imp$time.post.BMT,
               numPred = length(ts_pred),
               ts_pred = ts_pred,
               chiT2 = data_imp$fd
               ))

## create initial estimates
init <- function() list(
    chiEst = runif(1, 0, 1),
    qEstlog = rnorm(1, -2, 1),
    sigma = exp(rnorm(1,log(1.5), 1)))

## Specify the variables for which you want history and density plots
parametersToPlot <- c("chiEst","qEst","sigma")

## Additional variables to monitor
otherRVs <- c("ymean_pred","chipred")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
# run Stan

nChains <- 10
nPost <- 1000 ## Number of post-burn-in samples per chain after thinning
nBurn <- 1000 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

dir.create(outDir)

fit <- stan(file = file.path(modelDir, paste(modelName, ".stan", sep = "")),
            data = data,
            pars = parameters,
            iter = nIter,
            warmup = nBurnin,
            thin = nThin, 
            init = init,
            control = list(adapt_delta = 0.9),
            chains = nChains)

save(fit, file = file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))
##load(file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))

dir.create(figDir)
dir.create(tabDir)

## open graphics device
pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
    width = 8, height = 5, onefile = F)

pairs(fit, pars = parametersToPlot)

options(bayesplot.base_size = 15,
        bayesplot.base_family = "sans")
color_scheme_set(scheme = "viridis")
myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 22), axis.title =  element_text(size = 20, face = "bold"),
                 plot.title = element_text(size=20,  hjust = 0.5), legend.text = element_text(size=20), legend.title = element_text(size = 20))

# setting ggplot theme for rest fo the plots
theme_set(theme_bw())

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (e^+2 -> e^2)
  l <- gsub("e\\+","e",l)  
  # turn the 'e' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text=l)
}

log10minorbreaks=as.numeric(1:10 %o% 10^(3:6))

rhats <- rhat(fit, pars = parametersToPlot)
mcmc_rhat(rhats) + yaxis_text() + myTheme

ratios1 <- neff_ratio(fit, pars = parametersToPlot)
mcmc_neff(ratios1) + yaxis_text() + myTheme

posterior <- as.array(fit)
mcmc_acf(posterior, pars = parametersToPlot) + myTheme

mcmcHistory(fit, pars = parametersToPlot, nParPerPage = 4, myTheme = myTheme)

mcmc_dens_overlay(posterior, parametersToPlot)
mcmc_dens(posterior, parametersToPlot) + myTheme

ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
write.csv(ptable, file = file.path(tabDir, paste(modelName, "ParameterTable.csv", sep = "")))

################################################################################################
## posterior predictive distributions

Ypred <- as.data.frame(fit, pars = "ymean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = data$ts_pred)

fdpred <- as.data.frame(fit, pars = "chipred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = data$ts_pred)

ggplot() +
  geom_hline(aes(yintercept = 1), color = "#d11100", linetype = 2, size=1.2)+
  geom_ribbon(data = fdpred, aes(x = timeseries, ymin = lb, ymax=ub), fill = "#F05E82", alpha = 0.2) +
  geom_ribbon(data = Ypred, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#EC2F5E", alpha = 0.35)+
  geom_line(data = Ypred, aes(x = timeseries, y = median), size=1.5) +
  geom_point(data = data_imp, aes(x = time.post.BMT, y = fd), size=3) +
  labs(x = "Days post t0", y = NULL, title = "Normalised Donor fractions: TH.gdt") +
  scale_x_continuous(limits = c(0, 340), breaks = c(0,150,300,450))+
  scale_y_continuous(limits =c(-0.35, 1.25), breaks = c(0, 0.3, 0.6, 0.9, 1.2, 1.5))+ 
  guides(color = FALSE)+ theme(axis.text = element_text(size = 22),
                               axis.title =  element_text(size = 20, face = "bold"),
                               plot.title = element_text(size=20,  hjust = 0.5),
                               legend.text = element_text(size=20),
                               legend.title = element_text(size = 20))


dev.off()
