## Simple linear model example
rm(list = ls())
gc()
setwd("~/Desktop/Git_repos/GDT_dynamics")

library(rstan)
library(bayesplot)
library(tidyverse)
library(parallel)

modelName <- "source_counts_linear"
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

set.seed(1027) ## not required but assures repeatable results

################################################################################################

## import the data set
data_imp <- read_csv(dataDir)%>%
  arrange(age.at.S1K)

ts_pred <- seq(from= min(data_imp$age.at.S1K), to=max(data_imp$age.at.S1K), by=1)

## create data set
data <- list(
  numObs = nrow(data_imp),
  Time = data_imp$age.at.S1K,
  counts = data_imp$total_counts,
  numPred = length(ts_pred),
  ts_pred = ts_pred
  )

## create initial estimates
init <- function() list(
    mu = rnorm(1, 0.001, 0.02),
    y0Log = rnorm(1, 11, 0.1),
    sigma = exp(rnorm(1, log(1.5), 0.5)))

## Specify the variables for which you want history and density plots
parametersToPlot <- c("mu","y0Log","sigma")

## Additional variables to monitor
otherRVs <- c("countspred", "ymean_pred")

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

################################################################################################
## posterior distributions of parameters

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
# Total cell counts
Cpred <- as.data.frame(fit, pars = "countspred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)

Y1pred <- as.data.frame(fit, pars = "ymean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)


ggplot() +
  geom_ribbon(data = Cpred, aes(x = timeseries, ymin = lb, ymax= ub), fill = '#F05E82', alpha = 0.2) +
  geom_ribbon(data = Y1pred, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#EC2F5E", alpha = 0.3)+
  geom_line(data = Y1pred, aes(x = timeseries, y = median), size=1.5) +
  geom_point(data = data_imp, aes(x = age.at.S1K, y = total_counts), size=3) +
  labs(title=paste('Cell counts: TH.gdt' ),  y=NULL, x="Host age (days)") + 
  scale_x_continuous(limits = c(60, 460), trans = "log10",  breaks = c(10,30,100,300))+
  scale_y_continuous(limits = c(1e3, 2e5), trans="log10", breaks=c(1e4, 1e5, 1e3), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  guides(color = FALSE) + myTheme

dev.off()
