## clearing the environment
rm(list = ls())  
gc()    

require(rstan)
require(loo)
require(tidyverse)
####################################################################################

## model specific details that needs to be change for every run
modelName <- "SH_model"
data_derived1 <- "source_gdt.csv"    # name of the file for precursor pop
data_derived2 <- paste("counts_gdt.csv", sep="")
data_derived3 <- paste("NFd_gdt.csv", sep="")

## setting working dirctory
setwd("/opt/mesh/eigg/sanket/ki67_FM")

## Setting all the directories for opeartions
projectDir <- getwd()
scriptDir <- file.path(projectDir, "scripts")
modelDir <- file.path(projectDir, "models")
dataDir <- file.path(projectDir, "data")
toolsDir <- file.path(scriptDir, "tools")
outDir <- file.path(modelDir, paste(modelName, "_", substr(data_derived1, 1,2), sep=""))
outputDir <- file.path(projectDir, "output/T1_Norm")
saveDir <- file.path(outputDir, paste(modelName, "_", substr(data_derived1, 1,2), sep=""))
figDir <- file.path(projectDir, "deliv", "figures", paste(modelName, "_", substr(data_derived1, 1,2), sep=""))
tabDir <- file.path(projectDir, "deliv", "tables", paste(modelName, "_", substr(data_derived1, 1,2), sep=""))

# loadiong the scr# loadiong the script that contains functions for plotting stan parameters
source(file.path(toolsDir, "stanTools.R"))                # save results in new folder

# compiling multiple stan objects together that ran on different nodes
stanfit1 <- readRDS(file.path(saveDir, "3448_P0_ki67_INC_FM.rds"))
stanfit2 <- readRDS(file.path(saveDir, "3448_P1_ki67_INC_FM.rds"))
stanfit3 <- readRDS(file.path(saveDir, "3448_P2_ki67_INC_FM.rds"))
stanfit4 <- readRDS(file.path(saveDir, "3448_P3_ki67_INC_FM.rds"))

fit <- sflist2stanfit(list(stanfit1, stanfit2, stanfit4, stanfit3))

# finding the parameters used in the model 
# using the last parameter("sigma4") in the array to get the total number of parameters set in the model
num_pars <- which(fit@model_pars %in% "sigma4")      # the variable "sigma4" will change depdending on the data used
parametersToPlot <- fit@model_pars[1:num_pars]

# number of post-burnin samples that are used for plotting 
nPost <- nrow(fit)

################################################################################################
################################################################################################

## loading required datasets for plotting
source_sorted <- read_csv(file.path(dataDir, data_derived1))%>% arrange(time.post.BMT)
counts_sorted <- read_csv(file.path(dataDir, data_derived2))%>% arrange(age.at.S1K)
Nfd_sorted <- read_csv(file.path(dataDir, data_derived3))%>% arrange(time.post.BMT)
ki67_sorted <- read_csv(file.path(dataDir, data_derived4))%>% arrange(time.post.BMT)

# ################################################################################################
# calculating PSIS-L00-CV for the fit
combined_loglik <- extract_log_lik(fit, parameter_name = "log_lik", merge_chains = TRUE)

# optional but recommended
ll_array <- extract_log_lik(fit, merge_chains = FALSE)
r_eff <- relative_eff(exp(ll_array))

# loo-ic values
loo_loglik <- loo(combined_loglik, r_eff = r_eff, save_psis = FALSE, cores = 8)

# Widely applicable AIC
AICw_lok <- waic(combined_loglik)

# AIC from LLmax
#AIC_lok <-  -2 * max(combined_loglik)  + 2 * length(parametersToPlot)

ploocv <- rbind("loo-ic"=loo_loglik$estimates[3], "WAIC" = AICw_lok$estimates[3])  #, "AIC" = AIC_lok)
ploocv[1]

### posterior distributions of parameters
ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
ptable[,c(1,4,8)]

################################################################################################
dir.create(figDir)
dir.create(tabDir)

## writting the parameter values with CIs into a csv file
write.csv(ptable, file = file.path(tabDir, paste(modelName, "_ParameterTable.csv", sep = "")))
write.csv(ploocv, file = file.path(tabDir, paste(modelName, "_Informatiion_criteria.csv", sep = "")))

## open graphics device 
## saving  plots for quality control 
pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
    width = 12, height = 8, onefile = F)

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

log10minorbreaks=as.numeric(1:10 %o% 10^(4:8))

rhats <- rhat(fit, pars = parametersToPlot)
mcmc_rhat(rhats) + yaxis_text() + myTheme

ratios1 <- neff_ratio(fit, pars = parametersToPlot)
mcmc_neff(ratios1) + yaxis_text() + myTheme

posterior <- as.array(fit)
mcmc_acf(posterior, pars = parametersToPlot) + myTheme

mcmcHistory(fit, pars = parametersToPlot, nParPerPage = 4, myTheme = myTheme)

mcmc_dens_overlay(posterior, parametersToPlot)
mcmc_dens(posterior, parametersToPlot) + myTheme

################################################################################################
################################################################################################
## posterior predictive distributions
# time sequnce used for plotting 
ts_pred = seq(from = 0, to = 600, by = 1)

# Total cell counts
Cpred <- as.data.frame(fit, pars = "countspred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)

Y1pred <- as.data.frame(fit, pars = "y1_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)

Y1pred <- Y1pred%>%
  filter(timeseries >= 10)%>% filter(timeseries <= 350)

Cpred <- Cpred%>%
  filter(timeseries >= 10)%>% filter(timeseries <= 350)


ggplot() +
  geom_ribbon(data = Cpred, aes(x = timeseries, ymin = lb, ymax=ub), fill = '#0099cc', alpha = 0.2) +
  geom_ribbon(data = Y1pred, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#f23047", alpha = 0.3)+
  geom_line(data = Y1pred, aes(x = timeseries, y = median), size=1.5) +
  geom_point(data = counts_sorted, aes(x = time.post.BMT, y = total_counts), size=3) +
  labs(title=paste('Cell counts: gdt'),  y=NULL, x= "Time post BMT (days)") + 
  scale_x_continuous(limits = c(10, 400) , trans="log10", breaks=c(10, 30, 100, 300))+
  scale_y_continuous(limits = c(1e4, 2e6), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  guides(color = FALSE) + myTheme


# normalised donr fractions
fdpred <- as.data.frame(fit, pars = "fdpred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)

Y2pred <- as.data.frame(fit, pars = "y2_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)

Y2pred <- Y2pred%>%
  filter(timeseries >= 0)%>% filter(timeseries <= 340)

fdpred <- fdpred%>%
  filter(timeseries >= 0)%>% filter(timeseries <= 340)


ggplot() +
  geom_hline(aes(yintercept = 1), color = "#d11100", linetype = 2, size=1.2)+
  geom_ribbon(data = fdpred, aes(x = timeseries, ymin = lb, ymax=ub), fill = "#0099cc", alpha = 0.2) +
  geom_ribbon(data = Y2pred, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#f23047", alpha = 0.35)+
  geom_line(data = Y2pred, aes(x = timeseries, y = median), size=1.5) +
  geom_point(data = Nfd_sorted, aes(x = time.post.BMT, y = Nfd), size=3) +
  labs(x = "Days post BMT", y = NULL, title = "Normalised Donor fractions: gdt") +
  scale_x_continuous(limits = c(0, 350), breaks = c(0,150,300,450))+
  scale_y_continuous(limits =c(0, 1.2), breaks = c(0, 0.3, 0.6, 0.9, 1.2, 1.5))+ 
  guides(color = FALSE)+ myTheme


zdev.off()

