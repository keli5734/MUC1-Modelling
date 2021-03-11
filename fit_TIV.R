
rm(list=ls())  # clear memory
library(rstan)
library(tidyverse)
library(bayesplot)
library(deSolve)
library(readxl)
library(wesanderson)
library(latex2exp)

WT_viral_load_data <- read_excel("PR8-viral-load.xls", sheet = 1) %>% rename(Days = Days, ViralLoad = `Viral load`)
KO_viral_load_data <- read_excel("PR8-viral-load.xls", sheet = 2) %>% rename(Days = Days, ViralLoad = `Viral load`)

WT_macrophage_data <- read_excel("Macrophage.xls", sheet = 1) %>% rename(Days = Days, Macrophage = 'Macro level')
KO_macrophage_data <- read_excel("Macrophage.xls", sheet = 2) %>% rename(Days = Days, Macrophage = 'Macro level')


time_data <- c(1,2,3,5,7)
V_data_WT <- WT_viral_load_data$ViralLoad
V_data_KO <- KO_viral_load_data$ViralLoad


time_data_macrophage <- c(1,3,5,7)
Macrophage_data_WT <- WT_macrophage_data$Macrophage
Macrophage_data_KO <- KO_macrophage_data$Macrophage




data_combined_muc1_TIV <-  list(N_T_WT = length(time_data),
                                       N_V_WT = length(V_data_WT),
                                       time_data_WT = time_data,
                                       log_viral_load_data_WT = log(V_data_WT), # above V_WT
                                       N_T_KO = length(time_data),
                                       N_V_KO = length(V_data_KO),
                                       time_data_KO = time_data,
                                       log_viral_load_data_KO = log(V_data_KO), # above V_KO
                                       N_T_Macrophage_WT = length(time_data_macrophage),
                                       N_Macrophage_WT = length(Macrophage_data_WT),
                                       time_data_Macrophage_WT = time_data_macrophage,
                                       Macrophage_data_WT = log(Macrophage_data_WT), # above Macrophage_WT
                                       N_T_Macrophage_KO = length(time_data_macrophage),
                                       N_Macrophage_KO = length(Macrophage_data_KO),
                                       time_data_Macrophage_KO = time_data_macrophage,
                                       Macrophage_data_KO = log(Macrophage_data_KO), # above Macrophage_KO
                                       t0 = 0,
                                       T0 = 1e+7,
                                       I0 = 0,
                                       V0 = 30)





## =========== Viral parameters only ============== ##


init_condition_TIV1 <- list(
  log10_theta = c(log10(8e-6),
                  log10(20), log10(5e-5), 
                  log10(0.1),log10(3.3e+3)),
  theta = c(0.5, 3, 3, 0.5, 3.3e-3),
  sigma = c(1, 1))

init_condition_TIV2 <- list(
  log10_theta = c(log10(4e-6),
                  log10(30), log10(1e-5), 
                  log10(0.1),log10(3.3e+3)),
  theta = c(0.1, 3, 10, 0.5, 3.3e-3),
  sigma = c(1, 1))

init_condition_TIV3 <- list(
  log10_theta = c(log10(2e-6),
                  log10(10), log10(2e-5), 
                  log10(0.1), log10(3.3e+3)),
  theta = c(0.5, 3, 20, 0.1, 3.3e-3),
  sigma = c(1, 1))

init_condition_TIV4 <- list(
  log10_theta = c(log10(1e-6),
                  log10(40), log10(3e-5), 
                  log10(0.1),log10(3.3e+3)),
  theta = c(0.1, 3, 5, 0.1, 3.3e-3),
  sigma = c(1, 1))

options(mc.cores=parallel::detectCores()) # to utilise all cores available in your computer

fit_Model_TIV <- stan("Muc1_TIV.stan",
                                    data = data_combined_muc1_TIV,
                                    seed = 20202611,  # set random seed for reproducibility
                                    iter = 2000,
                                    chains = 4,
                                    init = list(init_condition_TIV1,init_condition_TIV2,init_condition_TIV3,init_condition_TIV4),
                                    warmup = 1000,
                                    control = list(adapt_delta = 0.99, max_treedepth = 15))



print(fit_Model_TIV , pars = c("theta_WT"))
stan_dens(fit_Model_TIV , pars = c("theta_WT"), separate_chains = TRUE,nrow = 6)
# change names ! 




# =============== Save the fitting results and generate outputs in CSV file  ===================
write.csv(posterior_samples_merged_after_burnin_TIV, file="Posteriors_TIV.csv")
save.image(file = "fitting_results_Mu1")












