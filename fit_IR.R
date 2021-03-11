
rm(list=ls())  # clear memory
library(rstan)
library(tidyverse)
library(bayesplot)
library(deSolve)
library(readxl)

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




data_combined_muc1_IR <-  list(N_T_WT = length(time_data),
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
                               V0 = 30,
                               E_naive0 = 100,
                               E0 = rep(0, 6),
                               B_naive0 = 100,
                               B0 = rep(0,6),
                               A0 = 0,
                               AL0 = 0)







init_condition_IR1 <- list(
  log10_theta = c(log10(8e-6),
                  log10(20), log10(5e-5), 
                  log10(0.1),log10(3.3e+3)),
  theta = c(0.5, 3, 3, 0.5, 3.3e-3),
  sigma = c(1, 1))

init_condition_IR2 <- list(
  log10_theta = c(log10(4e-6),
                  log10(30), log10(1e-5), 
                  log10(0.1),log10(3.3e+3)),
  theta = c(0.1, 3, 10, 0.5, 3.3e-3),
  sigma = c(1, 1))

init_condition_IR3 <- list(
  log10_theta = c(log10(2e-6),
                  log10(10), log10(2e-5), 
                  log10(0.1), log10(3.3e+3)),
  theta = c(0.5, 3, 20, 0.1, 3.3e-3),
  sigma = c(1, 1))

init_condition_IR4 <- list(
  log10_theta = c(log10(1e-6),
                  log10(40), log10(3e-5), 
                  log10(0.1),log10(3.3e+3)),
  theta = c(0.1, 3, 5, 0.1, 3.3e-3),
  sigma = c(1, 1))

options(mc.cores=parallel::detectCores()) # to utilise all cores available in your computer

fit_Model_IR <- stan("IR_prior.stan",
                     data = data_combined_muc1_IR,
                     seed = 20202611,  # set random seed for reproducibility
                     iter = 2000,
                     chains = 4,
                     init = list(init_condition_IR1, init_condition_IR2,init_condition_IR3,init_condition_IR4),
                     warmup = 1000,
                     control = list(adapt_delta = 0.99, max_treedepth = 15))



print(fit_Model_IR , pars = c("theta_WT"))
stan_dens(fit_Model_IR , pars = c("theta_WT"), separate_chains = TRUE,nrow = 6)
# change names ! 


#fit_M2 <- stan("New_model_Stan2_lpf.stan",
#               data = data_combined_muc1_V_W,
#               #pars = c("theta_KO","sigma"),
#               seed = 20201910,  # set random seed for reproducibility
#               iter = 2500,
#               chains = 1,
#               init = list(init5_model2),
#               warmup = 1250,
#               control = list(adapt_delta = 0.99, max_treedepth = 15))



# ASSESSING AND FIXING DIVERGENCES AND TREEDEPTH PROBLEMS

mack_diagnostics <- rstan::get_sampler_params(fit_Model_IR ) %>% 
  set_names(as.factor(1:4)) %>% 
  map_df(as_data_frame, .id = 'chain') %>% 
  group_by(chain) %>% 
  mutate(iteration = 1:length(chain)) %>% 
  mutate(warmup = iteration <= 1000)

mack_diagnostics %>% 
  group_by(warmup, chain) %>% 
  summarise(percent_divergent = mean(divergent__ > 0)) %>% 
  ggplot()+ 
  geom_col(aes(chain, percent_divergent, fill = warmup), position = 'dodge', color = 'black') + 
  scale_y_continuous(labels = scales::percent, name = "% Divergent Runs") + 
  theme_bw() # plot divergent rate for each chain 

mack_diagnostics %>% 
  ggplot(aes(iteration, treedepth__, color = chain)) + 
  geom_line() + 
  geom_hline(aes(yintercept = 15), color = 'red') + 
  theme_bw() # plot treedepth for each chain 


mack_diagnostics %>% 
  ggplot(aes(iteration, stepsize__, color = chain)) + 
  geom_line() + 
  lims(y = c(0,.1)) + 
  theme_bw() # plot stepsize for each chain 


# PARAMETER DIAGNOSTICS
para_summary <- summary(fit_Model_IR )$summary %>% 
  as.data.frame() %>% 
  mutate(variable = rownames(.)) %>% 
  select(variable, everything()) %>% 
  as_data_frame()

# n_eff 
para_summary %>% 
  ggplot(aes(n_eff)) + 
  geom_histogram(binwidth = 50) + 
  geom_vline(aes(xintercept = 4000), color = 'red')

# R_hat

para_summary %>% 
  ggplot(aes(Rhat)) + 
  geom_histogram() + 
  geom_vline(aes(xintercept = 1.01), color = 'red')

# parameter 

para_summary %>% 
  filter(variable %in% c('theta[1]', 'theta[2]', 'theta[3]','theta[4]', 'theta[5]', 'epsilon')) %>% 
  ggplot() + 
  geom_linerange(aes(variable, ymin = `2.5%`,ymax = `97.5%`)) + 
  geom_crossbar(aes(variable, mean, ymin = `25%`, ymax = `75%`), fill= 'grey') + 
  facet_wrap(~variable, scales = 'free') + 
  theme_bw()


# plot sampling parameters 
c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

params_cp <- as.data.frame(rstan::extract(fit_Model_TIV , permuted=FALSE))

names(params_cp) <- gsub("chain:1.", "", names(params_cp), fixed = TRUE)
names(params_cp) <- gsub("[", ".", names(params_cp), fixed = TRUE)
names(params_cp) <- gsub("]", "", names(params_cp), fixed = TRUE)
params_cp$iter <- 1:1000

plot(params_cp$iter, params_cp$theta.10, col=c_dark, pch=16, cex=0.8,
     xlab="Iteration", ylab="epsilon2 (Chain1)", ylim=c(0,1))
plot(params_cp$iter, params_cp$`chain:2.theta.10`, col=c_dark, pch=16, cex=0.8,
     xlab="Iteration", ylab="epsilon2 (Chain2)", ylim=c(0,1))

running_means <- sapply(params_cp$iter, function(n) mean(params_cp$theta.10[1:n]))
plot(params_cp$iter, running_means, col=c_dark, pch=16, cex=0.8, ylim=c(0, 1),
     xlab="Iteration", ylab="epsilon2")



divergent <- get_sampler_params(fit_Model_TIV , inc_warmup=FALSE)[[1]][,'divergent__']
params_cp$divergent <- divergent

div_params_cp <- params_cp[params_cp$divergent == 1,]
nondiv_params_cp <- params_cp[params_cp$divergent == 0,]

plot(div_params_cp$theta.9, div_params_cp$theta.10,
     col="green", pch=16, cex=0.8, xlab="phi", ylab="epsilon2",
     xlim=c(0,50), ylim=c(0,1))
points(nondiv_params_cp$theta.9, nondiv_params_cp$theta.10,
       col=c_dark, pch=16, cex=0.8)


# show mcmc_parcoord
draws <- as.array(fit_Model_IR ,pars = c("theta_WT"))
np <- nuts_params(fit_Model_IR )
str(np)

color_scheme_set("darkgray")
div_style <- parcoord_style_np(div_color = "green", div_size = 0.5, div_alpha = 0.9)
mcmc_parcoord(
  draws,
  transform = function(x) {(x - mean(x)) / sd(x)},
  size = 0.25,
  alpha = 0.3,
  np = np,
  np_style = div_style
)
d <- mcmc_parcoord_data(draws, np = np)
head(d)
tail(d)


# extract posterior samples for selected parameters
posterior_samples_all_IR = rstan::extract(fit_Model_IR , pars = c("theta_WT",'sigma'), inc_warmup = TRUE, permuted = FALSE)
posterior_samples_merged_after_burnin_IR = rstan::extract(fit_Model_IR , pars = c( "theta_WT",'sigma'))

color_scheme_set("brewer-Spectral")
mcmc_trace(posterior_samples_all_IR, n_warmup = 1000,
           facet_args = list(nrow = 4, labeller = label_parsed))

# show all marginal posterior distributions
posterior_sample_table_model_IR = data.frame(epsilon1 = posterior_samples_merged_after_burnin_IR$theta_WT[,1],
                                             beta = posterior_samples_merged_after_burnin_IR$theta_WT[,2],
                                             delta_I = posterior_samples_merged_after_burnin_IR$theta_WT[,3],
                                             p = posterior_samples_merged_after_burnin_IR$theta_WT[,4],
                                             delta_V = posterior_samples_merged_after_burnin_IR$theta_WT[,5],
                                             kappa_M = posterior_samples_merged_after_burnin_IR$theta_WT[,6],
                                             epsilon2 = posterior_samples_merged_after_burnin_IR$theta_WT[,7],
                                             delta_M = posterior_samples_merged_after_burnin_IR$theta_WT[,8],
                                             phi = posterior_samples_merged_after_burnin_IR$theta_WT[,9],
                                             s = posterior_samples_merged_after_burnin_IR$theta_WT[,10])

pairs(posterior_samples_merged_after_burnin_IR$theta_WT[,1:10], labels = c('epsilon1','beta','delta_I','p','delta_V','kappa_M','epsilon2','delta_M', 'phi', 's'), upper.panel = upper.panel, lower.panel = NULL)



# =============== Save the fitting results and generate outputs in CSV file  ===================
write.csv(posterior_samples_merged_after_burnin_IR, file="Posteriors_IR.csv")
save.image(file = "fitting_results_Mu1")
saveRDS(fit_Model_IR, "fit_IR.rds")







