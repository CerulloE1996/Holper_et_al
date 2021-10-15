#########################################################################################################
## Holper et al. analysis
## Code authors: Enzo Cerullo (using Stan code from Nyaga et al. - DOI: 10.1177/0962280216669182)
## Date: 02/10/2021
#########################################################################################################



# Set local working directory (change to your local working directoy) ---------------------------------------------------
setwd("/home/enzo/Documents/CRSU/Holper_et_al_collab") # linux 
setwd("C:/Users/Enzo/Documents/CRSU/Holper_et_al_collab") # windows 


# Source files needed -----------------------------------------------------
source('Holper_et_al_functions.R')

# Load R packages ---------------------------------------------------------
require(dplyr)
require(plyr)
require(ggplot2)
require(rstan)
require(patchwork)
require(data.table)
require(groupdata2)
require(bayesSurv)
require(scales)
require(RColorBrewer)
require(loo)
require(pbmcapply)
require(bayesSurv)
require(scales)



# Data for NMA and non-NMA multiple test models  --------------------------------------------------------

length(unique(X2$test_num))

mean_prev <- mean(c(0.11, 0.07, 0.10, 0.11, 0.08))
mean_prev

max_n_test_thr <- length(unique(X2$test_threshold_num))
max_n_test_thr

test_seq <- c(rep(1,3), rep(2,4), rep(3,5), rep(4,5), rep(5,5))
test_seq

data <-  list(N = nrow(X2),
              n_tests = length(unique(X2$test_num)) , 
              n_thresholds = c(3,4,5,5,5), # number of categories for each test
              n_test_thresholds = length(unique(X2$test_threshold_num)) ,
              fu = as.numeric(X2$follow_up_centered_5yrs),
              cts_cov_points = seq(from = min(X2$follow_up), to = max(X2$follow_up), length = 100),
              TP = X2$TP,
              FN = X2$FN, 
              TN = X2$TN,
              FP = X2$FP, 
              Test = X2$test_num,
              Threshold = X2$threshold_num,
              Test_Threshold = X2$test_threshold_num,
              Study = X2$study , 
              prev = rep(mean_prev, 5), 
              youden_weight = array(0.7), 
              holdout = rep(0,   nrow(X2) ))



#  ------------------------------------------------------------ ------------------------------------------------------------
# Model 0D - Non-NMA / direct comparison model w/ Compound symmetry (CS) var-cov structure  ------------------------------------------------------------
#  ------------------------------------------------------------ ------------------------------------------------------------
# CS i.e. implies the same SD and correlation between all 5 tests 

# compile model 
file_non_NMA <- stan_model(file = "./models/non_NMA_follow_up.stan") 

# prior set 1
data$prior_intercept_thr_1_diseased_mean <- 8
data$prior_intercept_thr_1_non_diseased_mean <- -8
data$prior_intercept_thr_1_sd <- 5
data$prior_intercept_thr_2_diseased_mean <- 1.5
data$prior_intercept_thr_2_non_diseased_mean <- -1.5
data$prior_intercept_thr_2_sd <- 2
data$prior_intercept_other_thr_sd <- 2
data$prior_coeff_sd <- 0.5
data$prior_sigma_sd <- 2

# prior set 2
data$prior_intercept_thr_1_diseased_mean <- 0
data$prior_intercept_thr_1_non_diseased_mean <- 0
data$prior_intercept_thr_1_sd <- 10
data$prior_intercept_thr_2_diseased_mean <- 0
data$prior_intercept_thr_2_non_diseased_mean <- 0
data$prior_intercept_thr_2_sd <- 3
data$prior_intercept_other_thr_sd <- 2
data$prior_coeff_sd <- 1
data$prior_sigma_sd <- 3


data$test_thr_param <- 0
data$fu_summary <- 0

data$diff_SDs_corr <- 1 # 1 for seperate Sds and corr between test, 0 for same

# compound symmetry (CS) var-cov matrix for the tets random effect (NMA, diff SDs/corr's between tests) 
model_non_NMA_CS <-     run_stan_model(model_file = file_non_NMA,
                                    diff_SDs_corr = 0,
                                    data = data, 
                                    warmup = 500,
                                    sampling = 500, 
                                    adapt_delta = 0.80, 
                                    max_treedepth = 10)

get_num_divergent(model_non_NMA_CS$model_output)
get_num_max_treedepth(model_non_NMA_CS$model_output)

#saveRDS(model_non_NMA_CS, file = "./outputs/model_non_NMA_CS_output_fu.rds") # save output

print(model_non_NMA_CS$model_output, pars = c("Se", "Sp"), probs = c(0.025, 0.5, 0.975))


print(model_non_NMA_CS$model_output, pars = c("Youden_index"), probs = c(0.025, 0.5, 0.975))
# youden - 3,3,4 (almost 3 ~ 0.02 off), 4 (almost 3 ~ 0.02 off), 4 (almost 3 ~0.04 off)

print(model_non_NMA_CS$model_output, pars = c("mu","sigma"), probs = c(0.025, 0.5, 0.975))
print(model_non_NMA_CS$model_output, pars = c("fu_coeff"), probs = c(0.025, 0.5, 0.975))

# Se and Sp 
model_non_NMA_CS$stan_list$stan_se
model_non_NMA_CS$stan_list$stan_sp

# Se and Sp only at the optimal thr (according to youden) 
model_non_NMA_CS$out_list$Se_youden
model_non_NMA_CS$out_list$Sp_youden

# pred intervals
model_non_NMA_CS$stan_list$stan_se_pred
model_non_NMA_CS$stan_list$stan_sp_pred

# differences 
model_non_NMA_CS$out_list$diff_Se_youden
model_non_NMA_CS$out_list$diff_Sp_youden

# trace plots
stan_trace(model_non_NMA_CS$model_output, pars = c("Se", "Sp"))
stan_trace(model_non_NMA_CS$model_output, pars = c("sigma", "Sigma"))

# posterior density plots
stan_dens(model_non_NMA_CS$model_output, pars = c("Se", "Sp"))
stan_dens(model_non_NMA_CS$model_output, pars = c("sigma", "Sigma"))
stan_dens(model_non_NMA_CS$model_output, pars = c("fu_coeff"))
stan_dens(model_non_NMA_CS$model_output, pars = c("mu"))



#  ------------------------------------------------------------ ------------------------------------------------------------
# Model 0E - Non-NMA / direct comparison model w/ Unstructued (UN) var-cov structure  ------------------------------------------------------------
#  ------------------------------------------------------------ ------------------------------------------------------------
# UN -> different SDs and corr for each of the 5 tests 

model_non_NMA_UN <-     run_stan_model(model_file = file_non_NMA,
                                    diff_SDs_corr = 1,
                                    data = data, 
                                    warmup = 1000,
                                    sampling = 1000, 
                                    adapt_delta = 0.80, 
                                    max_treedepth = 10)

get_num_divergent(model_non_NMA_UN$model_output)
get_num_max_treedepth(model_non_NMA_UN$model_output)

saveRDS(model_non_NMA_UN, file = "./outputs/model_non_NMA_UN_output_fu.rds") # save output

print(model_non_NMA_UN$model_output, pars = c("Se", "Sp"), probs = c(0.025, 0.5, 0.975))

print(model_non_NMA_UN$model_output, pars = c("mu", "sigma"), probs = c(0.025, 0.5, 0.975))

print(model_non_NMA_UN$model_output, pars = c("Youden_index"), probs = c(0.025, 0.5, 0.975))
# youden -  3,3 (almost 4, ~ 0.03 off),4,4,4 (almost 3, ~ 0.03 off)

# Se and Sp 
model_non_NMA_UN$stan_list$stan_se
model_non_NMA_UN$stan_list$stan_sp

# Se and Sp only at the optimal thr (according to youden) 
model_non_NMA_UN$out_list$Se_youden
model_non_NMA_UN$out_list$Sp_youden

# differences 
model_non_NMA_UN$out_list$diff_Se_youden
model_non_NMA_UN$out_list$diff_Sp_youden


# trace plots
stan_trace(model_non_NMA_UN$model_output, pars = c("Se", "Sp"))
stan_trace(model_non_NMA_UN$model_output, pars = c("sigma", "Sigma"))

# posterior density plots
stan_dens(model_non_NMA_UN$model_output, pars = c("Se", "Sp"))
stan_dens(model_non_NMA_UN$model_output, pars = c("sigma", "Sigma"))



# ---------------------------------------------------------------------------------------------------------------------
# Model comparison - Model 0D vs 0E ------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
model_0D <- readRDS(file = "./outputs/model_non_NMA_CS_output_fu.rds") # CS
model_0E <- readRDS(file = "./outputs/model_non_NMA_UN_output_fu.rds") # UN

# Using loo 
loo_function(model_0D$model_output) # loo = 1389.4
loo_function(model_0E$model_output) # loo = 1391.0

# suggests approx same fit, however high proportion of pareto-k > 0.7 (i.e. loo is not valid in this case)
# therefore would need to do full k-fold cross-validation





