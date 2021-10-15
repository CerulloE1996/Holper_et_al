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


# Set options -------------------------------------------------------
rstan_options(auto_write = TRUE)
options(mc.cores = 12)

options(scipen = 999)
options(max.print = 1000000000)

filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
rename <- dplyr::rename

# Load data --------------------------------------------------------------
X <- tibble(read.csv('table.csv'))
head(X)
length(unique(X$s))

# Prepare data for use with model -----------------------------------------
X2 <- X %>%
  mutate(N = TP+FP+FN+TN,
         test_num = as.numeric(as.factor(Tool)),
         se = round(TP/(TP+FN),2), 
         sp = round(TN/(TN+FP),2),
         follow_up = as.numeric(Followup),
         follow_up_centered_5yrs = follow_up - 5) %>%
  rename(study = s,
         test = Tool,
         threshold_num = Category_Num) %>%
  select(study, Category, test, test_num, threshold_num,  follow_up, follow_up_centered_5yrs,TP, FP, FN, TN, N, se, sp) %>% 
  dplyr::arrange(study,test_num, threshold_num,.by_group = TRUE) %>%
  mutate(test_threshold_num = case_when( 
    # STABLE-2007 (3 cats) - Low, Moderate, High
    test_num == 1 & threshold_num == 1 ~ 1, 
    test_num == 1 & threshold_num == 2 ~ 2, 
    test_num == 1 & threshold_num == 3 ~ 3, 
    # Static-99 (4 cats) - Low, Moderate-Low, Moderate-High, High
    test_num == 2 & threshold_num == 1 ~ 4, 
    test_num == 2 & threshold_num == 2 ~ 5, 
    test_num == 2 & threshold_num == 3 ~ 6, 
    test_num == 2 & threshold_num == 4 ~ 7,
    # Static-99/STABLE-2007 (5 cats)  - Low, Moderate-Low, Moderate-High, High, Very high
    test_num == 3 & threshold_num == 1 ~ 8, 
    test_num == 3 & threshold_num == 2 ~ 9, 
    test_num == 3 & threshold_num == 3 ~ 10, 
    test_num == 3 & threshold_num == 4 ~ 11,
    test_num == 3 & threshold_num == 5 ~ 12,
    # Static-99R (5 cats) - very low, below average, average, above average, well above average
    test_num == 4 & threshold_num == 1 ~ 13, 
    test_num == 4 & threshold_num == 2 ~ 14, 
    test_num == 4 & threshold_num == 3 ~ 15, 
    test_num == 4 & threshold_num == 4 ~ 16,
    test_num == 4 & threshold_num == 5 ~ 17,
    # Static-99R/STABLE-2007 (5 cats) - Low, Moderate-Low, Moderate-High, High, Very high
    test_num == 5 & threshold_num == 1 ~ 18, 
    test_num == 5 & threshold_num == 2 ~ 19, 
    test_num == 5 & threshold_num == 3 ~ 20, 
    test_num == 5 & threshold_num == 4 ~ 21,
    test_num == 5 & threshold_num == 5 ~ 22))


filter(X2, follow_up == 100) # 2 studies have missing FU

X2 <- filter(X2, follow_up != 100) # remove study w/ missing follow-up

unique(X2$study)
length(unique(X2$study))

# re-label study variable 
X2 <- X2 %>%  
  mutate(study = ifelse(study < 21 & study > 11, study-1, study))  %>%
  mutate(study = ifelse(study > 20, study-2, study))

unique(X2$study)
length(unique(X2$study))

head(X2)
#View(X2)



filter(X2, follow_up == 100) # 2 studies (7 obs.) have missing info for follow up time 

filter(X2, study == 3 & follow_up != 8.2)
filter(X2, study == 3 & follow_up == 8.2)

filter(X2, study == 4 & follow_up != 9)
filter(X2, study == 4 & follow_up == 9)




### Note that the variable test_num corresponds to:
# 1= STABLE-2007 (3 thr) - Low, Moderate, High
# 2= Static-99 (4 thr) - Low, Moderate-Low, Moderate-High, High
# 3= Static-99/STABLE-2007 (5 thr)  - Low, Moderate-Low, Moderate-High, High, Very high
# 4= Static-99R (5 thr) - very low, below average, average, above average, well above average
# 5= Static-99R/STABLE-2007 (5 thr) - Low, Moderate-Low, Moderate-High, High, Very high

unique(X2$follow_up)
length(unique(X2$study))


# Data for NMA  --------------------------------------------------------
mean_prev <- mean(c(0.11, 0.07, 0.10, 0.11, 0.08))
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
              Test = (X2$test_num),
              Test_Threshold = X2$test_threshold_num,
              Threshold = (X2$threshold_num),
              Study = X2$study , 
              prev = rep(mean_prev, 5), 
              youden_weight = array(0.6), 
              holdout = rep(0,   nrow(X2) ))




# see how many studies report each test 
for (t in 1:data$n_tests) { 
    dataset <- filter(X2, test_num == t)
    n_studies <- length(unique(dataset$study))
    print(paste( n_studies ,"studies report test", t))
}
# we can see that few studies report test 1,3 and 5 



# Exploring Priors  -----------------------------------------------------------------

prior_mean_function <- function(mean, sd) { 
  
  print(paste("intervals for mean=",mean,"sd=",sd))
  
  print(paste("50th %-ile",median(sort(plogis(rnorm(n=10000, mean = mean, sd = sd))))))
  print(paste("60th %-ile",sort(plogis(rnorm(n=10000, mean = mean, sd = sd)))[6000]))
  print(paste("75th %-ile",sort(plogis(rnorm(n=10000, mean = mean, sd = sd)))[7500]))
  print(paste("97.5th %-ile",sort(plogis(rnorm(n=10000, mean = mean, sd = sd)))[9750]))
  print(paste("99.8th %-ile",sort(plogis(rnorm(n=10000, mean = mean, sd = sd)))[9980]))
  print(paste("25th %-ile",sort(plogis(rnorm(n=10000, mean = mean, sd = sd)))[2500]))
  print(paste("2.5th %-ile",sort(plogis(rnorm(n=10000, mean = mean, sd = sd)))[250]))
  print(paste("0.02th %-ile",sort(plogis(rnorm(n=10000, mean = mean, sd = sd)))[20]))
  
  plot(density(plogis(rnorm(n=1000, mean = mean, sd = sd))))
  
}

par(mfrow = c(1,1))

# for threshold 1 Se
prior_mean_function(8,5)
prior_mean_function(0,10)

# for threshold 2 Se
prior_mean_function(-1.5, 2) # (0.084, 0.996), median = 0.817 for Se and median = 0.18, (0.004, 0.915) for Sp
prior_mean_function(0, 3)



prior_mean_function(0,1.5)
prior_mean_function(0,2)
prior_mean_function(0,3)
prior_mean_function(0,4)

plogis(-3+2*2) - plogis(-3+2*1) # 0.46
plogis(-3+2*3) - plogis(-3+2*2) # 0.22 - > i.e. the change in Se w/ 1 year follow-up decreases w/ increasing time 
plogis(-3+2*3) - plogis(-3+2*2)

plogis(-3+2*log(2)) - plogis(-3+2*log(1)) # 0.11
plogis(-3+2*log(3)) - plogis(-3+2*log(2)) # 0.14
plogis(-3+2*log(4)) - plogis(-3+2*log(3)) # 0.13




# ------------------------------------------------------------ ------------------------------------------------------------ ------------------------------------------------------------
#### Model 1A and 1B - NMA model w/ CS or UN var-cov structure   ---------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------ ------------------------------------------------------------ ------------------------------------------------------------

file_NMA <- stan_model(file = "./models/NMA_alt_param_2_follow_up_CS.stan")
file_NMA <- stan_model(file = "./models/NMA_alt_param_2_follow_up_UN.stan")

data$cts_cov_points <- seq(from = min(data$fu), to =  max(data$fu), length = 100)

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
data$prior_tau_sd <- 2


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
data$prior_tau_sd <- 3



data$fu_summary <- 0 # note that since using fu centered @5yrs, this means summary estimates are adjusted for 5 year follow-up



model_1A <- run_stan_model(    model_file = file_NMA,
                               diff_SDs_corr = 0,
                               data = data,
                               warmup = 500,
                               sampling = 500,
                               adapt_delta = 0.90,
                               max_treedepth = 10,
                               refresh = 100)


get_num_divergent(model_1A$model_output)
get_num_max_treedepth(model_1A$model_output)

saveRDS(model_1A, file = "./outputs/model_1A_prior_set_1_output.rds") # save output



print(model_1A$model_output, pars = c("Se", "Sp"), probs = c(0.025, 0.5, 0.975))
print(model_1A$model_output, pars = c("Youden_index"), probs = c(0.025, 0.5, 0.975))

# Se and Sp only at the optimal thr (according to youden) 
model_1A$out_list$Se_youden
model_1A$out_list$Sp_youden

print(model_1A$model_output, pars = c("mu"), probs = c(0.025, 0.5, 0.975))
print(model_1A$model_output, pars = c("tau", "sigma", "Sigma", "Omega"), probs = c(0.025, 0.5, 0.975))


# Se and Sp 
model_1A$stan_list$stan_se
model_1A$stan_list$stan_sp

# differences and ratios (according to thresholds using youden)
model_1A$out_list$diff_Se_youden
model_1A$out_list$diff_Sp_youden

model_1A$out_list$ratio_Se_youden
model_1A$out_list$ratio_Sp_youden


# Se/Sp at max according to DOR 
model_1A$out_list$Se_dor
model_1A$out_list$Sp_dor

# Se/Sp at max according to weighted youden  
model_1A$out_list$Se_w_youden
model_1A$out_list$Sp_w_youden


# trace plots
stan_trace(model_1A$model, pars = c("Se", "Sp"))
stan_trace(model_1A$model, pars = c("sigma", "Sigma"))

# posterior density plots
stan_dens(model_1A$model, pars = c("Se", "Sp"))
stan_dens(model_1A$model, pars = c("sigma", "Sigma"))



loo_function(model_1A$model_output) # loo = 1523.5




# ---------------------------------------------------------------------------------------------------------------------
# Model comparison - Model 1A vs 1B ------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
model_1A <- readRDS(file = "./outputs/model_1A_prior_set_1_output.rds") # CS
model_1B <- readRDS(file = "./outputs/model_1B_prior_set_1_output.rds") # UN

print(model_1A$model_output, pars = c("Se", "Sp"), probs = c(0.025, 0.5, 0.975))
print(model_1B$model_output, pars = c("Se", "Sp"), probs = c(0.025, 0.5, 0.975))


# Using loo 
loo_function(model_1A$model_output) # loo = 1379.4
loo_function(model_1B$model_output) # loo = 1376.5

# suggests UN gives better fit, however high proportion of pareto-k > 0.7 (i.e. loo is not valid in this case)

# therefore need to do full k-fold cross-validation




# Model comparison - K-fold cross-validation --------------------------------------------------------------------------------------------------------------------
set.seed(123)

# put data in tibble format
d <-  tibble( TP = X2$TP,
              FN = X2$FN, 
              TN = X2$TN,
              FP = X2$FP, 
              pos = X2$TP + X2$FN,
              neg = X2$TN + X2$FP,
              Test = X2$test_num,
              Threshold = X2$threshold_num,
              Test_Threshold = X2$test_threshold_num,
              fu = as.factor(X2$follow_up_centered_5yrs),
              Study = as.factor(X2$study))

d


# Assuming K = 10, the procedure is the following:
# 1) We extract 1/10 of the data and save it in the list G with k = 1;
# 2) We extract 1/9 of the remaining data and save it with k = 2;
# 3) We extract 1/8 of the remaining data and save it with k = 3;
# 4) ...;
# 10) We extract all the data the data and save it with k = 10

K <- 10
dK <- fold(d, k = K, 
           id_col = "Study",
           method = "n_rand")
dK

# Create list containing the K datasets
ldata <- plyr::llply(1:K, function(i) {
  list( N = nrow(X2),
        prior_intercept_thr_1_diseased_mean = 8,
        prior_intercept_thr_1_non_diseased_mean = 8,
        prior_intercept_thr_1_sd = 5,
        prior_intercept_thr_2_diseased_mean = 1.5,
        prior_intercept_thr_2_non_diseased_mean = -1.5,
        prior_intercept_thr_2_sd = 2,
        prior_intercept_other_thr_sd = 2,
        prior_coeff_sd = 0.5,
        prior_sigma_sd = 2,
        prior_tau_sd = 2,
        fu_summary = 0,
        diff_SDs_corr = 0,
        n_tests = data$n_tests,
        n_thresholds = c(3,4,5,5,5), # number of categories for each test
        n_test_thresholds = length(unique(X2$test_threshold_num)) ,
        fu = X2$follow_up_centered_5yrs,
        cts_cov_points = seq(from = min(X2$follow_up), to = max(X2$follow_up), length = 100),
        TP = X2$TP, 
        FN = X2$FN, 
        TN = X2$TN, 
        FP = X2$FP, 
        pos = X2$TP + X2$FN,
        neg = X2$TN + X2$FP,
        Test = X2$test_num,
        Threshold = X2$threshold_num,
        Test_Threshold = X2$test_threshold_num,
        Study = X2$study,
        prev = rep(mean_prev, 5), 
        youden_weight = array(0.7), 
        holdout = ifelse(dK$.folds == i, 1, 0))
})

#ldata

#  run k-fold  (with 1 chain only - may take a while ) -------------------------------------------------------------------------------------------------


# First model 
K <- length(ldata)
model <- stan_model(file="./models/NMA_alt_param_2_follow_up_CS.stan")
chains  <- 1
cores <- 1
iter <- 2000
warmup <- 1000

# First parallelize all chains:
sflist <-  pbmclapply(1:(K*chains), 
                      mc.cores = cores, 
                      function(i){
                        # Fold number:
                       # k <- round((i+1) / chains)
                        k <- i
                        s <- sampling(model,
                                      data = ldata[[k]], 
                                      chains = 1, 
                                      chain_id = i,  
                                      seed = 123, 
                                      iter = iter,
                                      warmup = warmup, 
                                      control=list(adapt_delta=0.80, 
                                                   max_treedepth = 10))
                        return(s)
                      })


stanfit <- list()
for(k in 1:K){
  stanfit[[k]] <- sflist2stanfit(sflist[k])
}  


holdout <- lapply(ldata, '[[', "holdout")

# We extract all the held_out log_lik of all the folds
log_lik_ab <- extract_log_lik_K(stanfit, holdout, "log_lik")
length <- length(log_lik_ab[1,][log_lik_ab[1,] != "NaN" & !is.na(log_lik_ab[1,]) ])
loglik2 <- array(data = log_lik_ab[log_lik_ab != "NaN" & !is.na(log_lik_ab) ], dim = c((iter-warmup), length ))

kfold_m1 <- kfold(loglik2)
kfold_m1$elpd_kfold  # -949.4073

#saveRDS(kfold_m1, file = "kfold_m1.rds")





# Second model 

model <- stan_model(file="./models/NMA_alt_param_2_follow_up_UN.stan")


# First parallelize all chains:
sflist <-  pbmclapply(1:(K*chains), 
                      mc.cores = cores, 
                      function(i){
                        # Fold number:
                        # k <- round((i+1) / chains)
                        k <- i
                        s <- sampling(model,
                                      data = ldata[[k]], 
                                      chains = 1, 
                                      chain_id = i,  
                                      seed = 123, 
                                      iter = iter,
                                      warmup = warmup, 
                                      control=list(adapt_delta=0.80, 
                                                   max_treedepth = 10))
                        return(s)
                      })


stanfit <- list()
for(k in 1:K){
  stanfit[[k]] <- sflist2stanfit(sflist[k])
}  

holdout <- lapply(ldata, '[[', "holdout")

# We extract all the held_out log_lik of all the folds
log_lik_ab <- extract_log_lik_K(stanfit, holdout, "log_lik")
length <- length(log_lik_ab[1,][log_lik_ab[1,] != "NaN" & !is.na(log_lik_ab[1,]) ])
loglik2 <- array(data = log_lik_ab[log_lik_ab != "NaN" & !is.na(log_lik_ab) ], dim = c((iter-warmup), length ))

kfold_m2 <- kfold(loglik2)

kfold_m2$elpd_kfold  # -972.8671

#saveRDS(kfold_m2, file = "kfold_m2.rds")


# Compare models ----------------------------------------------------------

kfold_m1_CS <- readRDS(file = "kfold_m1.rds")
kfold_m2_UN <- readRDS(file = "kfold_m2.rds")

kfold_m1_CS$elpd_kfold*(-2) # 1898.815
kfold_m2_UN$elpd_kfold*(-2) # 1945.734

## computer differences in elpd_kfold and se of differences 
kfold_m1_CS$elpd_kfold - kfold_m2_UN$elpd_kfold
sqrt(nrow(kfold_m1_CS$pointwise) * var(kfold_m1_CS$pointwise - kfold_m2_UN$pointwise))


#  run k-fold using 4 chains - faster but prone to errors on some computers / versions of R  b -------------------------------------------------------------------------------------------------


fits_m1 <- stan_kfold("./models/NMA_alt_param_2_follow_up_CS.stan",
                   list_of_datas = ldata,
                   chains = 4,
                   cores = 4,
                   seed = 123,
                   iter = 1000,
                   warmup = 500,
                   control=list(adapt_delta=0.80,
                                max_treedepth = 10))








