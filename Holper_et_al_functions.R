

# Function for model comparison (using loo) - for 4 chains  --------------------------------------------------------


loo_function <- function(model_object) { 
  
  log_lik <- extract_log_lik(model_object, parameter_name = "log_lik")
  str(log_lik) 
  
  draws <-  dim(log_lik)[1]/4
  length <- length(log_lik[1,][log_lik[1,] != "NaN" & !is.na(log_lik[1,]) ])
  loglik2 <- array(data = log_lik[log_lik != "NaN" & !is.na(log_lik) ], dim = c(draws*4, length ))
  
  r_eff <- relative_eff(exp(loglik2), cores = 4, chain_id = c(rep(1, draws), rep(2, draws), rep(3, draws), rep(4,draws)))
  
  mod_loo <- loo(loglik2, cores = 4, r_eff = r_eff)
  
  return(mod_loo)
}



# make function to run models with run w/ rstan  ----------------------------------------------------------------------------------


run_stan_model <- function(model_file, 
                           diff_SDs_corr,
                           data, 
                           warmup, 
                           sampling, 
                           adapt_delta, 
                           max_treedepth,
                           refresh) { 
  
  
  
  data$diff_SDs_corr <- diff_SDs_corr
  
  model_rstan <<-  rstan::sampling(
    object = model_file,
    data =  data,
    chains = 4,
    cores = 4,
    iter = (warmup + sampling),
    warmup = warmup,
    control=list(adapt_delta=adapt_delta,
                 max_treedepth = max_treedepth),
    seed= 123,
    refresh = refresh
  )
  
  
  stan_se <- round(rstan::summary(model_rstan,  pars= c("Se"),  probs = c(0.025,  0.5, 0.975))$summary, 2)
  stan_sp <- round(summary(model_rstan,  pars= c("Sp"),  probs = c(0.025,  0.5, 0.975))$summary, 2)
  stan_se_pred <- round(summary(model_rstan,   pars= c( "Se_pred"), probs = c(0.025,  0.5, 0.975))$summary, 2)
  stan_sp_pred <- round(summary(model_rstan,   pars= c( "Sp_pred"), probs = c(0.025,  0.5, 0.975))$summary, 2)
  stan_ppv <- round(summary(model_rstan,  pars= c("PPV"),  probs = c(0.025,  0.5, 0.975))$summary, 2)
  stan_npv <- round(summary(model_rstan,  pars= c("NPV"),  probs = c(0.025,  0.5, 0.975))$summary, 2)
  stan_LRpos <- round(summary(model_rstan,  pars= c("LRpos"), probs = c(0.025,  0.5, 0.975))$summary, 2)
  stan_LRneg <- round(summary(model_rstan,  pars= c("LRneg"), probs = c(0.025,  0.5, 0.975))$summary, 2)
  
  stan_var_cov <- round(summary(model_rstan, 
                               pars= c("sigma", 
                                       "Sigma",
                                       # NMA models only
                                       "tau",
                                       "sigmasq"
                               ),
                               probs = c(0.025,  0.5, 0.975))$summary, 2)
  
  stan_corr <- round(summary(model_rstan, 
                            pars= c("Omega",
                                    # NMA models only 
                                    "rho", "rho12"
                            ),
                            probs = c(0.025,  0.5, 0.975))$summary, 2)
  
  
  stan_youden <- round(summary(model_rstan,
                              pars=      c("Youden_index"),
                              probs = c(0.025,  0.5, 0.975))$summary, 2)
  
  stan_youden_weighted <- round(summary(model_rstan,
                                       pars=      c("Youden_index_weighted"),
                                       probs = c(0.025,  0.5, 0.975))$summary, 2)
  
  stan_DOR <- round(summary(model_rstan,
                           pars=      c("DOR"),
                           probs = c(0.025,  0.5, 0.975))$summary, 2)
  
  
  
  # extract full posterior distributions from model output to calculate other statistics of interest
  params <- rstan::extract(model_rstan)
  
  n_samps <- 4*sampling
  
  max_category <- max(data$n_thresholds)
  
  Se <- array(data = NA, dim = c(n_samps, data$n_tests, max_category))
  Sp <- array(data = NA, dim = c(n_samps, data$n_tests, max_category))
  Fp <- array(data = NA, dim = c(n_samps, data$n_tests, max_category))
  
  for (t in 1:data$n_tests){
    for (t_c in 1:data$n_thresholds[t]) { 
      Se[,t,t_c] <- params$Se[,t,t_c]
      Sp[,t,t_c] <- params$Sp[,t,t_c]
      Fp[,t,t_c] <- params$Fp[,t,t_c]
    }
  }
  
  
  
  # calculate AUC for each test
  AUC_height <- array(data = NA, dim = c(n_samps, data$n_tests,max_category+1))
  AUC_width <- array(data = NA, dim = c(n_samps, data$n_tests,max_category+1))
  AUC_height_x_width <- array(data = NA, dim = c(n_samps, data$n_tests,max_category+1))
  AUC_samps <- array(data = NA, dim = c(n_samps, data$n_tests))
  AUC <- tibble(Test = seq(1,data$n_tests),
                mean = rep(NA, data$n_tests), 
                sd =  rep(NA, data$n_tests),
                lower  = rep(NA, data$n_tests),
                median  = rep(NA, data$n_tests),
                upper  = rep(NA, data$n_tests))
  
  for (t in 1:data$n_tests) {
    
    
    for (n in 1:n_samps) {
      AUC_height[n,t,1:(data$n_thresholds[t]+1)] = ( c(1,Se[n,t,1:data$n_thresholds[t]],0)[-1] 
                                                       + c(1,Se[n,t,1:data$n_thresholds[t]],0)[-length(c(1,Se[n,t,1:data$n_thresholds[t]],0))]   ) / 2
      AUC_width[n,t,1:(data$n_thresholds[t]+1)] = -diff(c(1,Fp[n,t,1:data$n_thresholds[t]],0) )
    }
    
    
    for (i in 1:(data$n_thresholds[t]+1)) {
      AUC_height_x_width[,t,i]  <- AUC_height[,t,i] * AUC_width[,t,i]
    }
    
    for (i in 1:n_samps) { 
      AUC_samps[i,t] = sum(  AUC_height_x_width[i,t,] , na.rm = TRUE);
    }
    
    AUC$mean[t] <- mean( AUC_samps[ ,t] )
    AUC$sd[t] <- sd( AUC_samps[ ,t] )
    AUC$lower[t] <- sort( AUC_samps[ ,t] )[n_samps*0.025]
    AUC$median[t] <- median( AUC_samps[ ,t] )
    AUC$upper[t] <- sort( AUC_samps[ ,t] )[n_samps*0.975]
    
  }
  
  AUC
  
  # calculate the optimal threshold for each test (using se and sp summaries from stan model output)
  
  # according to youden index
  youden <- array(data = NA, dim = c(data$n_tests,max_category))
  max_youden <- c()
  optimal_threshold_youden <- c()
  
  for (t in 1:data$n_tests) { 
    youden[t,] <- stan_youden[,5][((t-1)*max_category+1):(t*max_category)]
    max_youden[t] <- max(youden[t,], na.rm = TRUE)
    optimal_threshold_youden[t] <-    match(max_youden[t], youden[t,])
  }
  
  
  # according to weighted youden index
  w_youden <- array(data = NA, dim = c(data$n_tests,max_category))
  max_w_youden <- c()
  optimal_threshold_w_youden <- c()
  
  for (t in 1:data$n_tests) { 
    w_youden[t,] <- stan_youden_weighted[,5][((t-1)*max_category+1):(t*max_category)]
    max_w_youden[t] <- max(w_youden[t,], na.rm = TRUE)
    optimal_threshold_w_youden[t] <-    match(max_w_youden[t], w_youden[t,])
  }
  
  # according to DOR
  dor <- array(data = NA, dim = c(data$n_tests,max_category))
  max_dor <- c()
  optimal_threshold_dor <- c()
  
  for (t in 1:data$n_tests) { 
    dor[t,] <- stan_DOR[,5][((t-1)*max_category+1):(t*max_category)]
    max_dor[t] <- max(dor[t,], na.rm = TRUE)
    optimal_threshold_dor[t] <-    match(max_dor[t], dor[t,])
  }
  
  # Comparative accuracy Se and Sp
  # Se and Sp  - calculated at the optimal threshold for each test
  Se_youden_samps<-list()
  Sp_youden_samps<-list()
  Se_w_youden_samps<-list()
  Sp_w_youden_samps<-list()
  Se_dor_samps<-list()
  Sp_dor_samps<-list()
  
  
  Se_youden <- tibble(Test = seq(1,data$n_tests),
                      Optimal_threshold = optimal_threshold_youden,
                      mean = rep(NA, data$n_tests),  sd =  rep(NA, data$n_tests),
                      lower  = rep(NA, data$n_tests), median  = rep(NA, data$n_tests), upper  = rep(NA, data$n_tests))
  Sp_youden <- Se_youden
  Se_w_youden <- tibble(Test = seq(1,data$n_tests),
                        Optimal_threshold = optimal_threshold_w_youden,
                        mean = rep(NA, data$n_tests),  sd =  rep(NA, data$n_tests),
                        lower  = rep(NA, data$n_tests), median  = rep(NA, data$n_tests), upper  = rep(NA, data$n_tests))
  Sp_w_youden <- Se_w_youden
  Se_dor <- tibble(Test = seq(1,data$n_tests),
                   Optimal_threshold = optimal_threshold_dor,
                   mean = rep(NA, data$n_tests),  sd =  rep(NA, data$n_tests),
                   lower  = rep(NA, data$n_tests), median  = rep(NA, data$n_tests), upper  = rep(NA, data$n_tests))
  Sp_dor <- Se_dor
  
  
  for (t in 1:data$n_tests) {
    # youden
    Se_youden_samps[[t]] <- Se[,t, optimal_threshold_youden[t]];
    Sp_youden_samps[[t]] <- Sp[,t, optimal_threshold_youden[t]];
    
    Se_youden$mean[t] <- mean( Se_youden_samps[[t]] )
    Se_youden$sd[t] <- sd( Se_youden_samps[[t]] )
    Se_youden$lower[t] <- sort( Se_youden_samps[[t]] )[n_samps*0.025]
    Se_youden$median[t] <- median( Se_youden_samps[[t]] )
    Se_youden$upper[t] <- sort( Se_youden_samps[[t]] )[n_samps*0.975]
    
    Sp_youden$mean[t] <- mean( Sp_youden_samps[[t]] )
    Sp_youden$sd[t] <- sd( Sp_youden_samps[[t]] )
    Sp_youden$lower[t] <- sort( Sp_youden_samps[[t]] )[n_samps*0.025]
    Sp_youden$median[t] <- median( Sp_youden_samps[[t]] )
    Sp_youden$upper[t] <- sort( Sp_youden_samps[[t]] )[n_samps*0.975]
    
    # weighted youden
    Se_w_youden_samps[[t]] <- Se[,t, optimal_threshold_w_youden[t]];
    Sp_w_youden_samps[[t]] <- Sp[,t, optimal_threshold_w_youden[t]];
    
    Se_w_youden$mean[t] <- mean( Se_w_youden_samps[[t]] )
    Se_w_youden$sd[t] <- sd( Se_w_youden_samps[[t]] )
    Se_w_youden$lower[t] <- sort( Se_w_youden_samps[[t]] )[n_samps*0.025]
    Se_w_youden$median[t] <- median( Se_w_youden_samps[[t]] )
    Se_w_youden$upper[t] <- sort( Se_w_youden_samps[[t]] )[n_samps*0.975]
    
    Sp_w_youden$mean[t] <- mean( Sp_w_youden_samps[[t]] )
    Sp_w_youden$sd[t] <- sd( Sp_w_youden_samps[[t]] )
    Sp_w_youden$lower[t] <- sort( Sp_w_youden_samps[[t]] )[n_samps*0.025]
    Sp_w_youden$median[t] <- median( Sp_w_youden_samps[[t]] )
    Sp_w_youden$upper[t] <- sort( Sp_w_youden_samps[[t]] )[n_samps*0.975]
    
    # dor
    Se_dor_samps[[t]] <- Se[,t, optimal_threshold_dor[t]];
    Sp_dor_samps[[t]] <- Sp[,t, optimal_threshold_dor[t]];
    
    Se_dor$mean[t] <- mean( Se_dor_samps[[t]] )
    Se_dor$sd[t] <- sd( Se_dor_samps[[t]] )
    Se_dor$lower[t] <- sort( Se_dor_samps[[t]] )[n_samps*0.025]
    Se_dor$median[t] <- median( Se_dor_samps[[t]] )
    Se_dor$upper[t] <- sort( Se_dor_samps[[t]] )[n_samps*0.975]
    
    Sp_dor$mean[t] <- mean( Sp_dor_samps[[t]] )
    Sp_dor$sd[t] <- sd( Sp_dor_samps[[t]] )
    Sp_dor$lower[t] <- sort( Sp_dor_samps[[t]] )[n_samps*0.025]
    Sp_dor$median[t] <- median( Sp_dor_samps[[t]] )
    Sp_dor$upper[t] <- sort( Sp_dor_samps[[t]] )[n_samps*0.975]
  }
  
  Se_youden
  
  # calculate pairwise accuracy differences and ratios for each test/threshold combo - total = 5C2 = 10 for Se and 10 for Sp
  
  diff_Se_youden_samps <- list()
  diff_Sp_youden_samps <- list()
  ratio_Se_youden_samps <- list()
  ratio_Sp_youden_samps <- list()
  
  diff_Se_w_youden_samps <- list()
  diff_Sp_w_youden_samps <- list()
  ratio_Se_w_youden_samps <- list()
  ratio_Sp_w_youden_samps <- list()
  
  diff_Se_dor_samps <- list()
  diff_Sp_dor_samps <- list()
  ratio_Se_dor_samps <- list()
  ratio_Sp_dor_samps <- list()
  
  diff_Se_youden <- tibble( Comparison = rep(NA, (choose(data$n_tests,2))),
                            mean = rep(NA, (choose(data$n_tests,2))), 
                            sd =  rep(NA, choose(data$n_tests,2)),
                            lower  = rep(NA, choose(data$n_tests,2)), 
                            median  = rep(NA,choose(data$n_tests,2)),
                            upper  = rep(NA, choose(data$n_tests,2)))
  
  for (i in 1:4) {
     diff_Se_youden$Comparison[i] <- paste("test", 1,"vs",i+1)
  }
  for (i in 5:7) { 
    diff_Se_youden$Comparison[i] <- paste("test", 2,"vs",i-2)
  }
  for (i in 8:9) { 
    diff_Se_youden$Comparison[i] <- paste("test", 2,"vs",i-4)
    }
  diff_Se_youden$Comparison[10] <- paste("test", 4,"vs",5)
  
  diff_Sp_youden <- diff_Se_youden
  diff_Se_w_youden <- diff_Se_youden
  diff_Sp_w_youden <- diff_Se_youden
  diff_Se_dor <- diff_Se_youden
  diff_Sp_dor <- diff_Se_dor
  
  ratio_Se_youden <- tibble( Comparison = rep(NA, (choose(data$n_tests,2))),
                            mean = rep(NA, (choose(data$n_tests,2))), 
                            sd =  rep(NA, choose(data$n_tests,2)),
                            lower  = rep(NA, choose(data$n_tests,2)), 
                            median  = rep(NA,choose(data$n_tests,2)),
                            upper  = rep(NA, choose(data$n_tests,2)))
  
  for (i in 1:4) {
    ratio_Se_youden$Comparison[i] <- paste("test", 1,"vs",i+1)
  }
  for (i in 5:7) { 
    ratio_Se_youden$Comparison[i] <- paste("test", 2,"vs",i-2)
  }
  for (i in 8:9) { 
    ratio_Se_youden$Comparison[i] <- paste("test", 2,"vs",i-4)
  }
  ratio_Se_youden$Comparison[10] <- paste("test", 4,"vs",5)
  
  ratio_Sp_youden <- ratio_Se_youden
  ratio_Se_w_youden <- ratio_Se_youden
  ratio_Sp_w_youden <- ratio_Se_youden
  ratio_Se_dor <- ratio_Se_youden
  ratio_Sp_dor <- ratio_Se_dor
  
  
  # using youden
  for (t in 1:4) {
    diff_Se_youden_samps[[t]] = Se_youden_samps[[1]] - Se_youden_samps[[t+1]];
    diff_Sp_youden_samps[[t]] = Sp_youden_samps[[1]] - Sp_youden_samps[[t+1]];
    ratio_Se_youden_samps[[t]] = Se_youden_samps[[1]] /  Se_youden_samps[[t+1]];
    ratio_Sp_youden_samps[[t]] = Sp_youden_samps[[1]] /  Sp_youden_samps[[t+1]];
  }
  
  for (t in 5:7) {
    diff_Se_youden_samps[[t]] =  Se_youden_samps[[2]] -  Se_youden_samps[[t-2]];
    diff_Sp_youden_samps[[t]] =  Sp_youden_samps[[2]] -  Sp_youden_samps[[t-2]];
    ratio_Se_youden_samps[[t]] = Se_youden_samps[[2]] /  Se_youden_samps[[t-2]];
    ratio_Sp_youden_samps[[t]] = Sp_youden_samps[[2]] /  Sp_youden_samps[[t-2]];
  }
  
  for (t in 8:9) {
    diff_Se_youden_samps[[t]] =  Se_youden_samps[[3]] -  Se_youden_samps[[t-4]];
    diff_Sp_youden_samps[[t]] =  Sp_youden_samps[[3]] -  Sp_youden_samps[[t-4]];
    ratio_Se_youden_samps[[t]] = Se_youden_samps[[3]] /  Se_youden_samps[[t-4]];
    ratio_Sp_youden_samps[[t]] = Sp_youden_samps[[3]] /  Sp_youden_samps[[t-4]];
  }
  
  diff_Se_youden_samps[[10]] =  Se_youden_samps[[4]] -  Se_youden_samps[[5]];
  diff_Sp_youden_samps[[10]] =  Sp_youden_samps[[4]] -  Sp_youden_samps[[5]];
  ratio_Se_youden_samps[[10]] = Se_youden_samps[[4]] /  Se_youden_samps[[5]];
  ratio_Sp_youden_samps[[10]] = Sp_youden_samps[[4]] /  Sp_youden_samps[[5]];
  
  for (t in 1:nrow(diff_Se_youden)) { 
    diff_Se_youden$mean[t] <- mean( diff_Se_youden_samps[[t]] )
    diff_Se_youden$sd[t] <- sd( diff_Se_youden_samps[[t]] )
    diff_Se_youden$lower[t] <- sort( diff_Se_youden_samps[[t]] )[n_samps*0.025]
    diff_Se_youden$median[t] <- median( diff_Se_youden_samps[[t]] )
    diff_Se_youden$upper[t] <- sort( diff_Se_youden_samps[[t]] )[n_samps*0.975]
    
    diff_Sp_youden$mean[t] <- mean( diff_Sp_youden_samps[[t]] )
    diff_Sp_youden$sd[t] <- sd( diff_Sp_youden_samps[[t]] )
    diff_Sp_youden$lower[t] <- sort( diff_Sp_youden_samps[[t]] )[n_samps*0.025]
    diff_Sp_youden$median[t] <- median( diff_Sp_youden_samps[[t]] )
    diff_Sp_youden$upper[t] <- sort( diff_Sp_youden_samps[[t]] )[n_samps*0.975]
    
    ratio_Se_youden$mean[t] <- mean( ratio_Se_youden_samps[[t]] )
    ratio_Se_youden$sd[t] <- sd( ratio_Se_youden_samps[[t]] )
    ratio_Se_youden$lower[t] <- sort( ratio_Se_youden_samps[[t]] )[n_samps*0.025]
    ratio_Se_youden$median[t] <- median( ratio_Se_youden_samps[[t]] )
    ratio_Se_youden$upper[t] <- sort( ratio_Se_youden_samps[[t]] )[n_samps*0.975]
    
    ratio_Sp_youden$mean[t] <- mean( ratio_Sp_youden_samps[[t]] )
    ratio_Sp_youden$sd[t] <- sd( ratio_Sp_youden_samps[[t]] )
    ratio_Sp_youden$lower[t] <- sort( ratio_Sp_youden_samps[[t]] )[n_samps*0.025]
    ratio_Sp_youden$median[t] <- median( ratio_Sp_youden_samps[[t]] )
    ratio_Sp_youden$upper[t] <- sort( ratio_Sp_youden_samps[[t]] )[n_samps*0.975]
    
  }
  
  
  # using weighted youden
  for (t in 1:4) {
    diff_Se_w_youden_samps[[t]] = Se_w_youden_samps[[1]] - Se_w_youden_samps[[t+1]];
    diff_Sp_w_youden_samps[[t]] = Sp_w_youden_samps[[1]] - Sp_w_youden_samps[[t+1]];
    ratio_Se_w_youden_samps[[t]] = Se_w_youden_samps[[1]] /  Se_w_youden_samps[[t+1]];
    ratio_Sp_w_youden_samps[[t]] = Sp_w_youden_samps[[1]] /  Sp_w_youden_samps[[t+1]];
  }
  
  for (t in 5:7) {
    diff_Se_w_youden_samps[[t]] =  Se_w_youden_samps[[2]] -  Se_w_youden_samps[[t-2]];
    diff_Sp_w_youden_samps[[t]] =  Sp_w_youden_samps[[2]] -  Sp_w_youden_samps[[t-2]];
    ratio_Se_w_youden_samps[[t]] = Se_w_youden_samps[[2]] /  Se_w_youden_samps[[t-2]];
    ratio_Sp_w_youden_samps[[t]] = Sp_w_youden_samps[[2]] /  Sp_w_youden_samps[[t-2]];
  }
  
  for (t in 8:9) {
    diff_Se_w_youden_samps[[t]] =  Se_w_youden_samps[[3]] -  Se_w_youden_samps[[t-4]];
    diff_Sp_w_youden_samps[[t]] =  Sp_w_youden_samps[[3]] -  Sp_w_youden_samps[[t-4]];
    ratio_Se_w_youden_samps[[t]] = Se_w_youden_samps[[3]] /  Se_w_youden_samps[[t-4]];
    ratio_Sp_w_youden_samps[[t]] = Sp_w_youden_samps[[3]] /  Sp_w_youden_samps[[t-4]];
  }
  
  diff_Se_w_youden_samps[[10]] =  Se_w_youden_samps[[4]] -  Se_w_youden_samps[[5]];
  diff_Sp_w_youden_samps[[10]] =  Sp_w_youden_samps[[4]] -  Sp_w_youden_samps[[5]];
  ratio_Se_w_youden_samps[[10]] = Se_w_youden_samps[[4]] /  Se_w_youden_samps[[5]];
  ratio_Sp_w_youden_samps[[10]] = Sp_w_youden_samps[[4]] /  Sp_w_youden_samps[[5]];
  
  
  for (t in 1:nrow(diff_Se_youden)) { 
    diff_Se_w_youden$mean[t] <- mean( diff_Se_w_youden_samps[[t]] )
    diff_Se_w_youden$sd[t] <- sd( diff_Se_w_youden_samps[[t]] )
    diff_Se_w_youden$lower[t] <- sort( diff_Se_w_youden_samps[[t]] )[n_samps*0.025]
    diff_Se_w_youden$median[t] <- median( diff_Se_w_youden_samps[[t]] )
    diff_Se_w_youden$upper[t] <- sort( diff_Se_w_youden_samps[[t]] )[n_samps*0.975]
    
    diff_Sp_w_youden$mean[t] <- mean( diff_Sp_w_youden_samps[[t]] )
    diff_Sp_w_youden$sd[t] <- sd( diff_Sp_w_youden_samps[[t]] )
    diff_Sp_w_youden$lower[t] <- sort( diff_Sp_w_youden_samps[[t]] )[n_samps*0.025]
    diff_Sp_w_youden$median[t] <- median( diff_Sp_w_youden_samps[[t]] )
    diff_Sp_w_youden$upper[t] <- sort( diff_Sp_w_youden_samps[[t]] )[n_samps*0.975]
    
    ratio_Se_w_youden$mean[t] <- mean( ratio_Se_w_youden_samps[[t]] )
    ratio_Se_w_youden$sd[t] <- sd( ratio_Se_w_youden_samps[[t]] )
    ratio_Se_w_youden$lower[t] <- sort( ratio_Se_w_youden_samps[[t]] )[n_samps*0.025]
    ratio_Se_w_youden$median[t] <- median( ratio_Se_w_youden_samps[[t]] )
    ratio_Se_w_youden$upper[t] <- sort( ratio_Se_w_youden_samps[[t]] )[n_samps*0.975]
    
    ratio_Sp_w_youden$mean[t] <- mean( ratio_Sp_w_youden_samps[[t]] )
    ratio_Sp_w_youden$sd[t] <- sd( ratio_Sp_w_youden_samps[[t]] )
    ratio_Sp_w_youden$lower[t] <- sort( ratio_Sp_w_youden_samps[[t]] )[n_samps*0.025]
    ratio_Sp_w_youden$median[t] <- median( ratio_Sp_w_youden_samps[[t]] )
    ratio_Sp_w_youden$upper[t] <- sort( ratio_Sp_w_youden_samps[[t]] )[n_samps*0.975]
    
  }
  
  
  # using dor
  for (t in 1:4) {
    diff_Se_dor_samps[[t]] = Se_dor_samps[[1]] - Se_dor_samps[[t+1]];
    diff_Sp_dor_samps[[t]] = Sp_dor_samps[[1]] - Sp_dor_samps[[t+1]];
    ratio_Se_dor_samps[[t]] = Se_dor_samps[[1]] /  Se_dor_samps[[t+1]];
    ratio_Sp_dor_samps[[t]] = Sp_dor_samps[[1]] /  Sp_dor_samps[[t+1]];
  }
  
  for (t in 5:7) {
    diff_Se_dor_samps[[t]] =  Se_dor_samps[[2]] -  Se_dor_samps[[t-2]];
    diff_Sp_dor_samps[[t]] =  Sp_dor_samps[[2]] -  Sp_dor_samps[[t-2]];
    ratio_Se_dor_samps[[t]] = Se_dor_samps[[2]] /  Se_dor_samps[[t-2]];
    ratio_Sp_dor_samps[[t]] = Sp_dor_samps[[2]] /  Sp_dor_samps[[t-2]];
  }
  
  for (t in 8:9) {
    diff_Se_dor_samps[[t]] =  Se_dor_samps[[3]] -  Se_dor_samps[[t-4]];
    diff_Sp_dor_samps[[t]] =  Sp_dor_samps[[3]] -  Sp_dor_samps[[t-4]];
    ratio_Se_dor_samps[[t]] = Se_dor_samps[[3]] /  Se_dor_samps[[t-4]];
    ratio_Sp_dor_samps[[t]] = Sp_dor_samps[[3]] /  Sp_dor_samps[[t-4]];
  }
  
  diff_Se_dor_samps[[10]] =  Se_dor_samps[[4]] -  Se_dor_samps[[5]];
  diff_Sp_dor_samps[[10]] =  Sp_dor_samps[[4]] -  Sp_dor_samps[[5]];
  ratio_Se_dor_samps[[10]] = Se_dor_samps[[4]] /  Se_dor_samps[[5]];
  ratio_Sp_dor_samps[[10]] = Sp_dor_samps[[4]] /  Sp_dor_samps[[5]];
  
  for (t in 1:nrow(diff_Se_youden)) { 
    diff_Se_dor$mean[t] <- mean( diff_Se_dor_samps[[t]] )
    diff_Se_dor$sd[t] <- sd( diff_Se_dor_samps[[t]] )
    diff_Se_dor$lower[t] <- sort( diff_Se_dor_samps[[t]] )[n_samps*0.025]
    diff_Se_dor$median[t] <- median( diff_Se_dor_samps[[t]] )
    diff_Se_dor$upper[t] <- sort( diff_Se_dor_samps[[t]] )[n_samps*0.975]
    
    diff_Sp_dor$mean[t] <- mean( diff_Sp_dor_samps[[t]] )
    diff_Sp_dor$sd[t] <- sd( diff_Sp_dor_samps[[t]] )
    diff_Sp_dor$lower[t] <- sort( diff_Sp_dor_samps[[t]] )[n_samps*0.025]
    diff_Sp_dor$median[t] <- median( diff_Sp_dor_samps[[t]] )
    diff_Sp_dor$upper[t] <- sort( diff_Sp_dor_samps[[t]] )[n_samps*0.975]
    
    ratio_Se_dor$mean[t] <- mean( ratio_Se_dor_samps[[t]] )
    ratio_Se_dor$sd[t] <- sd( ratio_Se_dor_samps[[t]] )
    ratio_Se_dor$lower[t] <- sort( ratio_Se_dor_samps[[t]] )[n_samps*0.025]
    ratio_Se_dor$median[t] <- median( ratio_Se_dor_samps[[t]] )
    ratio_Se_dor$upper[t] <- sort( ratio_Se_dor_samps[[t]] )[n_samps*0.975]
    
    ratio_Sp_dor$mean[t] <- mean( ratio_Sp_dor_samps[[t]] )
    ratio_Sp_dor$sd[t] <- sd( ratio_Sp_dor_samps[[t]] )
    ratio_Sp_dor$lower[t] <- sort( ratio_Sp_dor_samps[[t]] )[n_samps*0.025]
    ratio_Sp_dor$median[t] <- median( ratio_Sp_dor_samps[[t]] )
    ratio_Sp_dor$upper[t] <- sort( ratio_Sp_dor_samps[[t]] )[n_samps*0.975]
    
  }
  
  
  
  
  stan_list <- list(stan_se = stan_se,
                    stan_sp = stan_sp,
                    stan_se_pred = stan_se_pred,
                    stan_sp_pred = stan_sp_pred,
                    stan_ppv = stan_ppv,
                    stan_npv = stan_npv,
                    stan_LRpos = stan_LRpos,
                    stan_LRneg = stan_LRneg,
                    stan_var_cov = stan_var_cov,
                    stan_corr = stan_corr,
                    stan_youden = stan_youden,
                    stan_youden_weighted = stan_youden_weighted,
                    stan_DOR = stan_DOR)
  
  out_list <- list(AUC = AUC, 
                   Se_youden = Se_youden,
                   Sp_youden = Sp_youden,
                   Se_w_youden = Se_w_youden,
                   Sp_w_youden = Sp_w_youden,
                   Se_dor = Se_dor,
                   Sp_dor = Sp_dor,
                   diff_Se_youden = diff_Se_youden,
                   diff_Sp_youden = diff_Sp_youden,
                   ratio_Se_youden = ratio_Se_youden,
                   ratio_Sp_youden = ratio_Sp_youden,
                   diff_Se_w_youden = diff_Se_w_youden,
                   diff_Sp_w_youden = diff_Sp_w_youden,
                   ratio_Se_w_youden = ratio_Se_w_youden,
                   ratio_Sp_w_youden = ratio_Sp_w_youden,
                   diff_Se_dor = diff_Se_dor,
                   diff_Sp_dor = diff_Sp_dor,
                   ratio_Se_dor = ratio_Se_dor,
                   ratio_Sp_dor = ratio_Sp_dor)
                  
                  
  my_list <- list(model_output = model_rstan,
                  stan_list = stan_list, 
                  out_list = out_list)
  
  
  return(my_list)
  
  
}







# Functions for K-fold cross-validation  ----------------------------------------------------------------------------------------
#  found on GitHub from:   https://github.com/stan-dev/stancon_talks/blob/master/2017/Contributed-Talks/07_nicenboim/kfold.Rmd

# The following function can run all the chains of all the folds of the model in parallel:
stan_kfold <- function(file, list_of_datas, chains, cores,...){
  badRhat <- 1.1
  K <- length(list_of_datas)
  model <- stan_model(file=file)
  # First parallelize all chains:
  sflist <-  pbmclapply(1:(K*chains), mc.cores = cores, 
                       function(i){
                         # Fold number:
                         k <- round((i+1) / chains)
                         s <- sampling(model, data = list_of_datas[[k]], 
                                       chains = 1, chain_id = i,  ...)
                         return(s)
                       })
  # Then merge the K * chains to create K stanfits:
  stanfit <- list()
  for(k in 1:K){
    inchains <- (chains*k - 2):(chains*k)
    # Merge `chains` of each fold
    stanfit[[k]] <- sflist2stanfit(sflist[inchains])
  }  
  return(stanfit)
  # 
  # holdout <- lapply(list_of_datas, '[[', "holdout")
  # 
  # # We extract all the held_out log_lik of all the folds
  # log_lik_ab <- extract_log_lik_K(stanfit, holdout, "log_lik")
  # save(log_lik_ab, file = "log_lik_ab.Rda")
  # str(log_lik_ab)
  # 
  # length <- length(log_lik_ab[1,][log_lik_ab[1,] != "NaN" & !is.na(log_lik_ab[1,]) ])
  # loglik2 <- array(data = log_lik_ab[log_lik_ab != "NaN" & !is.na(log_lik_ab) ], dim = c(750, length ))
  # 
  # kfold_model <- kfold(loglik2)
  # 
  # my_list <- list(fits = stanfit, 
  #                 log_lik = loglik2,
  #                 kfold = kfold_model)
  # 
  # return(my_list)
  
}

# Wrapper function to extract the log_lik of the held-out data, given a list of stanfits, and a list which indicates with 1 and 0 whether the observation was held out or not:
extract_log_lik_K <- function(list_of_stanfits, list_of_holdout, ...){
  K <- length(list_of_stanfits)
  list_of_log_liks <- plyr::llply(1:K, function(k){
    extract_log_lik(list_of_stanfits[[k]], merge_chains = TRUE , ...)
  })
  # `log_lik_heldout` will include the loglike of all the held out data of all the folds.
  # We define `log_lik_heldout` as a (samples x N_obs) matrix
  # (similar to each log_lik matrix)
  log_lik_heldout <- list_of_log_liks[[1]] * NA
  for(k in 1:K){
    log_lik <- list_of_log_liks[[k]]
    samples <- dim(log_lik)[1] 
    N_obs <- dim(log_lik)[2]
    # This is a matrix with the same size as log_lik_heldout
    # with 1 if the data was held out in the fold k
    heldout <- matrix(rep(list_of_holdout[[k]], each = samples), nrow = samples)
    # Sanity check that the previous log_lik is not being overwritten:
    if(any(!is.na(log_lik_heldout[heldout==1]))){
      warning("Heldout log_lik has been overwritten!!!!")
    }
    # We save here the log_lik of the fold k in the matrix:
    log_lik_heldout[heldout==1] <- log_lik[heldout==1]
  }
  return(log_lik_heldout)
}


kfold <- function(log_lik_heldout)  {
  library(matrixStats)
  logColMeansExp <- function(x) {
    # should be more stable than log(colMeans(exp(x)))
    S <- nrow(x)
    colLogSumExps(x) - log(S)
  }
  # See equation (20) of @VehtariEtAl2016
  pointwise <-  matrix(logColMeansExp(log_lik_heldout), ncol= 1)
  colnames(pointwise) <- "elpd"
  # See equation (21) of @VehtariEtAl2016
  elpd_kfold <- sum(pointwise)
  se_elpd_kfold <-  sqrt(ncol(log_lik_heldout) * var(pointwise))
  
  out <- list(
    pointwise = pointwise,
    elpd_kfold = elpd_kfold,
    se_elpd_kfold = se_elpd_kfold)
  out
  #structure(out, class = "loo")
}







