

#  -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plots   -------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
# load models
mod <- readRDS("./outputs/model_1A_prior_set_1_output.rds") # CS NMA model (main model chosen according to full K-fold cross-validation) 




# Plot (1) ---------------------------------------------------------------------------------------------------------------------
# sROC plot of full ROC curves w/ points for Se/Sp for each test and threshold combo, using posterior medians 

AUC <- c()
AUC_L<- c()
AUC_U<- c()

for ( i in 1:5) { 
  AUC[(5*(i-1)+1):(5*i)] <- rep(mod$out_list$AUC$median[i], 5)
  AUC_L[(5*(i-1)+1):(5*i)] <- rep(mod$out_list$AUC$lower[i], 5)
  AUC_U[(5*(i-1)+1):(5*i)] <- rep(mod$out_list$AUC$upper[i], 5)
}

round(rep(mod$out_list$AUC$median, 5),2)

data_plot_1 <- tibble(Sensitivity = mod$stan_list$stan_se[,5], 
                      Sensitivity_L =  mod$stan_list$stan_se[,4],
                      Sensitivity_U =  mod$stan_list$stan_se[,6],
                      Specificity = mod$stan_list$stan_sp[,5],
                      Specificity_L = mod$stan_list$stan_sp[,4],
                      Specificity_U = mod$stan_list$stan_sp[,6],
                      AUC = round(AUC, 2),
                      AUC_L =  round(AUC_L, 2),
                      AUC_U = round(AUC_U, 2), 
                      AUC_lab = paste0("AUC=", ".",AUC*100," [",".",AUC_L*100,", ",".",AUC_U*100,"]"))

data_plot_1 <- data_plot_1 %>% 
  filter(!is.na(Sensitivity)) %>%
  mutate(test_threshold_num = seq(from = 1, to = 22))

# using data.table
setDT(X2); setDT(data_plot_1)  
data_plot_1 <- X2[data_plot_1, mult = "first", on = "test_threshold_num", nomatch=0L] 

data_plot_1 <- data_plot_1 %>% 
  rename(Threshold = Category,  Test = test) %>%
  select(Test, 
         test_num,
         Threshold,
         threshold_num,
         Sensitivity, Sensitivity_L, Sensitivity_U,
         Specificity, Specificity_L, Specificity_U,
         AUC, AUC_L, AUC_U,AUC_lab) %>%
  mutate(Test = factor(Test), 
         Threshold = factor(Threshold), 
         Sensitivity = as.numeric(Sensitivity), 
         Specificity = as.numeric(Specificity))

data_plot_1 <- tibble(data_plot_1)

# plot 
require(RColorBrewer)
require(ggrepel)


g <- ggplot(data = data_plot_1, 
            aes(y=Sensitivity, 
                x = 1 - Specificity, 
                colour = Test)) + 
  geom_point(size = 5)  +     
  geom_line() + 
  geom_text_repel(data = filter(data_plot_1, threshold_num  == 3),
                  aes(y=Sensitivity, 
                      x = 1 - Specificity,
                      label = AUC_lab),
                  nudge_x = 0.2,
                  nudge_y = - 0.05,
                  inherit.aes = FALSE) + 
  theme_bw() + 
  scale_x_continuous(breaks = seq(0,1,0.1), limits = c(0,1))  + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1))  + 
  scale_color_brewer(palette="Set1") + 
  coord_fixed()
g


#  produce high quality TIFF for publication
tiff("./plots/fig_1.tif",units = "in", width = 10, height=7, res=800, compression = "lzw") 
g
dev.off()




# Plot (2) --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# sROC plot of summary Se and Sp with corresponding credible and prediction 
# intervals for each test, only at the thresholds chosen according to Youden index 


## credible region
cred_1 <- list()

for (t in 1:data$n_tests) { 
  for (c in 1:data$n_thresholds[t]) { 
    cred_1[[c + sum(data$n_thresholds[1:t-1])]] <- tibble(y = (rstan::extract(mod$model_output, pars = "lSe")$lSe[,t,c]) , 
                                                          x = (rstan::extract(mod$model_output, pars = "lSp")$lSp[,t,c]))
  }
}

require(data.table)

cred <- rbindlist(cred_1, idcol = TRUE)



cred2 <- cred %>%
  mutate(test_threshold_num = .id)

setDT(X2); setDT(cred2)  
cred2 <- X2[cred2, mult = "first", on = "test_threshold_num", nomatch=0L] 

# Se and Sp only at the optimal thr (according to youden) 
ranks <- mod$out_list$Se_youden$Optimal_threshold # youden
#ranks <- mod$out_list$Se_w_youden$Optimal_threshold # weighted youden
#ranks <- mod$out_list$Se_dor$Optimal_threshold # DOR

cred2 <- cred2 %>%
  rename(Threshold = Category,  Test = test) %>%
  select(Test, 
         test_num,
         Threshold,
         threshold_num,
         x, y) %>%
  mutate(Test = factor(Test), 
         Threshold = factor(Threshold)) %>%
  # use thresholds according to Youden index - change accordingly if using weighted youden or DOR
  filter(test_num == 1 & threshold_num == ranks[1] |
           test_num == 2 & threshold_num == ranks[2] |
           test_num == 3 & threshold_num == ranks[3] |
           test_num == 4 & threshold_num == ranks[4] |
           test_num == 5 & threshold_num == ranks[5] )





# in inv_probit space
g <- ggplot(data = cred2, aes(x = x, y = y, colour = Test))  + 
  stat_ellipse()  

g

# Get ellipse coordinates from plot
pb <-  ggplot_build(g)
el = pb$data[[1]][c("x","y", "group")]

el$group
levels(el$group) <- list(1 = "hello")

levels(cred2$Test)[1]

el2 <- el %>% 
  mutate(group = case_when(group == "1" ~ levels(cred2$Test)[1],
                           group == "2" ~ levels(cred2$Test)[2],
                           group == "3" ~ levels(cred2$Test)[3],
                           group == "4" ~ levels(cred2$Test)[4],
                           group == "5" ~ levels(cred2$Test)[5]))

el2


credible_region <- tibble(x = plogis(el$x), y = plogis(el$y), Test = el2$group)


credible_region

g <- ggplot(data = credible_region, aes(x = x, y = y, colour = Test))  + 
  geom_polygon(data = credible_region, aes(x = 1  - x, y = y, colour = Test), alpha=0.05, size=0.4)  + 
  xlim(0,1) + 
  ylim(0,1)

g


## prediction region
pred_1 <- list()

for (t in 1:5) { 
  for (c in 1:data$n_thresholds[t]) { 
    pred_1[[c + sum(data$n_thresholds[1:t-1])]] <- tibble(y = (rstan::extract(mod$model_output, pars = "lSe_pred")$lSe_pred[,t,c]), 
                                                          x = (rstan::extract(mod$model_output, pars = "lSp_pred")$lSp_pred[,t,c]))
  }
}



pred <- rbindlist(pred_1, idcol = TRUE)


pred2 <- pred %>%
  mutate(test_threshold_num = .id)

setDT(X2); setDT(pred2)  
pred2 <- X2[pred2, mult = "first", on = "test_threshold_num", nomatch=0L] 

# Se and Sp only at the optimal thr (according to youden) 
ranks <- mod$out_list$Se_youden$Optimal_threshold # youden
#ranks <- mod$out_list$Se_w_youden$Optimal_threshold # weighted youden
#ranks <- mod$out_list$Se_dor$Optimal_threshold # DOR


pred2 <- pred2 %>%
  rename(Threshold = Category,  Test = test) %>%
  select(Test, 
         test_num,
         Threshold,
         threshold_num,
         x, y) %>%
  mutate(Test = factor(Test), 
         Threshold = factor(Threshold)) %>%
  # use thresholds according to Youden index - change accordingly if using weighted youden or DOR
  filter(test_num == 1 & threshold_num == ranks[1] |
           test_num == 2 & threshold_num == ranks[2] |
           test_num == 3 & threshold_num == ranks[3] |
           test_num == 4 & threshold_num == ranks[4] |
           test_num == 5 & threshold_num == ranks[5] )




# in inv_probit space
g <- ggplot(data = pred2, aes(x = x, y = y, colour = Test))  + 
  stat_ellipse()  

g

# Get ellipse coordinates from plot
pb <-  ggplot_build(g)
el = pb$data[[1]][c("x","y", "group")]


el2 <- el %>% 
  mutate(group = case_when(group == "1" ~ levels(pred2$Test)[1],
                           group == "2" ~ levels(pred2$Test)[2],
                           group == "3" ~ levels(pred2$Test)[3],
                           group == "4" ~ levels(pred2$Test)[4],
                           group == "5" ~ levels(pred2$Test)[5]))

pred_region <- tibble(x = plogis(el$x), y = plogis(el$y), Test = el2$group)


g <- ggplot(data = pred_region, aes(x = x, y = y, colour = Test))  + 
  geom_polygon(data = pred_region, aes(x = 1  - x, y = y, colour = Test), alpha=0.05, size=0.4)  + 
  xlim(0,1) + 
  ylim(0,1)

g



data_plot_1_subset <- data_plot_1   %>%     filter(test_num == 1 & threshold_num == ranks[1] |
                                                   test_num == 2 & threshold_num == ranks[2] |
                                                   test_num == 3 & threshold_num == ranks[3] |
                                                   test_num == 4 & threshold_num == ranks[4] |
                                                   test_num == 5 & threshold_num == ranks[5] )

g <- ggplot(data = data_plot_1_subset, 
            aes(y=Sensitivity, 
                x = 1 - Specificity, 
                colour = Test)) + 
  geom_point(size = 5)  +     
  geom_line() + 
  geom_path(data = pred_region, 
            aes(x= 1 - x, y= y, colour = Test), 
            linetype = 2, size = 0.4, inherit.aes = F) +                         # prediction region
  geom_polygon(data = credible_region, 
               aes(x= 1 - x, y= y, colour = Test), 
               alpha=0.05, size=0.4, linetype = 2, inherit.aes = F) +            # conf region
  theme_bw() + 
  scale_x_continuous(breaks = seq(0,1,0.1), limits = c(0,1))  + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1))  + 
  scale_color_brewer(palette="Set1")

g





tiff("./plots/fig_2.tif",units = "in", width = 7*1.25, height=5.2*1.25, res=500, compression = "lzw")
g
dev.off()





# Plots for follow-up vs accuracy -------------------------------------------------------------------------------------------------------------

# see if any coeffs cross zero line (i.e. sig @ 5% level)
print(mod$model_output, pars = c("fu_coeff"), probs = c(0.025, 0.50, 0.975))
# only test 2, thr 2 is sig 

sens_ests <- round(summary(mod$model_output, pars = c("Se_range"), probs = c(0.025, 0.50, 0.975))$summary,2)  ; sens_ests
spec_ests <- round(summary(mod$model_output, pars = c("Sp_range"), probs = c(0.025, 0.50, 0.975))$summary,2)  ; spec_ests


sens_vec <- sens_ests[,5][c(seq(from=1, to=300), seq(from=501, to=900),seq(from=1001, to=2500) )]
sens_vec_L <- sens_ests[,4][c(seq(from=1, to=300), seq(from=501, to=900),seq(from=1001, to=2500) )]
sens_vec_U <- sens_ests[,6][c(seq(from=1, to=300), seq(from=501, to=900),seq(from=1001, to=2500) )]
spec_vec <- spec_ests[,5][c(seq(from=1, to=300), seq(from=501, to=900),seq(from=1001, to=2500) )]
spec_vec_L <- spec_ests[,4][c(seq(from=1, to=300), seq(from=501, to=900),seq(from=1001, to=2500) )]
spec_vec_U <- spec_ests[,6][c(seq(from=1, to=300), seq(from=501, to=900),seq(from=1001, to=2500) )]


# put into a tibble 
data_mr <- tibble(test_threshold_num = rep(c(1:22), each = 100))
                  
data_mr <- tibble(data_mr) ; data_mr

View(data_mr)

# make seperate test and Thresholds variables 
setDT(X2); setDT(data_mr)  
data_mr_2 <- tibble( X2[data_mr, mult = "first", on = "test_threshold_num", nomatch=0L] )
data_mr_2

# dataset with all tests and thresholds
data_mr_3 <- tibble(data_mr_2) %>% 
              rename(Threshold = Category,  Test = test) %>%
              mutate(Se = sens_vec, 
                     Se_L = sens_vec_L, 
                     Se_U = sens_vec_U, 
                     Sp = spec_vec, 
                     Sp_L = spec_vec_L, 
                     Sp_U = spec_vec_U,
                     follow_up = rep(data$cts_cov_points, 22)) %>%
              select(Test, 
                     test_num,
                     Threshold,
                     threshold_num,
                     follow_up,
                     Se, Se_L,Se_U, 
                     Sp, Sp_L,Sp_U) %>%
              mutate(Test = factor(Test), 
                     Threshold = factor(Threshold)) 

data_mr_3_Se <- data_mr_3 %>% 
  select(-Sp, 
         -Sp_L, 
         -Sp_U) %>% 
  mutate(Measure = "Sensitivity") %>% 
  rename(Median = Se, 
         Lower = Se_L, 
         Upper = Se_U)
data_mr_3_Se

data_mr_3_Sp <- data_mr_3 %>% 
  select(-Se, 
         -Se_L, 
         -Se_U) %>% 
  mutate(Measure = "Specificity") %>%
  rename(Median = Sp, 
         Lower = Sp_L, 
         Upper = Sp_U)
data_mr_3_Sp

data_mr_4 <- rbind(data.frame(data_mr_3_Se), data.frame(data_mr_3_Sp))

data_mr_5 <- tibble(data_mr_4)
data_mr_5





View(data_mr_3)
# data set containing only test/thr combos with max rank 

data_mr_ranks <- data_mr_3 %>%
              # use thresholds according to Youden index - change accordingly if using weighted youden or DOR
              filter(  test_num == 1 & threshold_num == ranks[1] |
                       test_num == 2 & threshold_num == ranks[2] |
                       test_num == 3 & threshold_num == ranks[3] |
                       test_num == 4 & threshold_num == ranks[4] |
                       test_num == 5 & threshold_num == ranks[5] )

data_mr_ranks_Se <- data_mr_ranks %>% 
                    select(-Sp, 
                           -Sp_L, 
                           -Sp_U) %>% 
                    mutate(Measure = "Sensitivity") %>% 
                    rename(Median = Se, 
                           Lower = Se_L, 
                           Upper = Se_U)
data_mr_ranks_Se

data_mr_ranks_Sp <- data_mr_ranks %>% 
                    select(-Se, 
                           -Se_L, 
                           -Se_U) %>% 
                    mutate(Measure = "Specificity") %>%
                    rename(Median = Sp, 
                           Lower = Sp_L, 
                           Upper = Sp_U)
data_mr_ranks_Sp

data_mr_ranks2 <- rbind(data.frame(data_mr_ranks_Se), data.frame(data_mr_ranks_Sp))

data_mr_ranks3 <- tibble(data_mr_ranks2)
data_mr_ranks3






levels(data_mr_ranks3$Test)

data_mr_ranks4 <- filter(data_mr_ranks3, Test == "STABLE-2007") 
data_mr_ranks4 <- filter(data_mr_ranks3, Test == "Static-99") 
data_mr_ranks4 <- filter(data_mr_ranks3, Test == "Static-99/STABLE-2007") 
data_mr_ranks4 <- filter(data_mr_ranks3, Test == "Static-99R") 
data_mr_ranks4 <- filter(data_mr_ranks3, Test == "Static-99R/STABLE-2007") 




# Plot (3) - follow up vs accuracy @ optimal thresholds ----------------------------------------
g <- ggplot(data = data_mr_ranks3,
            aes(y=Median,
                x =  follow_up, 
                colour = Test, 
                fill = Test)) + 
  geom_line(size=2) + 
  theme_bw() + 
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))  + 
  theme(text = element_text(size=15), 
        legend.position = "bottom") + 
  ylab("Accuracy") + 
  xlab("Follow up time (years)") + 
  facet_wrap(~Measure)


g


tiff("./plots/fig_3.tif",units = "in", width = 10, height=6, res=500, compression = "lzw")
g
dev.off()


# Plot (4) - follow-up vs accuracy @optimal thresholds w/ credible bands --------------
g1 <- ggplot(data = filter(data_mr_ranks3, Measure == "Sensitivity"),
            aes(y=Median,
                x =  follow_up, 
                colour = Test, 
                fill = Test)) + 
  geom_line(size=2) + 
  geom_ribbon(data = filter(data_mr_ranks3, Measure == "Sensitivity"), 
              aes(ymin = Lower, 
                  ymax = Upper, 
                  fill = Test), 
              colour = NA, 
              alpha = 0.25,
              size = 0) + 
  theme_bw() + 
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))  + 
  theme(text = element_text(size=15), 
        legend.position = "none") + 
  ylab("Accuracy") + 
  xlab("Follow up time (years)") + 
  facet_wrap(~ Test, ncol = 1)


g1

g2 <- ggplot(data = filter(data_mr_ranks3, Measure == "Specificity"),
             aes(y=Median,
                 x =  follow_up, 
                 colour = Test, 
                 fill = Test)) + 
  geom_line(size=2) + 
  geom_ribbon(data = filter(data_mr_ranks3, Measure == "Specificity"), 
              aes(ymin = Lower, 
                  ymax = Upper, 
                  fill = Test), 
              colour = NA, 
              alpha = 0.25,
              size = 0) + 
  theme_bw() + 
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))  + 
  theme(text = element_text(size=15), 
        legend.position = "right") + 
  ylab("Accuracy") + 
  xlab("Follow up time (years)") + 
  facet_wrap(~ Test, ncol = 1)


g2

g1 + g2

tiff("./plots/fig_4.tif",units = "in", width = 15, height=15, res=500, compression = "lzw")
g1 + g2
dev.off()





# plot (5) - Plot showing relationship between follow-up time and Se and Sp at all tests/thresholds --------------

data_mr_5_t1 <- filter(data_mr_5, Test == "STABLE-2007") 
data_mr_5_t2 <- filter(data_mr_5, Test == "Static-99") 
data_mr_5_t3 <- filter(data_mr_5, Test == "Static-99/STABLE-2007") 
data_mr_5_t4 <- filter(data_mr_5, Test == "Static-99R") 
data_mr_5_t5 <- filter(data_mr_5, Test == "Static-99R/STABLE-2007") 

gg_data <- list(data_mr_5_t1,
                data_mr_5_t2,
                data_mr_5_t3,
                data_mr_5_t4,
                data_mr_5_t5)

gg_se <- list()

for (i in 1:5) { 
gg_se[[i]] <- ggplot(data = filter(gg_data[[i]], Measure == "Sensitivity"),
             aes(y=Median,
                 x =  follow_up, 
                 colour = Test, 
                 fill = Test)) + 
            geom_line(size=2) + 
            geom_ribbon(data = filter(gg_data[[i]], Measure == "Sensitivity"), 
                        aes(ymin = Lower, 
                            ymax = Upper, 
                            fill = Test), 
                        colour = NA, 
                        alpha = 0.25,
                        size = 0) + 
            theme_bw() + 
            scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))  + 
            theme(text = element_text(size=15), 
                  legend.position = "none") + 
            ylab("Accuracy") + 
            xlab("Follow up time (years)") + 
            ggtitle(paste0(as.character(gg_data[[i]]$Test[1]))) + 
            facet_wrap(~ Test + Threshold,dir = "h", ncol = length(unique(gg_data[[i]]$threshold_num)))

}

gg_se[[1]] + gg_se[[2]] + gg_se[[3]] + gg_se[[4]] + gg_se[[5]]+ 
  plot_layout(nrow =5)


tiff("./plots/fig_5_se.tif",units = "in", width = 15, height=15, res=500, compression = "lzw")
gg_se[[1]] + gg_se[[2]] + gg_se[[3]] + gg_se[[4]] + gg_se[[5]]+ 
  plot_layout(nrow =5)
dev.off()


gg_sp <- list()

for (i in 1:5) { 
  gg_sp[[i]] <- ggplot(data = filter(gg_data[[i]], Measure == "Specificity"),
                       aes(y=Median,
                           x =  follow_up, 
                           colour = Test, 
                           fill = Test)) + 
    geom_line(size=2) + 
    geom_ribbon(data = filter(gg_data[[i]], Measure == "Specificity"), 
                aes(ymin = Lower, 
                    ymax = Upper, 
                    fill = Test), 
                colour = NA, 
                alpha = 0.25,
                size = 0) + 
    theme_bw() + 
    scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))  + 
    theme(text = element_text(size=15), 
          legend.position = "none") + 
    ylab("Accuracy") + 
    xlab("Follow up time (years)") + 
    ggtitle(paste0(as.character(gg_data[[i]]$Test[1]))) + 
    facet_wrap(~ Test + Threshold,dir = "h", ncol = length(unique(gg_data[[i]]$threshold_num)))
  
}

gg_sp[[1]] + gg_sp[[2]] + gg_sp[[3]] + gg_sp[[4]] + gg_sp[[5]] + 
  plot_layout(nrow =5)


tiff("./plots/fig_5_sp.tif",units = "in", width = 15, height=15, res=500, compression = "lzw")
gg_sp[[1]] + gg_sp[[2]] + gg_sp[[3]] + gg_sp[[4]] + gg_sp[[5]] + 
  plot_layout(nrow =5)
dev.off()
























