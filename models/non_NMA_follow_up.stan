// Model is from Nyaga et al - arm-based NMA model, "traditional" (non-copula) normal-binomial

// Code based on that provided from Nyaga et al. supp. material - with the some small modifications:
// using more efficient non-centered param for between-study model, uses LKJ corr priors, and calculation
// of parameters of interest in the generated quantities block, allows calculation of K-fold CV which
// is more reliable than WAIC and can be used when LOO fails (i.e. many pareto-k's > 0.7)
 
functions {  // make function to count # of unique elements in a vector - equiv to unique() in R. 
            int unique(int[] x) {
                int res[num_elements(x)] = sort_desc(x);   // Sort the elements of x in descending orde
                int count = 1;
                for (i in 2:num_elements(res))  {
                   if ( res[i] != res[i-1])  {
                     count = count + 1;
                    } else { 
                      count = count;
                    }
                }
                return(count); 
            }
}
                

data {
              int N; // number of obs. / rows 
              int n_tests; // total number of tests 
              int n_thresholds[n_tests]; // number of categories (i.e. # of cut-offs + 1) for each ordinal-scale test
              int n_test_thresholds;
              int test_thr_param;
              int fu_summary;
              int TP[N];
              int FN[N];
              int FP[N];
              int TN[N];
              int Study[N];
              int Test[N];
              int Threshold[N];
              int Test_Threshold[N];
              vector[N] fu;
              real prior_intercept_thr_1_sd;
              real prior_intercept_thr_1_diseased_mean;
              real prior_intercept_thr_1_non_diseased_mean;
              real prior_intercept_thr_2_sd;
              real prior_intercept_thr_2_diseased_mean;
              real prior_intercept_thr_2_non_diseased_mean;
              real prior_intercept_other_thr_sd;
              real prior_coeff_sd;
              real prior_sigma_sd;
              vector[100] cts_cov_points;
              vector<lower=0,upper=1>[n_tests] prev;
              vector<lower=-1,upper=1>[1] youden_weight;
              int<lower=0,upper=1> diff_SDs_corr; //  indicator for whether to share SDs and Corr between tests (unstructured) or not (homogenous)
              int<lower=0,upper=1> holdout[N]; //index whether the observation should be used (for K-fold CV for model comparison)
}

transformed data { 
              int max_Threshold;     // Stan does not support staggered arrays, so make a variable make a variable "max_Threshold"  and for each test t, when e.g. looping, only index the categories up to category[t], not max_Threshold
              int n_studies; // total number of studies 
              int pos[N];
              int neg[N];
              
                for (n in 1:N) { 
                  pos[n] = TP[n] + FN[n];
                  neg[n] = TN[n] + FP[n];
                }
              
              max_Threshold = max(n_thresholds);
              n_studies = unique(Study);
}

parameters {
              matrix<lower=-20, upper = 20>[max_Threshold, 2] mu[n_tests];  
              cholesky_factor_corr[2] L_Omega[n_tests]; // prior for cholesky corr
              matrix<lower=0>[n_tests, 2] sigma; 
              matrix[N, 2] z;
              matrix[max_Threshold, 2] fu_coeff[n_tests];  // coefficient for cts. covariate for follow-up time 
}

transformed parameters {
              vector[2] logit_pi[N];
              vector[2] log_lik[N];

                           
  for (n in 1:N) {

       if (diff_SDs_corr == 1)    logit_pi[n, ] =  to_vector(mu[Test[n], Threshold[n], ]) + // mu is intercept for test t at threshold t(c) (since have coeff for FU)
                                                   to_vector(fu_coeff[Test[n], Threshold[n], ]*fu[n] )  +
                                                   diag_pre_multiply(sigma[Test[n], ], L_Omega[Test[n],,]) * to_vector(z[n, ]);
       else                       logit_pi[n, ] =  to_vector(mu[Test[n], Threshold[n], ]) + 
                                                   to_vector(fu_coeff[Test[n], Threshold[n], ]*fu[n] ) + 
                                                   diag_pre_multiply(sigma[1, ], L_Omega[1,,]) *  to_vector(z[n, ]) ;
       
              
                           // Likelihood Model -  Pointwise (i.e. observation level) log-Likelihood 
                        log_lik[n, 1]  =  binomial_logit_lpmf(TP[n]  | pos[n] , logit_pi[n, 1] );
                        log_lik[n, 2]  =  binomial_logit_lpmf(TN[n]  | neg[n] , logit_pi[n, 2] );
      }
}

model {

           //Prior Model
    for (n in 1:N) {
       if ( Threshold[n] == 1) { 
         mu[Test[n], Threshold[n], 1] ~ normal(prior_intercept_thr_1_diseased_mean, prior_intercept_thr_1_sd); 
         mu[Test[n], Threshold[n], 2] ~ normal(prior_intercept_thr_1_non_diseased_mean, prior_intercept_thr_1_sd); 
       }    
      else if ( Threshold[n] == 2) { 
         mu[Test[n], Threshold[n], 1] ~ normal(prior_intercept_thr_2_diseased_mean, prior_intercept_thr_2_sd); 
         mu[Test[n], Threshold[n], 2] ~ normal(prior_intercept_thr_2_non_diseased_mean, prior_intercept_thr_2_sd);
      }
       else {                
         mu[Test[n], Threshold[n], 1] ~ normal(0, prior_intercept_other_thr_sd); 
         mu[Test[n], Threshold[n], 2] ~ normal(0, prior_intercept_other_thr_sd);
         }                                                                                           
    }
    

 
       to_vector(sigma) ~ normal(0, prior_sigma_sd);
       to_vector(z) ~ std_normal();
       
     for (t in 1:n_tests) 
      to_vector(fu_coeff[t,,]) ~ normal(0, prior_coeff_sd); 
       
    for (t in 1:n_tests) 
       L_Omega[t,,] ~ lkj_corr_cholesky(2);

  // Likelihood Model
    for (n in 1:N) {
      if(holdout[n] == 0) {
        target += log_lik[n, 1];
        target += log_lik[n, 2];
       }
    }
    
}

generated quantities { 
                 matrix[n_tests, max_Threshold] Se; 
                 matrix[n_tests, max_Threshold] Sp;  
                 matrix[n_tests, max_Threshold] lSe; 
                 matrix[n_tests, max_Threshold] lSp; 
                 matrix[max_Threshold, 2] pred[n_tests]; 
                 matrix[n_tests, max_Threshold] Se_pred; 
                 matrix[n_tests, max_Threshold] Sp_pred;  
                 matrix[n_tests, max_Threshold] lSe_pred; 
                 matrix[n_tests, max_Threshold] lSp_pred; 
                 corr_matrix[2] Omega[n_tests]; 
                 matrix[2,2] Sigma[n_tests]; 
                 matrix[n_tests, max_Threshold] Fp; 
                 matrix[n_tests, max_Threshold] LRpos;  
                 matrix[n_tests, max_Threshold] LRneg; 
                 matrix[n_tests, max_Threshold] DOR; 
                 matrix[n_tests, max_Threshold] PPV; 
                 matrix[n_tests, max_Threshold] NPV; 
                 matrix[n_tests, max_Threshold] Youden_index; 
                 matrix[n_tests, max_Threshold] Youden_index_weighted; 
                 matrix[n_tests, 2] tausq;
                 vector[2] sigmabsq[n_tests];
                 matrix[n_tests, n_tests] sigmasq[2];
                 matrix[n_tests, n_tests] rho[2];
                 vector[n_tests] rho12;
                 matrix[2,2] Sigma_bs; 
                 matrix[n_tests, 2] tau;
                 matrix[max_Threshold,100] Se_cov[n_tests];
                 matrix[max_Threshold,100] Sp_cov[n_tests];

  
        // Between-study variance/cor parameters  
      if (diff_SDs_corr == 1) { 
        for (t in 1:n_tests) { 
                  Omega[t,] = multiply_lower_tri_self_transpose(L_Omega[t,]); // correlation matrix 
                  Sigma[t,] = quad_form_diag(Omega[t,], sigma[t,]);               // var-cov matrix 
        }
      }
      else { 
                for (t in 1:n_tests) { 
                  Omega[t,] = multiply_lower_tri_self_transpose(L_Omega[1,]); // correlation matrix 
                  Sigma[t,] = quad_form_diag(Omega[1,], sigma[1,]);               // var-cov matrix
        }
      }

         for (t in 1:n_tests) { 
           for (c in 1:n_thresholds[t]) { 
                Se[t,c] = inv_logit(mu[t,c,1] + fu_summary*fu_coeff[t,c,1]) ; 
                Sp[t,c] = inv_logit(mu[t,c,2] + fu_summary*fu_coeff[t,c,2]); 
                lSe[t,c] = mu[t,c,1];
                lSp[t,c] = mu[t,c,2];
                pred[t,c,] = to_row_vector(multi_normal_rng( to_row_vector(mu[t,c,] + fu_coeff[t,c,] )  , Sigma[t,] )); 
                lSe_pred[t,c] = pred[t,c,1];
                lSp_pred[t,c] = pred[t,c,2];
                Se_pred[t,c] = inv_logit(pred[t,c,1]);
                Sp_pred[t,c] = inv_logit(pred[t,c,2]);
          }
         }
          
            // summary Se and Sp for meta-regression (within range of min and max observed point)
      for (t in 1:n_tests) { 
        for (c in 1:n_thresholds[t]) { 
           for (i in 1:100) { 
                    Se_cov[t,c,i] =  inv_logit(mu[t,c,1] + fu_coeff[t,c,1]*cts_cov_points[i]);
                    Sp_cov[t,c,i] =  inv_logit(mu[t,c,2] + fu_coeff[t,c,2]*cts_cov_points[i]);
           }
        }
      }
                
                


        // Other summary non-comparative estimates 
                 for (t in 1:n_tests) { 
                     for (c in 1:n_thresholds[t]) { 
                          Fp[t,c] = 1 - Sp[t,c];
                          
                        // LR's
                         	LRpos[t,c] = Se[t,c]/(1-Sp[t,c]);
          	              LRneg[t,c] = (1-Se[t,c])/Sp[t,c];
          	              
          	            // PPV and NPV (using inputted prev data vector)
          	            PPV[t,c] = Se[t,c]*prev[t]/( Se[t,c]*prev[t] + (1-Sp[t,c])*(1-prev[t]) );
          	            NPV[t,c] = Sp[t,c]*prev[t]/( (1-Se[t,c])*prev[t] + Sp[t,c]*(1-prev[t]) );
                       }
	              
                     }
                 
        // Rank statistics
                 for (t in 1:n_tests) { 
                     for (c in 1:n_thresholds[t]) { 
                       
                        // Youden index
                        Youden_index[t,c] =  Se[t,c] + Sp[t,c] - 1;
                        
                        // Weighted Youden index (weight = 0.5 -> same as youden, weight > 0.5 -> more importance on Se)
                        Youden_index_weighted[t,c] = 2 * ( youden_weight[1]*Se[t,c] + (1 - youden_weight[1])*Sp[t,c] ) ;
                        
                        // DOR
                      	DOR[t,c] = (Se[t,c]*Sp[t,c])/((1-Se[t,c])*(1-Sp[t,c]));
    
                     }
                 }
                 
                 

}



