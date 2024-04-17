
### These are the three working models, plus a data generating function

library(rstan)
library(ggplot2)

generate_yearly_simp = function(
    K,   # Number of years
    J,   # Number of modules
    Nm,  # Number of men module places (same each year)
    Nw,  # Number of women module places (same each year)
    alpha,   # J vector alpha  
    beta,    # J vector of beta
    sig_alpha, # J vector of variances for alpha
    sig_beta, # J vector of variances for beta
    seed = 2863737    # seed for reproducibility
){
  
  nm_out = matrix(NA, nrow=K, ncol=J)
  nw_out = matrix(NA, nrow=K, ncol=J)
  pjm_out = matrix(NA, nrow=K, ncol=J)
  pjw_out = matrix(NA, nrow=K, ncol=J)
  
  
  ### For a particular year I need alpha and beta vectors
  for (i in 1:K){
    alpha_yearly = sapply(1:J, function(j){rnorm(1, mean=alpha[j], sd = sig_alpha[j])})
    beta_yearly = sapply(1:J, function(j){rnorm(1, mean=beta[j], sd = sig_beta[j])})
    denom_m = sum(sapply(1:J, function(j){exp(alpha_yearly[j])}))
    denom_w = sum(sapply(1:J, function(j){exp(alpha_yearly[j] + beta_yearly[j])}))
    pjm = sapply(1:J, function(j){exp(alpha_yearly[j])/denom_m})
    pjw = sapply(1:J, function(j){exp(alpha_yearly[j]+beta_yearly[j])/denom_w})
    
    Mj = rmultinom(n=1, size=Nm, prob = pjm)
    Wj = rmultinom(n=1, size=Nw, prob = pjw)
    
    nm_out[i,] = Mj
    nw_out[i,] = Wj
    pjm_out[i,] = pjm
    pjw_out[i,] = pjw
  }
  
  list(nm = nm_out, nw = nw_out, pjm = pjm_out, pjw = pjw_out
  )
}


## Model for aggregated data (not yearly)

cs_model_agg = 
  "
  data {
    int<lower=1> J;             // Number of modules
    int<lower=0> Mj[J];         // Total men students on module
    int<lower=0> Wj[J];         // Total of women on module
    real<lower=0> a_sigalp[J];   // sd for sig_alpha
    real<lower=0> b_sigalp[J];
    real<lower=0> a_sigbet[J];   // sd for sig_beta
    real<lower=0> b_sigbet[J];
  }
  

  parameters {
    vector[J] alpha; // vector of intercept parameters
    vector[J] beta;   // vector of women coefficient/parameters
    vector<lower=0>[J] sig_alpha;
    vector<lower=0>[J] sig_beta;
  }
  
  transformed parameters {
    vector[J] pjm;
    vector[J] pjw;
    vector[J] denom_m_vec;
    vector[J] denom_w_vec;
    real denom_m;
    real denom_w;
    
    
    for (i in 1:J){
      denom_m_vec[i] = exp(alpha[i]);
      denom_w_vec[i] = exp(alpha[i] + beta[i]);
    }
     
    denom_m = sum(denom_m_vec);
    denom_w = sum(denom_w_vec);
    
    for (j in 1:J){
      pjm[j] = exp(alpha[j])/denom_m;
      pjw[j] = exp(alpha[j] + beta[j])/denom_w;
    }
  }
  
  model {
    target += inv_gamma_lpdf(sig_alpha|a_sigalp,b_sigalp);
    target += inv_gamma_lpdf(sig_beta|a_sigbet,b_sigbet);


    alpha ~ normal(0, sig_alpha);         // Priors for each module
//    beta ~ normal(beta0, sig_beta);
    beta ~ normal(0, sig_beta);
    Mj ~ multinomial(pjm);          // Number of men students per module
    Wj ~ multinomial(pjw);
  }
"

## Yearly data
## With no constraint on beta

cs_yearly = 
  "
  data {
    int<lower=1> J;             // Number of modules
    int<lower=1> K;             // Number of years
    int<lower=0> Mj[K,J];         // Total men students on J modules over K years (years as rows)
    int<lower=0> Wj[K,J];         // Total of women on J modules over K years
    real<lower=0> a_sigalp[J];   // sd for sig_alpha
    real<lower=0> b_sigalp[J];
    real<lower=0> a_sigbet[J];   // sd for sig_beta
    real<lower=0> b_sigbet[J];
  }
  

  parameters {
    real alpha[J]; // vector of intercept parameters
    real beta[J];   // vector of women coefficient/parameters
    real<lower=0> sig_alpha[J];
    real<lower=0> sig_beta[J];
  }
  
  transformed parameters {
    simplex[J] pjm;
    simplex[J] pjw;
    vector[J] denom_m_vec;
    vector[J] denom_w_vec;
    real denom_m;
    real denom_w;
    
    for (i in 1:J){
      denom_m_vec[i] = exp(alpha[i]);
      denom_w_vec[i] = exp(alpha[i] + beta[i]);
    }
     
    denom_m = sum(denom_m_vec);
    denom_w = sum(denom_w_vec);
    
    for (j in 1:J){
      pjm[j] = exp(alpha[j])/denom_m;
      pjw[j] = exp(alpha[j] + beta[j])/denom_w;
    }
  }

  model {
    

    target += inv_gamma_lpdf(sig_alpha|a_sigalp,b_sigalp);
    target += inv_gamma_lpdf(sig_beta|a_sigbet,b_sigbet);

//    sig_beta ~ gamma(2,2);
  
    alpha ~ normal(0, sig_alpha);         // Priors for each module
    beta ~ normal(0, sig_beta);

    for (k in 1:K){
      Mj[k] ~ multinomial(pjm);          // Number of men students per module
      Wj[k] ~ multinomial(pjw);
    }
  }
"


## Model for yearly data
## More conservative version with beta0 term

cs_yearly_beta0 = 
  "
  data {
    int<lower=1> J;             // Number of modules
    int<lower=1> K;             // Number of years
    int<lower=0> Mj[K,J];         // Total men students on J modules over K years (years as rows)
    int<lower=0> Wj[K,J];         // Total of women on J modules over K years
    real<lower=0> a_sigalp[J];   // sd for sig_alpha
    real<lower=0> b_sigalp[J];
    real<lower=0> a_sigbet[J];   // sd for sig_beta
    real<lower=0> b_sigbet[J]; 
    real<lower=0> sig_beta0;
  }
  

  parameters {
    real alpha[J]; // vector of intercept parameters
    real beta[J];   // vector of women coefficient/parameters
    real beta0;
    real<lower=0> sig_alpha[J];
    real<lower=0> sig_beta[J];
  }
  
  transformed parameters {
    simplex[J] pjm;
    simplex[J] pjw;
    vector[J] denom_m_vec;
    vector[J] denom_w_vec;
    real denom_m;
    real denom_w;
    
    for (i in 1:J){
      denom_m_vec[i] = exp(alpha[i]);
      denom_w_vec[i] = exp(alpha[i] + beta[i]);
    }
     
    denom_m = sum(denom_m_vec);
    denom_w = sum(denom_w_vec);
    
    for (j in 1:J){
      pjm[j] = exp(alpha[j])/denom_m;
      pjw[j] = exp(alpha[j] + beta[j])/denom_w;
    }
  }

  model {
    

    target += inv_gamma_lpdf(sig_alpha|a_sigalp,b_sigalp);
    target += inv_gamma_lpdf(sig_beta|a_sigbet,b_sigbet);

//    sig_beta ~ gamma(2,2);
  
    alpha ~ normal(0, sig_alpha);         // Priors for each module
    beta ~ normal(beta0, sig_beta);
    beta0 ~ normal(0, sig_beta0);
    
    for (k in 1:K){
      Mj[k] ~ multinomial(pjm);          // Number of men students per module
      Wj[k] ~ multinomial(pjw);
    }
  }
"

### Testing

## compare this to what the aggregate data would produce from the aggregate model

gen10s = generate_yearly_simp(K=10, J=6, Nm=200, Nw=50, 
                              alpha = c(-2, -1, 0, 0, 1, 2), beta = c(2, 0, -1, 0, -2, 1),
                              sig_alpha = c(0.5, 0.01, 1, 0.5, 0.01, 1), sig_beta= c(0.01, 1, 0.5, 0.5, 0.01, 1))

gen100s = generate_yearly_simp(K=100, J=6, Nm=200, Nw=50, 
                               alpha = c(-2, -1, 0, 0, 1, 2), beta = c(2, 0, -1, 0, -2, 1),
                               sig_alpha = c(0.5, 0.01, 1, 0.5, 0.01, 1), sig_beta= c(0.01, 1, 0.5, 0.5, 0.01, 1))


gen10s_yearly= stan(
  model_code= cs_yearly, 
  data=list(K=10, J=6, Mj = gen10s$nm, Wj = gen10s$nw, 
            a_sigalp = rep(3,6), b_sigalp = rep(1.5,6),
            a_sigbet = rep(3,6), b_sigbet = rep(1.5,6)), 
  iter=5000)
stan_plot(gen10s_yearly, pars = c("alpha", "beta"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest
stan_plot(gen10s_yearly, pars = c("pjw", "pjm"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest
print(gen10s_yearly, probs=c(0.1, 0.9))

gen10s_agg = stan(
  model_code= cs_model_agg, 
  data=list(J=6, Mj = colSums(gen10s$nm), Wj = colSums(gen10s$nw), 
            a_sigalp = rep(3,6), b_sigalp = rep(1.5,6),
            a_sigbet = rep(3,6), b_sigbet = rep(1.5,6)), 
  iter=5000)
stan_plot(gen10s_agg, pars = c("alpha", "beta"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest


gen100s_yearly= stan(
  model_code= cs_yearly, 
  data=list(K=100, J=6, Mj = gen100s$nm, Wj = gen100s$nw, 
            a_sigalp = rep(3,6), b_sigalp = rep(1.5,6),
            a_sigbet = rep(3,6), b_sigbet = rep(1.5,6)), 
  iter=5000)
stan_plot(gen100s_yearly, pars = c("alpha", "beta"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest
stan_plot(gen100s_yearly, pars = c("sig_alpha", "sig_beta"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest
print(gen100s_yearly, probs=c(0.1, 0.9))

## There's no difference here, which makes me think things are going wrong!

gen100s_yearly= stan(
  model_code= cs_yearly, 
  data=list(K=100, J=6, Mj = gen100s$nm, Wj = gen100s$nw, 
            a_sigalp = rep(3,6), b_sigalp = rep(1.5,6),
            a_sigbet = rep(3,6), b_sigbet = rep(1.5,6)), 
  iter=5000)
stan_plot(gen100s_yearly, pars = c("alpha", "beta"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest

gen100s_agg = stan(
  model_code= cs_model_agg, 
  data=list(J=6, Mj = colSums(gen100s$nm), Wj = colSums(gen100s$nw), 
            a_sigalp = rep(3,6), b_sigalp = rep(1.5,6),
            a_sigbet = rep(3,6), b_sigbet = rep(1.5,6)), 
  iter=5000)
stan_plot(gen100s_agg, pars = c("alpha", "beta"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest



## Actual data

cs_df = data.frame(
  modID = 0:27,
  n_women = c(12, 17, 15, 16, 1, 21, 17, 24, 11, 11, 13, 1, 
              23, 11, 6, 5, 4, 3, 4, 5, 6, 0, 5, 5, 11, 7, 1, 3),
  n_men = c(85, 118, 94, 125, 40, 74, 102, 151, 80, 111, 90, 9, 122, 80,
            8, 21, 32, 18, 33, 40, 35, 6, 45, 29, 53, 47, 16, 23)
)

cs_agg= stan(
  model_code= cs_model_agg, 
  data=list(J=28, Mj = cs_df$n_men, Wj = cs_df$n_women, 
            a_sigalp = rep(3,28), b_sigalp = rep(1.5,28),
            a_sigbet = rep(3,28), b_sigbet = rep(1.5,28)), 
  iter=5000)
stan_plot(cs_agg, pars = c("alpha", "beta"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest
stan_plot(cs_agg, pars = c("pjw", "pjm"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest


## A slightly more complicated yearly model, which hopefully allows the variance to
## reflect the actual variance

cs_logit_years = 
  "
  data {
    int<lower=1> J;             // Number of modules
    int<lower=1> K;             // Number of years
    int<lower=0> Mj[K,J];         // Total men students on module over K years (years as rows)
    int<lower=0> Wj[K,J];         // Total of women on module over K years
    real<lower=0> a_sigalp[J];   // sd for sig_alpha
    real<lower=0> b_sigalp[J];
    real<lower=0> a_sigbet[J];   // sd for sig_beta
    real<lower=0> b_sigbet[J]; 
  }
  

  parameters {
    real alpha[J]; // vector of intercept parameters
    real alpha_jk[K,J];  // intercept parameters for each module and year
    real beta[J];   // vector of women coefficient/parameters
    real beta_jk[K,J]; // women coefficients for each module and year
    real<lower=0> sig_alpha[J];
    real<lower=0> sig_beta[J];
    real<lower=0> sig_alpjk[J];
    real<lower=0> sig_betjk[J];
  }
  
  transformed parameters {
    matrix[K,J] pjm;
    matrix[K,J] pjw;
    matrix[K,J] denom_m_vec;
    matrix[K,J] denom_w_vec;
    vector[K] denom_m;
    vector[K] denom_w;
    
    for (k in 1:K){
      for (j in 1:J){
        denom_m_vec[k,j] = exp(alpha_jk[k,j]);
        denom_w_vec[k,j] = exp(alpha_jk[k,j] + beta_jk[k,j]);
      }
      denom_m[k] = sum(denom_m_vec[k]);
      denom_w[k] = sum(denom_w_vec[k]);
    }

    for (k in 1:K){
      for (j in 1:J){
        pjm[k,j] = exp(alpha_jk[k,j])/denom_m[k];
        pjw[k,j] = exp(alpha_jk[k,j] + beta_jk[k,j])/denom_w[k];
      }
    }
    
  }
  model {
// The sig_alpha and sig_alpjk (ditto sig_beta, sig_betjk) have the same prior, but this seems reasonable
// as a first try
    target += inv_gamma_lpdf(sig_alpha|a_sigalp,b_sigalp);
    target += inv_gamma_lpdf(sig_beta|a_sigbet,b_sigbet);

    target += inv_gamma_lpdf(sig_alpjk|a_sigalp,b_sigalp);
    target += inv_gamma_lpdf(sig_betjk|a_sigbet,b_sigbet);


    alpha ~ normal(0, sig_alpha);         // Priors for each module
    beta ~ normal(0, sig_beta);
    
    for (k in 1:K){
      alpha_jk[k] ~ normal(alpha, sig_alpjk);  // Priors for each module and year
      beta_jk[k] ~ normal(beta, sig_betjk);
      Mj[k] ~ multinomial(to_vector(pjm[k]));          // Number of men students per module
      Wj[k] ~ multinomial(to_vector(pjw[k]));
    }
  }
"

gen10s_yearly_logit = stan(
  model_code= cs_logit_years, 
  data=list(K=10, J=6, Mj = gen10s$nm, Wj = gen10s$nw, 
            a_sigalp = rep(3,6), b_sigalp = rep(1.5,6),
            a_sigbet = rep(3,6), b_sigbet = rep(1.5,6)), 
  iter=5000)

stan_plot(gen10s_yearly_logit, pars = c("alpha", "beta"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest

gen100s_yearly_logit = stan(
  model_code= cs_logit_years, 
  data=list(K=100, J=6, Mj = gen100s$nm, Wj = gen100s$nw, 
            a_sigalp = rep(3,6), b_sigalp = rep(1.5,6),
            a_sigbet = rep(3,6), b_sigbet = rep(1.5,6)), 
  iter=5000)

stan_plot(gen100s_yearly_logit, pars = c("alpha", "beta"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest
stan_plot(gen100s_yearly_logit, pars = c("sig_alpha", "sig_beta"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest


## Slightly less outlandish variances 

gen100s2 = generate_yearly_simp(K=100, J=6, Nm=200, Nw=50, 
                               alpha = c(-2, -1, 0, 0, 1, 2), beta = c(2, 0, -1, 0, -2, 1),
                               sig_alpha = c(0.5, 0.25, 1, 0.5, 0.25, 1), sig_beta= c(0.25, 1, 0.5, 0.5, 0.25, 1))


gen100s2_yearly= stan(
  model_code= cs_yearly, 
  data=list(K=100, J=6, Mj = gen100s2$nm, Wj = gen100s2$nw, 
            a_sigalp = rep(3,6), b_sigalp = rep(1.5,6),
            a_sigbet = rep(3,6), b_sigbet = rep(1.5,6)), 
  iter=10000)

stan_plot(gen100s2_yearly, pars = c("alpha", "beta"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest
stan_plot(gen100s2_yearly, pars = c("sig_alpha", "sig_beta"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest



gen100s2_yearly_logit = stan(
  model_code= cs_logit_years, 
  data=list(K=100, J=6, Mj = gen100s2$nm, Wj = gen100s2$nw, 
            a_sigalp = rep(3,6), b_sigalp = rep(1.5,6),
            a_sigbet = rep(3,6), b_sigbet = rep(1.5,6)), 
  iter=10000)

stan_plot(gen100s2_yearly_logit, pars = c("alpha", "beta"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest
stan_plot(gen100s2_yearly_logit, pars = c("sig_alpha", "sig_beta"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest
