## CS women model
## Symmetric version, where the men are alpha - beta/2 
## and the women are alpha + beta/2

## Yearly data
## With no constraint on beta
## We do need to allow the yearly varying of the alpha and beta
## because otherwise it is the same as the aggregate data
## As shown (sort of) in cs_women_binomial.R

library(rstan)
library(ggplot2)

generate_yearly = function(
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




cs_sym = 
  "
  data {
    int<lower=1> J;             // Number of modules
    int<lower=1> K;             // Number of years
    int<lower=0> Mj[K,J];         // Total men on module over K years (years as rows)
    int<lower=0> Wj[K,J];         // Total women on module over K years (years as rows)
    real<lower=0> a_sigalp;
    real<lower=0> b_sigalp;
    real<lower=0> a_sigbet;
    real<lower=0> b_sigbet;
    real<lower=0> a_sigalpjk;   // sd for sig_alpjk (same for all J?)
    real<lower=0> b_sigalpjk; 
    real<lower=0> a_sigbetjk;   // sd for sig_alpjk (same for all J?)
    real<lower=0> b_sigbetjk; 
  }
  
  parameters {
    real alpha[J]; // vector of intercept parameters
    real alpha_jk[K,J];  // intercept parameters for each module and year
    real<lower=0> sig_alpjk[J];
    real<lower=0> sig_alpha;
    real beta[J]; // vector of parameters
    real beta_jk[K,J];  // parameters for each module and year
    real<lower=0> sig_betjk[J];
    real<lower=0> sig_beta;
  }
  
  transformed parameters {
    matrix[K,J] pjm;
    matrix[K,J] pjw;
    matrix[K,J] denom_m_vec;
    vector[K] denom_m;
    matrix[K,J] denom_w_vec;
    vector[K] denom_w;

    for (k in 1:K){
      for (j in 1:J){
        denom_m_vec[k,j] = exp(alpha_jk[k,j] - 0.5*beta_jk[k,j]);
        denom_w_vec[k,j] = exp(alpha_jk[k,j] + 0.5*beta_jk[k,j]);
      }
      denom_m[k] = sum(denom_m_vec[k]);
      denom_w[k] = sum(denom_w_vec[k]);
    }

    for (k in 1:K){
      for (j in 1:J){
        pjm[k,j] = exp(alpha_jk[k,j] - 0.5*beta_jk[k,j])/denom_m[k];
        pjw[k,j] = exp(alpha_jk[k,j] + 0.5*beta_jk[k,j])/denom_w[k];
      }
    }
    
  }
  model {
    target += inv_gamma_lpdf(sig_alpha|a_sigalp,b_sigalp);
    target += inv_gamma_lpdf(sig_alpjk|a_sigalpjk,b_sigalpjk);
    target += inv_gamma_lpdf(sig_beta|a_sigbet,b_sigbet);
    target += inv_gamma_lpdf(sig_betjk|a_sigbetjk,b_sigbetjk);

    alpha ~ normal(0, sig_alpha);         // Priors for each module
    beta ~ normal(0, sig_beta);         // Priors for each module

    for (k in 1:K){
      alpha_jk[k] ~ normal(alpha, sig_alpjk);  // Priors for each module and year
      beta_jk[k] ~ normal(beta, sig_betjk);  // Priors for each module and year
      Mj[k] ~ multinomial(to_vector(pjm[k]));          // Number of men students per module
      Wj[k] ~ multinomial(to_vector(pjw[k]));          // Number of men students per module
    }
  }
"


gen20_sig0_same = generate_yearly(
  K=20,  J=5,  Nm=100, Nw = 50,  
  alpha = rep(0,5), beta = rep(0,5),  
  sig_alpha = rep(0,5), sig_beta = rep(0,5),
  seed = 2863737    # seed for reproducibility
)

gen20_sig03_same = generate_yearly(
  K=20,  J=5,  Nm=100, Nw = 50,  
  alpha = rep(0,5), beta = rep(0,5),  
  sig_alpha = rep(0.3,5), sig_beta = rep(0.3,5),
  seed = 2863737    # seed for reproducibility
)

gen20_sig05_same = generate_yearly(
  K=20,  J=5,  Nm=100, Nw = 50,  
  alpha = rep(0,5), beta = rep(0,5),  
  sig_alpha = rep(0.5,5), sig_beta = rep(0.5,5),
  seed = 2863737    # seed for reproducibility
)


cs20_sig0_same= stan(
  model_code= cs_sym, 
  data=list(K=20, J=5, Mj = gen20_sig0_same$nm, Wj = gen20_sig0_same$nw, 
            a_sigalp=3, b_sigalp = 0.5, a_sigbet=3, b_sigbet = 0.5, 
            a_sigalpjk = 3, b_sigalpjk = 0.5,
            a_sigbetjk = 3, b_sigbetjk = 0.5,
  iter=80000))

cs20_sig03_same= stan(
  model_code= cs_sym, 
  data=list(K=20, J=5, Mj = gen20_sig03_same$nm, Wj = gen20_sig03_same$nw, 
            a_sigalp=3, b_sigalp = 0.5, a_sigbet=3, b_sigbet = 0.5, 
            a_sigalpjk = 3, b_sigalpjk = 0.5,
            a_sigbetjk = 3, b_sigbetjk = 0.5,
            iter=80000))

cs20_sig05_same= stan(
  model_code= cs_sym, 
  data=list(K=20, J=5, Mj = gen20_sig05_same$nm, Wj = gen20_sig05_same$nw, 
            a_sigalp=3, b_sigalp = 0.5, a_sigbet=3, b_sigbet = 0.5, 
            a_sigalpjk = 3, b_sigalpjk = 0.5,
            a_sigbetjk = 3, b_sigbetjk = 0.5,
            iter=80000))


stan_plot(cs20_sig0_same, pars = c("alpha", "beta"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest
stan_plot(cs20_sig03_same, pars = c("alpha", "beta"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest


print(cs20_sig0_same, probs=c(0.1, 0.9))
print(cs20_sig03_same, probs=c(0.1, 0.9))

## Start smaller

gen10_sig0_same = generate_yearly(
  K=10,  J=2,  Nm=100, Nw = 50,  
  alpha = rep(0,2), beta = rep(0,2),  
  sig_alpha = rep(0,2), sig_beta = rep(0,2),
  seed = 2863737    # seed for reproducibility
)

gen10_sig03_same = generate_yearly(
  K=10,  J=2,  Nm=100, Nw = 50,  
  alpha = rep(0,2), beta = rep(0,2),  
  sig_alpha = rep(0.3,2), sig_beta = rep(0.3,2),
  seed = 2863737    # seed for reproducibility
)
cs10_sig0_same= stan(
  model_code= cs_yearly_sym, 
  data=list(K=10, J=2, Mj = gen10_sig0_same$nm, Wj = gen10_sig0_same$nw, 
            sig_alpha=0.25, sig_beta=0.25, 
            a_sigalpjk = 3, b_sigalpjk = 0.5,
            a_sigbetjk = 3, b_sigbetjk = 0.5,
            iter=40000))

cs10_sig03_same= stan(
  model_code= cs_yearly_sym, 
  data=list(K=10, J=2, Mj = gen10_sig03_same$nm, Wj = gen10_sig03_same$nw, 
            sig_alpha=0.25, sig_beta=0.25, 
            a_sigalpjk = 3, b_sigalpjk = 0.5,
            a_sigbetjk = 3, b_sigbetjk = 0.5,
            iter=40000))



print(cs10_sig0_same, pars = c("alpha", "beta", "alpha_jk", "beta_jk", "sig_alpjk", "sig_betjk"), probs=c(0.1, 0.9))
print(cs10_sig03_same, pars = c("alpha", "beta", "alpha_jk", "beta_jk", "sig_alpjk", "sig_betjk"), probs=c(0.1, 0.9))

stan_plot(cs10_sig0_same, pars = c("alpha", "beta"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest
stan_plot(cs10_sig03_same, pars = c("pjw", "pjm"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest





### For document
gen100_sig01_same = generate_yearly(
  K=100,  J=5,  Nm=200, Nw = 100,  
  alpha = rep(0,5), beta = rep(0,5),  
  sig_alpha = rep(0.1,5), sig_beta = rep(0.1,5),
  seed = 2863737    # seed for reproducibility
)

gen100_sig05_same = generate_yearly(
  K=100,  J=5,  Nm=200, Nw = 100,  
  alpha = rep(0,5), beta = rep(0,5),  
  sig_alpha = rep(0.5,5), sig_beta = rep(0.5,5),
  seed = 2863737    # seed for reproducibility
)

gen100_sig1_same = generate_yearly(
  K=100,  J=5,  Nm=200, Nw = 100,  
  alpha = rep(0,5), beta = rep(0,5),  
  sig_alpha = rep(1,5), sig_beta = rep(1,5),
  seed = 2863737    # seed for reproducibility
)


test_sym01= stan(
  model_code= cs_sym, 
  data=list(K=100, J=5, Mj = gen100_sig01_same$nm, Wj = gen100_sig01_same$nw, 
            a_sigalp=3, b_sigalp = 0.5, a_sigbet=3, b_sigbet = 0.5, 
            a_sigalpjk = 3, b_sigalpjk = 0.5,
            a_sigbetjk = 3, b_sigbetjk = 0.5,
            iter=20000))

test_sym05= stan(
  model_code= cs_sym, 
  data=list(K=100, J=5, Mj = gen100_sig05_same$nm, Wj = gen100_sig05_same$nw, 
            a_sigalp=3, b_sigalp = 0.5, a_sigbet=3, b_sigbet = 0.5, 
            a_sigalpjk = 3, b_sigalpjk = 0.5,
            a_sigbetjk = 3, b_sigbetjk = 0.5,
            iter=10000))

test_sym1= stan(
  model_code= cs_sym, 
  data=list(K=100, J=5, Mj = gen100_sig1_same$nm, Wj = gen100_sig1_same$nw, 
            a_sigalp=3, b_sigalp = 0.5, a_sigbet=3, b_sigbet = 0.5, 
            a_sigalpjk = 3, b_sigalpjk = 0.5,
            a_sigbetjk = 3, b_sigbetjk = 0.5,
            iter=10000))

plot01 = stan_plot(test_sym01, pars = c("alpha", "beta"))+
  xlim(-0.4,0.4) +
  geom_vline(xintercept=0, linewidth=1, col="blue")+
  ggtitle("Sigma = 0.1") # this is probably the clearest
plot05 = stan_plot(test_sym05, pars = c("alpha", "beta")) +
  xlim(-0.4,0.4) +
  geom_vline(xintercept=0, linewidth=1, col="blue")+
  ggtitle("Sigma = 0.5") # this is probably the clearest
plot1 = stan_plot(test_sym1, pars = c("alpha", "beta")) + 
  xlim(-0.4,0.4) +
  geom_vline(xintercept=0, linewidth=1, col="blue")+
  ggtitle("Sigma = 1") # this is probably the clearest

pdf(file=  "plot1.pdf")
plot1
dev.off()
