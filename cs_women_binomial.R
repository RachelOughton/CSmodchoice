### Binomial version, just to test that I'm not going mad
### This works exactly as I'd hoped!

### Then building up the multinomial model with just one group (ie. no genders just alpha)


library(ggplot2)
library(gridExtra)


cs_binom = 
  "
  data {
    int<lower=1> K;             // Number of years
    int<lower=0> Y[K];         // Binomial outcomes
    int<lower=0> N[K];         // Yearly totals
    real<lower=0> a_sigpi;   // IG prior for sig_pi
    real<lower=0> b_sigpi;
    real<lower=0> a_p;   // beta prior for p
    real<lower=0> b_p; 
  }

  parameters {
    real<lower=0, upper=1> p; // vector of intercept parameters
    real<lower=0, upper=1> pi[K];   // vector of women coefficient/parameters
    real<lower=0> sig_pi;           // Should this be K dimensions? No I don't think so. J when I go multinomial  
  }
  
/// Check to here!
  model {
    target += inv_gamma_lpdf(sig_pi|a_sigpi,b_sigpi);

    p ~ beta(a_p,b_p);
    for (k in 1:K){
      pi[k] ~ normal(p, sig_pi);          // Number of men students per module
      Y[k] ~ binomial(N[k], pi[k]);
    }
  }
"

## I want some where the pi are almost (or exactly) identical (ie. sig_pi=0) 
## and some where it varies a lot
## Plan of attack
### Specify sig_pi and p
### Generate the pi from N(p, sig_pi)
### Generate the Y from Binom(N, pi[k])

binom_generate = function(
  k, # years
  p,   # base p
  sig_pi,   # SD of pi)
  N         # K-vector of Ns for each year
){
  pi = rnorm(k, mean = p, sd=sig_pi)
  pi[which(pi>1)] = 1
  pi[which(pi<0)] = 0
  
  Y = sapply(
    1:k, 
    function(k){
      rbinom(1, size=N[k], prob = pi[k])
    }
  )
  
  out_df = data.frame(Y=Y, pi=pi, N = N)
  return(out_df)
}

gen20_sig0 = binom_generate(k=20, p=0.6, sig_pi=0, N=rep(100,20))
gen20_sig02 = binom_generate(k=20, p=0.6, sig_pi=0.2, N=rep(100,20))
gen20_sig04 = binom_generate(k=20, p=0.6, sig_pi=0.4, N=rep(100,20))

bin20_sig0= stan(
  model_code= cs_binom, 
  data=list(K=20, J=6, Y = gen20_sig0$Y, N= gen20_sig0$N, a_sigpi=3, b_sigpi=1.5, a_p=1, b_p=1), 
  iter=5000)
stan_plot(bin20_sig0, pars = c("p", "pi", "sig_pi"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest

bin20_sig02= stan(
  model_code= cs_binom, 
  data=list(K=20, J=6, Y = gen20_sig02$Y, N= gen20_sig02$N, a_sigpi=3, b_sigpi=1.5, a_p=1, b_p=1), 
  iter=5000)
stan_plot(bin20_sig02, pars = c("p", "pi", "sig_pi"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest

bin20_sig04= stan(
  model_code= cs_binom, 
  data=list(K=20, J=6, Y = gen20_sig04$Y, N= gen20_sig04$N, a_sigpi=3, b_sigpi=1.5, a_p=1, b_p=1), 
  iter=5000)
stan_plot(bin20_sig04, pars = c("p", "pi", "sig_pi"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest

### This shows exactly what I'd hoped
### Maybe I need to look at the posteriors for the pjw and pjm? 
## Ie. not by year, but by module? But allowing an individual p for each year

grid.arrange(
  stan_plot(bin20_sig0, pars = c("p", "sig_pi"))+xlim(0,1), 
  stan_plot(bin20_sig02, pars = c("p", "sig_pi"))+xlim(0,1), 
  stan_plot(bin20_sig04, pars = c("p", "sig_pi"))+xlim(0,1)
  )


#### Now multinomial, but with only one group

cs_yearly_one = 
  "
  data {
    int<lower=1> J;             // Number of modules
    int<lower=1> K;             // Number of years
    int<lower=0> Nj[K,J];         // Total on module over K years (years as rows)
//    real<lower=0> sig_alpha;   // sd for alpha - choose reasonably large, like 1
    real<lower=0> a_sigalp;
    real<lower=0> b_sigalp;
    real<lower=0> a_sigalpjk;   // sd for sig_alpjk (same for all J?)
    real<lower=0> b_sigalpjk; 
  }
  

  parameters {
    real alpha[J]; // vector of intercept parameters
    real alpha_jk[K,J];  // intercept parameters for each module and year
    real<lower=0> sig_alpjk[J];
    real<lower=0> sig_alpha;
  }
  
  transformed parameters {
    matrix[K,J] pjn;
    matrix[K,J] denom_n_vec;
    vector[K] denom_n;

    for (k in 1:K){
      for (j in 1:J){
        denom_n_vec[k,j] = exp(alpha_jk[k,j]);
      }
      denom_n[k] = sum(denom_n_vec[k]);
    }

    for (k in 1:K){
      for (j in 1:J){
        pjn[k,j] = exp(alpha_jk[k,j])/denom_n[k];
      }
    }
    
  }
  model {
    target += inv_gamma_lpdf(sig_alpha|a_sigalp,b_sigalp);
    target += inv_gamma_lpdf(sig_alpjk|a_sigalpjk,b_sigalpjk);


    alpha ~ normal(0, sig_alpha);         // Priors for each module

    for (k in 1:K){
      alpha_jk[k] ~ normal(alpha, sig_alpjk);  // Priors for each module and year
      Nj[k] ~ multinomial(to_vector(pjn[k]));          // Number of men students per module
    }
  }
"

cs20one_sig0_same= stan(
  model_code= cs_yearly_one, 
  data=list(K=20, J=5, Nj = gen20_sig0_same$nm, 
            a_sigalp=3, b_sigalp=1.5,
            a_sigalpjk = 3, b_sigalpjk = 1.5,
            iter=80000))

cs20one_sig03_same= stan(
  model_code= cs_yearly_one, 
  data=list(K=20, J=5, Nj = gen20_sig03_same$nm,
            a_sigalp=3, b_sigalp=1.5,
            a_sigalpjk = 3, b_sigalpjk = 1.5,
            iter=80000))


gen20_sig05_same = generate_yearly(
  K=20,  J=5,  Nm=100, Nw = 50,  
  alpha = rep(0,5), beta = rep(0,5),  
  sig_alpha = rep(0.5,5), sig_beta = rep(0.5,5),
  seed = 2863737    # seed for reproducibility
)



print(cs20one_sig0_same, pars = c("alpha", "alpha_jk", "sig_alpha", "sig_alpjk"), probs=c(0.1, 0.9))
print(cs20one_sig03_same, pars = c("alpha", "alpha_jk", "sig_alpha", "sig_alpjk"),probs=c(0.1, 0.9))

## This works much better (more restrictive IG distribution)

cs20onetp_sig0_same= stan(
  model_code= cs_yearly_one, 
  data=list(K=20, J=5, Nj = gen20_sig0_same$nm, 
            a_sigalp=3, b_sigalp=0.5,
            a_sigalpjk = 3, b_sigalpjk = 0.5,
            iter=80000))

cs20onetp_sig03_same= stan(
  model_code= cs_yearly_one, 
  data=list(K=20, J=5, Nj = gen20_sig03_same$nm,
            a_sigalp=3, b_sigalp=0.5,
            a_sigalpjk = 3, b_sigalpjk = 0.5,
            iter=80000))


cs20onetp_sig05_same= stan(
  model_code= cs_yearly_one, 
  data=list(K=20, J=5, Nj = gen20_sig05_same$nm,
            a_sigalp=3, b_sigalp=1,
            a_sigalpjk = 3, b_sigalpjk = 1,
            iter=80000))




print(cs20one_sig0_same, pars = c("alpha", "alpha_jk", "sig_alpha", "sig_alpjk"), probs=c(0.1, 0.9))
print(cs20one_sig03_same, pars = c("alpha", "alpha_jk", "sig_alpha", "sig_alpjk"),probs=c(0.1, 0.9))

stan_plot(cs20one_sig0_same, pars = c("alpha", "alpha_jk", "sig_alpha", "sig_alpjk"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest
stan_plot(cs20one_sig03_same, pars = c("alpha", "alpha_jk", "sig_alpha", "sig_alpjk"))+geom_vline(xintercept=0, linewidth=1, col="blue") # this is probably the clearest


#### This is working too!!!!

grid.arrange(
  stan_plot(cs20onetp_sig0_same, pars = c("alpha")) + xlim(-0.5,0.5), 
  stan_plot(cs20onetp_sig03_same, pars = c("alpha")) + xlim(-0.5,0.5),
  stan_plot(cs20onetp_sig05_same, pars = c("alpha")) + xlim(-0.5,0.5)
)

