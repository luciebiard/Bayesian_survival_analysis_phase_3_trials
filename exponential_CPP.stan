data {
    int<lower=1> N_uncensored;                                      
    int<lower=1> N_uncensored_ext;  
    int<lower=1> N_censored;                                        
    int<lower=1> N_censored_ext;
    int<lower=0> NC;                                    
    matrix[N_censored, NC] X_censored;
    matrix[N_censored_ext, NC] X_censored_ext;
    matrix[N_uncensored, NC] X_uncensored;                       
    matrix[N_uncensored_ext,NC] X_uncensored_ext;
    vector<lower=0>[N_censored] times_censored;                        
    vector<lower=0>[N_censored_ext] times_censored_ext;
    vector<lower=0>[N_uncensored] times_uncensored;
    vector<lower=0>[N_uncensored_ext] times_uncensored_ext;
    real<lower=0, upper=1> a0;
}
parameters {
    vector[NC] betas;
    real intercept;                            
}
model {
    betas ~ normal(0,2);
    intercept ~ normal(-5,2);

    target += a0*exponential_lpdf(times_uncensored_ext | exp(intercept+X_uncensored_ext*betas));
    target += a0*exponential_lccdf(times_censored_ext | exp(intercept+X_censored_ext*betas));

    target += exponential_lpdf(times_uncensored | exp(intercept+X_uncensored*betas));
    target += exponential_lccdf(times_censored | exp(intercept+X_censored*betas));
}

generated quantities {
    vector[N_uncensored] times_uncensored_sampled;
    for(i in 1:N_uncensored) {
        times_uncensored_sampled[i] = exponential_rng(exp(intercept+X_uncensored[i,]*betas));
    }
}
