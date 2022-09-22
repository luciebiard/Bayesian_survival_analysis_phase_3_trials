In this respository, you may find the code to estimate models such as those presented in “Bayesian survival analysis for early detection of treatment effects in phase 3 clinical trials” paper in Contemporary Clinical Trials Communications (2020).

Here are some information on how to use these files:

    cll_github.R: Code for analyses corresponding to motivating example 1 in the above mentioned paper
    
    functions_historical_data: data-management and plotting functions 

    Pooled_piecewise.stan: Stan code for the piecewise model in motivating example 1 in the above mentioned paper

    exponential.stan: Stan code for the exponential model used in motivating example 2 (without power prior)

    exponential_CPP.stan: Stan code for the power prior exponential model used in motivating example 2

    bazpower.txt: Simulated dataset mimicking final data from current data in motivating example 2

    bazpoweri.txt: Simulated dataset mimicking interim data from current data in motivating example 2

    prima_data.txt: Reconstructed dataset from PRIMA trial publication (Salles et al - Lancet 2010;377:42-51) using Guyot et al's method (Guyot et al - BMCMRM 2012,12:9) 

The work was performed using R 3.5.1 with rstan version 2.18.2, and rstanarm development version inclusing stan_surv() function. 
For any question, please contact lucie.biard@-paris.fr.
