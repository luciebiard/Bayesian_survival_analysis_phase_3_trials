library(rstan)
library(plyr)
library(R2jags)
library(MCMCpack)
library(coda)
library(bayesplot)
library(survival)
library(purrr)
library(survminer)
library(tidyr)
library(broom)
library(dplyr)
source("functions_historical_data.R")


#########################################################
# Load all datasets
#########################################################
#reconstructed PRIMA historical dataset
dat0 <- read.table("prima_data.txt", header=T)
head(dat0)

# Current trial: final dataset 
# Due to confidentiality, the original data cannot be shared via open platfor
# Shared dataset are simulated based on published results of CLL7-SA trial
dat <-read.table("bazpower.txt", header=T)
head(dat)

# Current trial: interim dataset
dati <- read.table("bazpoweri.txt", header=T)
head(dati)

km_prima <- survfit(Surv(time, evt) ~ arm, dat0)
km_filo <- survfit(Surv(time, evt) ~ arm, dati)

par(mai=c(0.5, 0.6, 0.2,0.1))
plot(km_prima , las=1, bty="n", xaxt="n", yaxt="n", lty=1:2, col="gray")
axis(1,pos=0, at=seq(0, 60,6), tcl=-0.3, cex.axis=0.9, mgp=c(0.5, 0.5, 0.5))
axis(2, pos=0,at=seq(0,1, 0.2), tcl=-0.3, cex.axis=0.9, las=1, mgp=c(0.5, 0.5, 0.5))
points(c(0,summary(km_filo)$time[ summary(km_filo)$strata=="arm=0"]),  c(1,summary(km_filo)$surv[ summary(km_filo)$strata=="arm=0"]), type="s", lty=1)
points(c(0, summary(km_filo)$time[ summary(km_filo)$strata=="arm=1"]),c(1 ,summary(km_filo)$surv[ summary(km_filo)$strata=="arm=1"]), type="s", lty=2)
legend(0.1, 0.4, lty=c(1,2,1,2), col=c(1,1,"gray", "gray"), leg=c("CLL7-SA SoC", "CLL7-SA RTX", "External SoC", "External RTX"), bty="n", cex=0.9)
mtext(side=1, at=27, text="Months since randomization", cex=0.9, line=1)
mtext(side=2, at=0.5, text="PFS probability", cex=0.9, line=1.8)




##################################################################
##################################################################
# Early interim analysis after publication of PRIMA data: cut off 2011/1/1
##################################################################
##################################################################
# scaling times

echelle <- max(dat0$time)
dati$delai <- dati$time/max(dat0$time)
dat0$delai <- dat0$time/max(dat0$time)


# number of intervals and cutpoints
K <- max(5, min(sum(dati$evt)/8,20))
K

##################################################################
# Create dataset for piecewise-constant exponential model
##################################################################

test <- function_datalong(dati)
K
#Knot Placement/Location Calculations


cutpoints <-c(0,quantile(c(dati$delai[dati$evt==1]),seq(0,1,by=1/K))[-c(1,K+1)],1)
cutpoints <-c(0,quantile(c(dati$delai[dati$evt==1]),seq(0,1,by=1/K))[-c(1,K+1)],max(dati$delai))

dataset <- data.frame(trial=rep(1:2, c(nrow(dat0),nrow(dati))),
                       delai=c(dat0$delai, dati$delai),
                       evt=c(dat0$evt, dati$evt),
                       arm=c(dat0$arm, dati$arm))


dataset$TimeIntervals<- matrix(0,nrow=length(dataset$delai),K)  #amount of time spent in each interval
dataset$EventsIntervals<-matrix(0,nrow=length(dataset$delai),K)  #did a death occur during this interval for this patient?
dataset$CensIntervals<-matrix(0,nrow=length(dataset$delai),K)  #did censoring occur during this interval for this patient?
dataset$TrialArray<-matrix(0,nrow=length(dataset$delai),K)  
dataset$IntervalArray<-matrix(0,nrow=length(dataset$delai),K)  
dataset$InterventionArray<-matrix(0,nrow=length(dataset$delai),K)  

#Create the cutpoints, using the quantiles of the observed distribution of deaths (see Murray et al. 2014)
for (ii in 1:length(dataset$delai)){
  for (kk in 1:K){
    dataset$TimeIntervals[ii,kk]<-min(dataset$delai[ii],cutpoints[kk+1])-min(dataset$delai[ii],cutpoints[kk])
    if (dataset$evt[ii]==1){
      if(dataset$delai[ii]>cutpoints[kk] & dataset$delai[ii]<=cutpoints[kk+1]){
        dataset$EventsIntervals[ii,kk]<-1
      }
    }
    
    if (dataset$evt[ii]==0){
      if(dataset$delai[ii]>cutpoints[kk] & dataset$delai[ii]<=cutpoints[kk+1]){
        dataset$CensIntervals[ii,kk]<-1
      }
    }
    
    dataset$TrialArray[ii,kk]<-as.numeric(dataset$trial[ii])
    dataset$IntervalArray[ii,kk]<-kk
    dataset$InterventionArray[ii,kk]<-as.numeric(dataset$arm[ii])

  }  
}

#Create dataset in long format without observations for which interval is 0.
dataset_long<-data.frame(TimeIntervals=as.vector(dataset$TimeIntervals),
                         EventsIntervals=as.vector(dataset$EventsIntervals),
                         CensIntervals=as.vector(dataset$CensIntervals),
                         TrialArray=as.vector(dataset$TrialArray),
                         IntervalArray=as.vector(dataset$IntervalArray),
                         InterventionArray=as.vector(dataset$InterventionArray))
dataset_long <- dataset_long[dataset$TimeIntervals>0,]
dataset_long <- (ddply(dataset_long,.variables=.(TrialArray,IntervalArray,InterventionArray),summarise, TimeIntervals=sum(TimeIntervals),EventsIntervals=sum(EventsIntervals),CensIntervals=sum(CensIntervals)))





##########################################################################
# EXPONENTIAL MODEL
##########################################################################

##########################################################################
# Exclude historical data (keep current data only), a_0=0
##########################################################################

N <- nrow(dati)
X <- as.matrix(pull(dati, arm))
is_censored <- pull(dati,evt)==0
times <- pull(dati,delai)
msk_censored <- is_censored == 1
N_censored <- sum(msk_censored)

stan_data <- list(N_uncensored=N-N_censored, 
                  N_censored=N_censored, 
                  X_censored=as.matrix(X[msk_censored,]),
                  X_uncensored=as.matrix(X[!msk_censored,]),
                  times_censored=times[msk_censored],
                  times_uncensored = times[!msk_censored],
                  NC=ncol(X)
)
rstan_options(auto_write = TRUE)
expo_simple <- stan_model("CLL7SA/exponential.stan")
fit_0 <- sampling(expo_simple, data=stan_data, seed=42, chains=4, iter=10000, thin=5)


##########################################################################
# Pool historical data, a_0=1
##########################################################################

N <- nrow(dati)+nrow(dat0)
X <- as.matrix(c(dati$arm, dat0$arm))
is_censored <- I(c(dati$evt, dat0$evt)==0)
times <- c(dati$delai, dat0$delai)
msk_censored <- is_censored == 1
N_censored <- sum(msk_censored)

stan_data <- list(N_uncensored=N-N_censored, 
                  N_censored=N_censored, 
                  X_censored=as.matrix(X[msk_censored,]),
                  X_uncensored=as.matrix(X[!msk_censored,]),
                  times_censored=times[msk_censored],
                  times_uncensored = times[!msk_censored],
                  NC=ncol(X)
)


fit_1 <- sampling(expo_simple, data=stan_data, seed=42, chains=4, iter=10000, thin=5)




##########################################################################
# Conditional power prior
##########################################################################

N <- nrow(dati)
X <- as.matrix(pull(dati, arm))
is_censored <- pull(dati,evt)==0
times <- pull(dati,delai)
msk_censored <- is_censored == 1
N_censored <- sum(msk_censored)

N_ext <- nrow(dat0)
X_ext <- as.matrix(pull(dat0, arm))
is_censored_ext <- pull(dat0,evt)==0
times_ext <- pull(dat0,delai)
msk_censored_ext <- is_censored_ext == 1
N_censored_ext <- sum(msk_censored_ext)

stan_data <- list(N_uncensored=N-N_censored, 
                  N_uncensored_ext=N_ext-N_censored_ext, 
                  N_censored=N_censored, 
                  N_censored_ext=N_censored_ext, 
                  X_censored=as.matrix(X[msk_censored,]),
                  X_censored_ext=as.matrix(X_ext[msk_censored_ext,]),
                  X_uncensored=as.matrix(X[!msk_censored,]),
                  X_uncensored_ext=as.matrix(X_ext[!msk_censored_ext,]),
                  times_censored=times[msk_censored],
                  times_censored_ext=times_ext[msk_censored_ext],
                  times_uncensored = times[!msk_censored],
                  times_uncensored_ext = times_ext[!msk_censored_ext],
                  NC=ncol(X),
                  a0=0.1
)

rstan_options(auto_write = TRUE)
expo_cpp <- stan_model("exponential_CPP.stan") # compile stan model, may take a while


stan_data$a0 <- 0.75 # define desired power prior a_0 value
stan_data
fit_0.75 <- sampling(expo_cpp, data=stan_data, seed=42, chains=4, iter=10000, thin=5)



#diagnoses
check_hmc_diagnostics(fit)
monitor(fit)
obj_mcmc <- As.mcmc.list(fit)
np_fit <- nuts_params(fit)
posterior <- as.array(fit)

lapply(obj_mcmc, effectiveSize)
ratios_cp <- neff_ratio(fit)[1:2]
rhats <- rhat(fit)[1:2]

par(mfrow=c(1,4))
color_scheme_set("mix-brightblue-gray")
p1 <- mcmc_trace(posterior, pars = "betas[1]" , np = np_fit) + 
  xlab("Post-warmup iteration")+ggtitle("Exponential")
p2 <- mcmc_dens_overlay(posterior, pars = c("betas[1]"))
p3 <- mcmc_rhat(rhats) + yaxis_text(hjust = 1)
p4 <- mcmc_neff(ratios_cp, size = 2) + yaxis_text(hjust = 1)
multiplot(p1, p2, p3,p4,cols=4)



######################################################################################
# RESULTS DISPLAY
###################################################################################
library(tidybayes)

res_fit0 <- fit_0 %>%
  spread_draws(betas[1]) %>%
  dplyr::mutate(type = 'a0=0')
  
res_fit10 <- fit_0.1 %>%
  spread_draws(betas[1]) %>%
  dplyr::mutate(type = 'a0=0.1')

res_fit25 <- fit_0.25 %>%
  spread_draws(betas[1]) %>%
  dplyr::mutate(type = 'a0=0.25')

res_fit50 <- fit_0.5 %>%
  spread_draws(betas[1]) %>%
  dplyr::mutate(type = 'a0=0.5')

res_fit75 <- fit_0.75 %>%
  spread_draws(betas[1]) %>%
  dplyr::mutate(type = 'a0=0.75')


res_fit100 <- fit_1%>%
  spread_draws(betas[1]) %>%
  dplyr::mutate(type = 'a0=1')


all_cpp_data <-  res_fit0 %>%
  bind_rows(res_fit10, res_fit10, res_fit25, res_fit50, res_fit75, res_fit100)

all_cpp_data <- all_cpp_data %>%
  dplyr::mutate(betas = unlist(betas))


means = all_cpp_data %>%
  dplyr::group_by(type) %>%
  mean_qi(betas, .width = c(0.95,0.8, 0.5))
medians = all_cpp_data %>%
  dplyr::group_by(type) %>%
  median_qi(betas, .width = c(0))


  
  
all_cpp_data %>%
  ggplot(aes(y = betas, x = type)) +
  stat_interval() +
  scale_color_brewer(5)+
  geom_hline(yintercept=0, color='red') +
  geom_pointinterval(aes(y = betas, x=type), position = position_nudge(x = -0.1), data = means)+
  geom_pointinterval(aes(y = betas, x=type), position = position_nudge(x = 0), data = medians)+
 theme( legend.title = element_blank(), 
       legend.text = element_text(size=18), 
       plot.subtitle = element_text(size=14), 
       plot.title = element_text(size=18), 
       axis.title.x = element_text(size=20),
       axis.title.y = element_text(size=20),
       axis.text.x = element_text(size=20),
       axis.text.y = element_text(size=16),
       legend.position = "top") +
   xlab( expression(paste("Power prior parameter ", a[0])))+ylab("Posterior log(HR)")

 




xx <- as.data.frame(all_cpp_data)

tab <- cbind(tapply(xx$betas, xx$type, mean),
             tapply(xx$betas, xx$type, mean),
             tapply(xx$betas, xx$type, function(x){quantile(x, 0.025)}),
             tapply(xx$betas, xx$type, function(x){quantile(x, 0.975)}),
             tapply(xx$betas, xx$type, function(x){quantile(x, 0.9)}),
             tapply(xx$betas, xx$type, function(x){sum(x<log(0.6))/length(x)}))
             
xtable(round(tab,3), digits=3)




