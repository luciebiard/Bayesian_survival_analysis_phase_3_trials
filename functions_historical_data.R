function_datalong <- function(dati){
  
  K <- max(5, min(sum(dati$evt)/8,20))
  K <- floor(K)
  
  
  
  cutpoints <- c(0,quantile(dati$time[dati$evt==1],seq(0,1,by=1/K)[-c(1,K+1)]),max(dati$time))
  cutpoints <- unique(cutpoints)
  K <- length(cutpoints)-1
  
  dataset <- dati
  dataset$TimeIntervals<- matrix(0,nrow=length(dataset$time),K)  #amount of time spent in each interval
  dataset$EventsIntervals<-matrix(0,nrow=length(dataset$time),K)  #did a death occur during this interval for this patient?
  dataset$CensIntervals<-matrix(0,nrow=length(dataset$time),K)  #did censoring occur during this interval for this patient?
  dataset$IntervalArray<-matrix(0,nrow=length(dataset$time),K)  
  dataset$InterventionArray<-matrix(0,nrow=length(dataset$time),K)  
  
  
  
  for (ii in 1:length(dataset$time)){
    for (kk in 1:K){
      dataset$TimeIntervals[ii,kk]<-min(dataset$time[ii],cutpoints[kk+1])-min(dataset$time[ii],cutpoints[kk])
      if (dataset$evt[ii]==1){
        if(dataset$time[ii]>cutpoints[kk] & dataset$time[ii]<=cutpoints[kk+1]){
          dataset$EventsIntervals[ii,kk]<-1
        }
      }
      
      if (dataset$evt[ii]==0){
        if(dataset$time[ii]>cutpoints[kk] & dataset$time[ii]<=cutpoints[kk+1]){
          dataset$CensIntervals[ii,kk]<-1
        }
      }
      
      
      dataset$IntervalArray[ii,kk]<-kk
      dataset$InterventionArray[ii,kk]<-as.numeric(dataset$arm[ii])
      
    }  
  }
  
  #Create dataset in long format without observations for which interval is 0.
  dataset_long<-data.frame(TimeIntervals=as.vector(dataset$TimeIntervals),
                           EventsIntervals=as.vector(dataset$EventsIntervals),
                           CensIntervals=as.vector(dataset$CensIntervals),
                           IntervalArray=as.vector(dataset$IntervalArray),
                           InterventionArray=as.vector(dataset$InterventionArray))
  dataset_long <- dataset_long[dataset$TimeIntervals>0,]
  dataset_long <- (ddply(dataset_long,.variables=.(IntervalArray,InterventionArray),summarise, TimeIntervals=sum(TimeIntervals),EventsIntervals=sum(EventsIntervals),CensIntervals=sum(CensIntervals)))
  
  
  
  return(list(K=K, cutpoints=cutpoints, dataset_long))
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



get_sim_pem_one <- function(zz , alpha, beta, theta=0, cutpoints, maxtime=NULL, event.only=T){
  
  y_hat_time <- y_hat_event <- vector("numeric", length=length(zz))
  
  K <- length(alpha)
  
  t_dur <- diff(cutpoints)
  
  if(event.only) {
    borne_temps <- maxtime+1
    
    while(borne_temps>maxtime){
      u <- runif(1, 0, 1)
      uu <- -log(1-u)
      
      ft_dur <- c(0, exp(as.numeric(alpha)+ beta*zz + theta)*t_dur)
      
      ii <- ( uu - cumsum(ft_dur) > 0)
      
      if(all(ii)){
        
        y_hat_time <- borne_temps <- maxtime+1
        y_hat_event <- 0
        
      }else{
        ii <- which(ii==F)[1] - 1
        y_hat_time <- borne_temps <- cutpoints[ii] + (uu -  cumsum(ft_dur)[ii])/exp(as.numeric(alpha[ii])+ beta*zz + theta)
        y_hat_event <- 1
      }
    }
  } else {
  
   
      u <- runif(1, 0, 1)
      uu <- -log(1-u)
      cc <- runif(1, 6, maxtime) # administrative censoring mimicking uniform trial accrual
    
      ft_dur <- c(0, exp(as.numeric(alpha)+ as.numeric(beta)*zz + theta)*t_dur)
    
      ii <- ( uu - cumsum(ft_dur) > 0)
    
      if(all(ii)){
    
        y_hat_time <- min(maxtime, cc)
        y_hat_event <- 0
    
      }else{
        ii <- which(ii==F)[1] - 1
        y_hat_time <- cutpoints[ii] + (uu -  cumsum(ft_dur)[ii])/exp(as.numeric(alpha[ii])+ beta*zz + theta)
        
        y_hat_event <- ifelse(y_hat_time<cc, 1,0)
        
        y_hat_time <- min(cc, y_hat_time)
        
      }
    }
    
  return(c(y_hat_time=y_hat_time, y_hat_event=y_hat_event, arm=zz))

}


get_sim_pem <- function(zz , alpha, beta, theta=0, cutpoints, maxtime, event.only=T){
  
  tmp <- lapply(zz, function(x){get_sim_pem_one(x, alpha, beta, theta, cutpoints, maxtime, event.only)})
  

  return(data.frame(y_hat_time=map_dbl(tmp, 1), y_hat_event=map_dbl(tmp, 2), arm=map_dbl(tmp, 3)))
  
}



get_pp_data <- function(mat_mcmc, cutpoints, arm, maxtime, modele=NULL, event.only=T){
  K <- length(cutpoints)-1
  
  if(is.null(modele)){
    res <- apply(mat_mcmc,1, function(x){get_sim_pem(arm, alpha=x[1:K], beta=x[K+1], cutpoints=cutpoints, maxtime=maxtime, event.only=event.only)})
  } else {
    res <- apply(mat_mcmc,1, function(x){get_sim_pem(arm, alpha=x[1:K], beta=x[K+1], theta=x[K+5], cutpoints=cutpoints, maxtime=maxtime, event.only=event.only)})
  }
  return(res)
}


plot_km_pp_check <- function(km, km0, mat_mcmc, cutpoints,arm=dataset$arm, maxtime, ncurves=100, niter=4000, titre, cex.title=0.8, modele=NULL){
  library(RColorBrewer)
  coul <- brewer.pal(n = 11, name = "RdBu")[c(2,5, 7, 11)]

  res <- get_pp_data(mat_mcmc, cutpoints,arm=dataset$arm, modele=modele, maxtime=maxtime)
  
  res2 <- lapply(res, function(x){survfit(Surv(x$y_hat_time, x$y_hat_event)~x$arm)})


  plot(km, col=coul[c(1,4)], 
       ylab="Survival probability", xlab="Time from randomization", cex.lab=0.8, cex.axis=0.8, main=titre, cex.main=cex.title, 
       tcl=-0.2, mgp=c(1.5, 0.3, 0))
  
  lapply(res2, function (x){
    
    points(c(0,summary(x)$time[summary(x)$strata=='x$arm=0']),c(1, summary(x)$surv[summary(x)$strata=='x$arm=0']), type="s",  col=coul[2])
    points(c(0,summary(x)$time[summary(x)$strata=='x$arm=1']), c(1,summary(x)$surv[summary(x)$strata=='x$arm=1']), type="s",  col=coul[3])
    
    
  })

   points(c(0,summary(km0)$time[summary(km0)$strata=='arm=0']), c(1, summary(km0)$surv[summary(km0)$strata=='arm=0']), type="s", col=coul[1], lwd=1.5, lty=2)
  points(c(0,summary(km0)$time[summary(km0)$strata=='arm=1']), c(1, summary(km0)$surv[summary(km0)$strata=='arm=1']), type="s", col=coul[4], lwd=1.5, lty=2)
  
  points(c(0,summary(km)$time[summary(km)$strata=='arm=0']), c(1, summary(km)$surv[summary(km)$strata=='arm=0']), type="s", col=coul[1], lwd=1.5)
  points(c(0,summary(km)$time[summary(km)$strata=='arm=1']), c(1, summary(km)$surv[summary(km)$strata=='arm=1']), type="s", col=coul[4], lwd=1.5)
  
  
  
  legend(0.05, 0.4, c("Ctrl arm, observed","Ctrl arm, historical","Ctrl arm, posterior predicted", "RTX arm, observed","RTX arm, historical", "RTX arm, posterior predicted"), col=coul[c(1,1,2,4,4,3)],lty=c(1,2,1,1,2,1),  bty="n", cex=0.8, lwd=1.5)

}


get_sim_expo_one <- function(zz , alpha, beta,maxtime=NULL){
  
  y_hat_time <- y_hat_event <- vector("numeric", length=length(zz))
  
    
    tt <- rexp(1, exp(alpha + beta*zz))
    
    
    cc <- runif(1, 0, maxtime) # administrative censoring mimicking uniform trial accrual
    
    
    y_hat_time <- pmin(tt, cc)
    y_hat_event <- ifelse(tt<cc, 1, 0)
    
  return(c(y_hat_time=y_hat_time, y_hat_event=y_hat_event, arm=zz))
  
}


get_sim_expo <- function(zz , alpha, beta,  maxtime, event.only=T){
  
  tmp <- lapply(zz, function(x){get_sim_expo_one(x, alpha, beta, maxtime)})
  
  
  return(data.frame(y_hat_time=map_dbl(tmp, 1), y_hat_event=map_dbl(tmp, 2), arm=map_dbl(tmp, 3)))
  
}



get_pp_data_expo <- function(mat_mcmc, arm, maxtime, modele=NULL){
  
  
  
  
    res <- apply(mat_mcmc,1, function(x){get_sim_expo(arm, alpha=x[1], beta=x[2], maxtime=maxtime)})
  
  return(res)
}

