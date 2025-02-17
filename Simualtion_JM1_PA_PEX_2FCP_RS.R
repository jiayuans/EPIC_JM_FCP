#########################################################################################
# Joint model for PA and PE (Left truncation) - Simulation Code
# Author: Jiayuan Shi
#
# Description:
#   This script simulates data and implements a joint model (JM1) for Pseudomonas aeruginosa (PA) and pulmonary
#   exacerbations (PE) in cystic fibrosis patients, incorporating left truncation.
#   This model includes two fixed change points and with random intercept and slopes.
#
# Components:
#   - Data preprocessing for PA and PE cohorts
#   - Simulation of event times for a non-homogeneous Poisson process (NHPP)
#   - Specification of a Bayesian joint model in JAGS
#   - MCMC execution and results analysis
#   - Calculation of WAIC, DIC, and convergence diagnostics
#
# Input:
#   - Data-PA-cohort.csv
#   - Data-multiple-final-cohort-EPICstart-Bt-LengthPEx.csv
#
# Output:
#   - Model estimates and convergence diagnostics
#   - WAIC and DIC values for model comparison
#   - Results saved to CSV files and trace plots
#########################################################################################

# Load required libraries
library(coda)
library(rjags)
library(runjags)
library(tidyverse)
library(mcmcplots)

#########################################################################################
# Resembled the EPIC study data
#########################################################################################

# Define directory
dirg <- "C:/UCHealth/RA/Project/EPIC-CF/Analysis_Jiayuan/EPIC/"
setwd(dirg)

# Load PA data
dat.pa <- read.csv(file=paste("C:/UCHealth/RA/Project/EPIC-CF/Data/Source","/Data-PA-cohort.csv",sep=""))
dat.pa <- dat.pa[,-1] # Remove first column (unnecessary index)

sexf <- dat.pa[!duplicated(dat.pa$cffidno,dat.pa$sexf),5]

# Use similar follow-up time from the first 400 subjects
dat.pa1 <-aggregate(dat.pa$VisitAge, by=list(dat.pa$cffidno),
                    FUN=max, na.rm=TRUE)
names(dat.pa1) <- c('cffidno','age.max')

dat.pa2 <-aggregate(dat.pa$age.min, by=list(dat.pa$cffidno),
                    FUN=max, na.rm=TRUE)
names(dat.pa2) <- c('cffidno','age.min')

first.t<-dat.pa2$age.min
last.t<-dat.pa1$age.max
first.tt<-first.t[c(1:401)]
last.tt<-last.t[c(1:401)]
sexf<-sexf[c(1:401)]

first.tt<-first.tt[-359]
last.tt<-last.tt[-359]

long.time <- cbind(first.tt,last.tt)
first.tt <- long.time[,1]
last.tt <- long.time[,2]

# Number of unique individuals
N<-length(last.tt)
# Participant ID
id<-rep(1:N)
length(id)

# Time of first visit and last visit
t<-round(first.tt)
tt<-round(last.tt)
k.pa<-(tt-t)*4
kk=max(k.pa)

# set number of iterations
I=201

# Set true values
c0=-4.5
c1=0.1
c2=0.2 
c3=0.1
c4=0.1
Verror=1
cp1.true=4.5
cp2.true=14.5

# Set seed for reproducible code
set.seed(123)

#########################################################################################
# Function that generates observations from a NHPP- returns event times
# Input: parameters for the mean of a poisson process: a(shape parameter),b, T (exposure time)
# Output: Returns the event times and a variable that indicates whether the observation is an event or a 
#		  censoring time (no events observed in the whole interval); status=1 indicates event and 0 censoring
#########################################################################################

# Function to generate event times for NHPP
NHPP<-function(a,b,T){
  mu <- b*T^a  # Mean of the Poisson process up to time T
  n <-rpois(1, mu)  # Number of events poisson
  if (n!=0) {
    u <- runif(n,0,1)
    u <- sort(u)
    y <- T*u^(1/a) 
    y[length(y)+1] <- T
    y_0 <- rep(NA,length(y))
    for (i in 2:length(y_0)){
      y_0[i] <- y[i-1]
    }
    y_0[which(is.na(y_0)==TRUE)] <- 0
    return(cbind(y_0,y,c(rep(1,length(y)-1),0),n))    #returns n event times
  } else 
    return(cbind(0,T,0,n)) 
}

#########################################################################################
# Function that creates an event times dataset for a poisson process (continuous data )
# Input: parameters for the intensity function alpha; beta; beta0; x; ga (association parameter); Tei 
# Output: A dataset with variables
# 		 id, xi (treatment),Tei, time, status
# -------------- Building the simulated poisson data -----
#########################################################################################

poisson.d <- function(alpha,beta,beta0,x,ga,TTei){
  le <- length(x)
  c_0i <- rnorm(le,0,1) #1.3
  vi <- exp(ga*b_0i+c_0i)

  times <- NHPP(b=vi[1]*exp(beta*x[1])*exp(beta0),a=alpha,T=TTei[1])
  start <-  times[,1]
  stop <- times[,2]
  status <- times[,3]
  n.rec <- times[,4]
  id <- rep(1,length(stop))
  xi <- rep(x[1],length(stop))
  Tei <- rep(TTei[1],length(stop))
  for (i in 2:length(x)){
    times2 <- NHPP(b=vi[i]*exp(beta0+beta*x[i]),a=alpha,T=TTei[i]) 
    start2 <-  times2[,1]
    stop2 <- times2[,2]
    status2 <- times2[,3]
    n.rec2 <- times2[,4]
    id <- c(id,rep(i,length(stop2)))
    xi <- c(xi,rep(x[i],length(stop2)))
    Tei <- c(Tei,rep(TTei[i],length(stop2)))
    
    start <- c(start,start2)
    stop <- c(stop,stop2)
    status <- c(status,status2)
    n.rec <- c(n.rec,n.rec2)
  }
  return(data.frame(id,xi,Tei,n.rec,start,stop,status))
}

# This code is only for one simulation
  b_0i<-rnorm(N,0,1) 
  X1=c(rep(1,N/2),rep(0,N/2))

  I1<-matrix(NA, nrow=N, ncol=kk, byrow=TRUE)
  I2<-matrix(NA, nrow=N, ncol=kk, byrow=TRUE)
  p2<-matrix(NA, nrow=N, ncol=kk, byrow=TRUE)
  Y<-matrix(NA, nrow=N, ncol=kk, byrow=TRUE)
  X<-matrix(NA, nrow=N, ncol=kk, byrow=TRUE)
 
# Simulation of X: visit age for PA in matrix
   for (i in 1:N){
    X[i,1:k.pa[i]]<-c(seq(t[i],tt[i]-0.25,0.25))
  }
  
# Simulation of Y: binary PA outcome in matrix 
  for (i in 1:N){
    for (j in 1:k.pa[i]){
      I1[i,j]<-ifelse(X[i,j]< cp1.true,-1,1)
      I2[i,j]<-ifelse(X[i,j]< cp2.true,-1,1)
      p2[i,j]=exp(c0+c1*(X[i,j]-cp1.true)+c2*(X[i,j]-cp1.true)*I1[i,j]+c3*(X[i,j]-cp2.true)*I2[i,j]+c4*X1[i]+b_0i[i])/(1+exp(c0+c1*(X[i,j]-cp1.true)+c2*(X[i,j]-cp1.true)*I1[i,j]+c3*(X[i,j]-cp2.true)*I2[i,j]+c4*X1[i]+b_0i[i]))
      Y[i,j]=rbinom(Verror, 1, p2[i,j])
    }
  }

# Simulation of recurrent PE events
  simdat.pe00 <- poisson.d(alpha=1.8,beta=0.25,beta0=-4.5,x=X1,ga=0.25,TTei=tt-0.25)
  simdat.pe00 <- as.data.frame(simdat.pe00)
 
  timeS <- as.data.frame(cbind(id,t)) ## left truncation time
  timeE <- as.data.frame(cbind(id,tt))
  
  simdat.pe0 <- merge(simdat.pe00, timeS,all=TRUE)
  simdat.pe <- subset(simdat.pe0, stop >= t)
  
  time <- subset(simdat.pe,status==1)
  time1 <- time[,c("id","stop")]  
  simdat.pe1 <- merge(timeS,timeE,all=TRUE)
  simdat.pe2 <- merge(simdat.pe1,time1,all=TRUE)
  simdat.pe2$stop[which(is.na(simdat.pe2$stop))] <- 0
  
  count <- simdat.pe2 %>% count(id)
  max.count <- max(count$n) 
  
  # Assign unique time points per subject
  simdat.pe3 <- simdat.pe2 %>% group_by(id) %>% mutate(time = c(1:length(id)))
  
  # Create a complete time grid for all subjects
  Yd.temp <- data.frame(id = rep(unique(simdat.pe00$id),each=max.count), time = 1:max.count) 
  Y.epic <- merge(simdat.pe3,Yd.temp,by=c('id','time'),all.y=TRUE)
  
  # Convert PE time data to matrix form
  Ti <- matrix(Y.epic$stop, N, max.count, byrow=TRUE)
  
  # Start time and end time for PE
  time.t0 <- t
  time.tau <- tt
  
  # Get number of observations per subject
  sum.na <- rep(NA,N)
  k.pe=rep(NA,N)
  
  ids <- unique(Y.epic$id) ## 103104 103125 103129 103145 103147
  for (i in 1:N){
    na.indices <- which(Y.epic$t[Y.epic$id==ids[i]] %in% NA)
    if (length(na.indices)==0){
      k.pe[i] <- max.count} else{
        k.pe[i] <- min(na.indices)-1}
  }
  
  
  #########################################################################################
  # Model Specification: Joint Model (PA and PE) with 2 fixed change point in PA
  #########################################################################################
  
  modelfixcp <- "
data { 
  for(i in 1:N){
       zeros[i]<- 0
  }
}
model { 
  for(i in 1:N){ 
        for(j in 1:k.pa[i]){
  ### PA model
        Y[i,j] ~ dbin(p2[i,j],1)
        logit(p2[i,j]) <- c0 + c[1] * (X[i,j]-cp1) + c[2] * (X[i,j]-cp1) * (2*step(X[i,j]-cp1)-1) + c[3] * (X[i,j]-cp2) * (2*step(X[i,j]-cp2)-1) + c[4] * X1[i] + u[i]
        }
        for(j in 1:k.pe[i]){
  ### PE model
       ## Weibull baseline
        lambda0[i,j] <- a*(Ti[i,j])^(a-1)
        lambda[i,j] <- lambda0[i,j]*v[i]*exp(b0+b*X1[i])
       }
        u[i] ~ dnorm(0,u.tau)
        L.a[i] <- prod(((p2[i,1:k.pa[i]])^(Y[i,1:k.pa[i]]))*((1-p2[i,1:k.pa[i]])^(1-Y[i,1:k.pa[i]])))
        ll.a[i] <- log(L.a[i])
        w[i] ~ dnorm(0,w.tau)
        v[i] <- exp(ga*u[i]+w[i])
        L.e[i] <- ifelse(Ti[i,1]!=0, prod(lambda[i,1:k.pe[i]]) * exp(v[i]*exp(b0+b*X1[i])*(time.t0[i]^a-time.tau[i]^a)), exp(v[i]*exp(b0+b*X1[i])*(time.t0[i]^a-time.tau[i]^a)))
        ll.e[i] <- log(L.e[i])
        phi[i] <- -log(L.e[i]) + 1000
        zeros[i] ~ dpois(phi[i])
  }
  log_lik0.a <- sum(ll.a[])
  log_lik0.e <- sum(ll.e[]) 
  dev.a <- -2*log_lik0.a
  dev.e <- -2*log_lik0.e
  c0 ~ dnorm(0,0.0001)
	for (k in 1:4){
	      c[k] ~ dnorm(0,0.0001)	
	}
  ## prior distributions
	u.tau ~ dgamma(0.001,0.001)
	cp1 ~ dnorm(cp1.mu,cp1.tau)	
	cp2.temp ~ dunif(0,max)
	cp2 <- cp1 + cp2.temp
	cp1.mu ~ dnorm(0,0.001)
	cp1.tau ~ dgamma(0.001,0.001)
	B1 <-c[1]-c[2]-c[3]
  B2 <-c[1]+c[2]-c[3]
  B3 <-c[1]+c[2]+c[3]
  u.tau.inv <- 1/u.tau  ## variance 
  a ~ dgamma(0.01,0.01)
  b0 ~ dnorm(0,0.0001)	
  b ~ dnorm(0,0.0001)		
	ga ~ dnorm(0,0.0001)
	w.tau ~ dgamma(0.001,0.001)
	w.tau.inv <- 1/w.tau  ## variance 
}"

#########################################################################################
# Run Bayesian Model using JAGS
#########################################################################################

# Observed DATA
data <- dump.format(list(X=X, Y=Y, N=N, k.pa=k.pa, max=max(tt),
                         X1=X1, k.pe=k.pe, time.t0=time.t0, time.tau=time.tau, Ti=Ti)) 
# Define initial Values
inits1 <- dump.format(list(c0=-4.6, c=c(0.1,0.17,0.1,0.1), u.tau=1, cp1=4.5, cp2.temp=10,
                           b0=-4.5, b=0.25, a=1.8, w.tau=1, ga=0.25,
                           .RNG.name="base::Super-Duper", .RNG.seed=1))
inits2 <- dump.format(list(c0=-4.61, c=c(0.1,0.17,0.1,0.1)+0.01, u.tau=1, cp1=4.6, cp2.temp=10,
                           b0=-4.51, b=0.26, a=1.81, w.tau=1, ga=0.26,
                           .RNG.name="base::Super-Duper", .RNG.seed=2))

# Run MCMC
res <- run.jags(model=modelfixcp, burnin=20000, sample=1000, 
                monitor=c("B1","B2","B3","cp1","cp2","c0","c","u.tau.inv",
                          "b0","b","a","ga","w.tau.inv","u","v","w",
                          "u.tau","w.tau","cp1.mu","cp1.tau","cp2.temp","ll.a","ll.e","dev.a","dev.e","dic"), 
                data=data, n.chains=2, inits=c(inits1,inits2), thin=10, module='dic')

# Results and traceplots
summary <- summary(res)
result_df <- as.data.frame(summary)

res_jm <- res$mcmc
vars<-mcmc.list(res_jm[[1]][,c(1:16)],res_jm[[2]][,c(1:16)])
traplot(vars)



