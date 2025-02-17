#########################################################################################
# Joint model for PA+ and PEX (Left truncation)
# Author: Jiayuan Shi
#
# Description:
#   This script implements a joint model (JM1) for Pseudomonas aeruginosa (PA) and pulmonary
#   exacerbations (PE) in cystic fibrosis patients, incorporating left truncation.
#   This model includes two fixed change points and with random intercept and slopes.
#
# Components:
#   - WAIC computation functions
#   - Data preprocessing for PA and PE cohorts
#   - Bayesian joint model specification using JAGS
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
#   - Results saved to CSV files
#########################################################################################

# Clear environment
rm(list=ls())
# Load required libraries
library(coda)
library(rjags)
library(runjags)
library(mcmcplots)
library(tidyverse)

#########################################################################################
# Functions to calculate WAIC
#########################################################################################

# Function to compute variance of individual likelihood contributions
colVars <- function (a){
  diff <- a - matrix (colMeans(a), nrow(a), ncol(a), byrow=TRUE)
  vars <- apply(diff^2,2,sum)/(nrow(a)-1)
  return (vars)
}

# Function to compute WAIC based on likelihood samples
# Return lppd, p_waic_1, p_waic_2, and waic, which we define
# as 2*(lppd - p_waic_2), as recommmended in BDA
waic <- function (log_lik){
  lppd <- sum (log (colMeans(exp(log_lik))))
  p_waic_1 <- 2*sum (log(colMeans(exp(log_lik))) - colMeans(log_lik))
  p_waic_2 <- sum (colVars(log_lik))
  waic_2 <- -2*lppd + 2*p_waic_2
  return (list (waic=waic_2, p_waic=p_waic_2, lppd=lppd, p_waic_1=p_waic_1))
}

# Function to compute WAIC from a CODA object
# input: CODA object; from = index where the iterations for the likelihood start; 
#         to = index where the iterations for the likelihood end
# output: WAIC
waicf <- function(mcmc.obj,from,to){
  log.like1 <- mcmc.obj[[1]][,c(from:to)]
  log.like2 <- mcmc.obj[[2]][,c(from:to)]
  log.like <- rbind(log.like1,log.like2)
  
  waico <- waic(log.like)
  return(waico)
}

#########################################################################################
# Data Processing: PA cohort
#########################################################################################

# Define directory
dirg <- "Z:/EJCStudents/ShiJ/EPIC-CF/Data/Source/"

# Load PA data
dat.pa <- read.csv(file=paste(dirg,"/Data-PA-cohort.csv",sep=""))
dat.pa <- dat.pa[,-1] # Remove first column (unnecessary index)
head(dat.pa, n=10)

# Number of unique individuals
M <- length(unique(dat.pa$cffidno))
length(dat.pa$cffidno)

count.pa <- dat.pa %>% count(cffidno)

# Assign unique time points per subject
dat.pa1 <- dat.pa[order(dat.pa$cffidno) , ]
dat.pa2 <- dat.pa1 %>% group_by(cffidno) %>% mutate(time = c(1:length(cffidno)))

# Create a complete time grid for all subjects
Yd.temp.pa <- data.frame(cffidno = rep(unique(dat.pa$cffidno),each=126), time = 1:126)
head(Yd.temp.pa,n=12)

Y.epic.pa <- merge(dat.pa2,Yd.temp.pa,by=c('cffidno','time'),all.y=TRUE)
nrow(dat.pa)
nrow(Y.epic.pa)
str(Y.epic.pa)
head(Y.epic.pa[,1:6],150 )

# Prepare covariates
X.dat.pa <- dat.pa1[!duplicated(dat.pa1$cffidno,dat.pa1$sexf,dat.pa1$mut), ]
X.dat.pa1 <- model.matrix(cltpa ~ sexf + as.factor(mut), data = X.dat.pa)
X1_pa <- as.numeric(X.dat.pa1[,2]) ## sexf: female
X2_pa <- as.numeric(X.dat.pa1[,3]) ## mut2: 1 alleles D F508
X3_pa <- as.numeric(X.dat.pa1[,4]) ## mut3: 0 alleles D F508

# Convert longitudinal PA data to matrix form
X <- matrix(Y.epic.pa$VisitAge, M, 126, byrow=TRUE)
Y <- matrix(Y.epic.pa$cltpa, M, 126, byrow=TRUE)

# Get number of observations per subject
k.pa=rep(NA,M)

ids.pa <- unique(Y.epic.pa$cffidno) ## 103104 103125 103129 103145 103147
for (i in 1:M){
  na.indices.pa <- which(Y.epic.pa$VisitAge[Y.epic.pa$cffidno==ids.pa[i]] %in% NA)
  if (length(na.indices.pa)==0){
    k.pa[i] <- 126} else{
      k.pa[i] <- min(na.indices.pa)-1}
}

#########################################################################################
# Data Processing: PE cohort
#########################################################################################

# Load PEX data
dat.pe0 <- read.csv(file=paste(dirg,"/Data-multiple-final-cohort-EPICstart-Bt-LengthPEx.csv",sep=""))
dat.pe0 <- dat.pe0[with(dat.pe0, order(cffidno, tstart)), ]
length(unique(dat.pe0$cffidno)) #1734

# Extract relevant columns
dat.pe <- dat.pe0[, c(1,5:7,21,24,25)]
summary(dat.pe)
head(dat.pe, n=10)

# Compute time intervals for each subject
timeS <- aggregate(dat.pe$tstart, by=list(dat.pe$cffidno),
                   FUN=min, na.rm=TRUE)
colnames(timeS) <- c("cffidno","time.t0")
timeE <-aggregate(dat.pe$tend, by=list(dat.pe$cffidno),
                  FUN=max, na.rm=TRUE)
colnames(timeE) <- c("cffidno","time.tau")

X.dat.pe <- dat.pe[!duplicated(dat.pe$cffidno,dat.pe$sexf,dat.pe$mut), ]
X.dat.pe1 <- model.matrix(status ~ sexf + as.factor(mut) , data = X.dat.pe)

time <- subset(dat.pe,status==1)
time$t <- time$tend
time1 <- time[,c(1,8)]  
dat.pe1 <- merge(timeS,timeE)
dat.pe2 <- merge(dat.pe1,time1,all=TRUE)
dat.pe2$t[which(is.na(dat.pe2$t))] <- 0

length(which(dat.pe2$t== 0))

# Number of unique individuals
M <- length(unique(dat.pe2$cffidno)) 
count.pe <- dat.pe2 %>% count(cffidno)
max(count.pe$n)

# Assign unique time points per subject
dat.pe3 <- dat.pe2 %>% group_by(cffidno) %>% mutate(time = c(1:length(cffidno)))

# Create a complete time grid for all subjects
Yd.temp.pe <- data.frame(cffidno = rep(unique(dat.pe$cffidno),each=43), time = 1:43) 
head(Yd.temp.pe,n=12)

Y.epic.pe <- merge(dat.pe3,Yd.temp.pe,by=c('cffidno','time'),all.y=TRUE)
nrow(dat.pe)
nrow(Y.epic.pe)
str(Y.epic.pe)

# Convert PE time data to matrix form
Ti <- matrix(Y.epic.pe$t, M, 43, byrow=TRUE)
Ti[1:2,]##need to remove first row of column headers
dim(Ti)
Ti[1,]

# Prepare covariates
X1 <- as.numeric(X.dat.pe1[,2]) ## sexf: female
X2 <- as.numeric(X.dat.pe1[,3]) ## mut2: 1 alleles D F508
X3 <- as.numeric(X.dat.pe1[,4]) ## mut3: 0 alleles D F508
##X4 <- as.numeric(X.dat.pe1[,5]) ## rcaucasian: white

# Store survival times
time.t0 <- timeS$time.t0
time.tau <- timeE$time.tau

# Get number of observations per subject
k.pe=rep(NA,M)

ids.pe <- unique(Y.epic.pe$cffidno) ## 103104 103125 103129 103145 103147
for (i in 1:M){
  na.indices.pe <- which(Y.epic.pe$t[Y.epic.pe$cffidno==ids.pe[i]] %in% NA)
  if (length(na.indices.pe)==0){
    k.pe[i] <- 43} else{
      k.pe[i] <- min(na.indices.pe)-1}
}


#########################################################################################
# Model Specification: Joint Model (PA and PE) with 2 fixed change point in PA
#########################################################################################

modelfixcp <- "
data { 
  for(i in 1:M){
       zeros[i]<- 0
  }
}
model { 
  for(i in 1:M){ 
        for(j in 1:k.pa[i]){
  ### PA model
        Y[i,j] ~ dbin(p2[i,j],1)
        logit(p2[i,j]) <- c0 + (c[1]+u1[i]) * (X[i,j]-cp1) + (c[2]+u2[i]) * (X[i,j]-cp1) * (2*step(X[i,j]-cp1)-1) + (c[3]+u3[i]) * (X[i,j]-cp2) * (2*step(X[i,j]-cp2)-1) + c[4] * X1_pa[i] + c[5] * X2_pa[i] + c[6] * X3_pa[i] + u[i]
        }
        for(j in 1:k.pe[i]){
  ### PE model
       ## Weibull baseline
        lambda0[i,j] <- a*(Ti[i,j])^(a-1)
        lambda[i,j] <- lambda0[i,j]*v[i]*exp(b0+b[1]*X1[i]+b[2]*X2[i]+b[3]*X3[i])
       }
        u[i] ~ dnorm(0,u.tau)
        u1[i] ~ dnorm(0,u.tau1)
        u2[i] ~ dnorm(0,u.tau2)
        u3[i] ~ dnorm(0,u.tau3)
        L.a[i] <- prod(((p2[i,1:k.pa[i]])^(Y[i,1:k.pa[i]]))*((1-p2[i,1:k.pa[i]])^(1-Y[i,1:k.pa[i]])))
        ll.a[i] <- log(L.a[i])
        w[i] ~ dnorm(0,w.tau)
        v[i] <- exp(ga*u[i]+w[i]+ga1*u1[i]+ga2*u2[i]+ga3*u3[i])
        L.e[i] <- ifelse(Ti[i,1]!=0, prod(lambda[i,1:k.pe[i]]) * exp(v[i]*exp(b0+b[1]*X1[i]+b[2]*X2[i]+b[3]*X3[i])*(time.t0[i]^a-time.tau[i]^a)), exp(v[i]*exp(b0+b[1]*X1[i]+b[2]*X2[i]+b[3]*X3[i])*(time.t0[i]^a-time.tau[i]^a)))
        ll.e[i] <- log(L.e[i])
        phi[i] <- -log(L.e[i]) + 1000
        zeros[i] ~ dpois(phi[i])
  }
  log_lik0.a <- sum(ll.a[]) 
  log_lik0.e <- sum(ll.e[]) 
  dev.a <- -2*log_lik0.a
  dev.e <- -2*log_lik0.e
  ## prior distributions
  c0 ~ dnorm(0,0.0001)
	for (k in 1:6){
	      c[k] ~ dnorm(0,0.0001)	
	}
	u.tau ~ dgamma(0.001,0.001)
	u.tau1 ~ dgamma(0.001,0.001)
	u.tau2 ~ dgamma(0.001,0.001)
	u.tau3 ~ dgamma(0.001,0.001)
	cp1 ~ dnorm(cp1.mu,cp1.tau)	
	cp2.temp ~ dunif(0,max)
	cp2 <- cp1 + cp2.temp
	cp1.mu ~ dnorm(0,0.001)
	cp1.tau ~ dgamma(0.001,0.001)
	B1 <-c[1]-c[2]-c[3]
  B2 <-c[1]+c[2]-c[3]
  B3 <-c[1]+c[2]+c[3]
  u.tau.inv <- 1/u.tau  ## variance 
  u.tau1.inv <- 1/u.tau1  ## variance 
  u.tau2.inv <- 1/u.tau2  ## variance 
  u.tau3.inv <- 1/u.tau3  ## variance 
  a ~ dgamma(0.01,0.01)
  b0 ~ dnorm(0,0.0001)	
  for (p in 1:3){
	     b[p] ~ dnorm(0,0.0001)		
  }
	ga ~ dnorm(0,0.0001)
	ga1 ~ dnorm(0,0.0001)
	ga2 ~ dnorm(0,0.0001)
	ga3 ~ dnorm(0,0.0001)
	w.tau ~ dgamma(0.001,0.001)
	w.tau.inv <- 1/w.tau  ## variance 
}"

#########################################################################################
# Run Bayesian Model using JAGS
#########################################################################################

# Observed DATA 
data <- dump.format(list(X=X, Y=Y, M=M, max=max(dat.pa$VisitAge), k.pa=k.pa, X1_pa=X1_pa,X2_pa=X2_pa,X3_pa=X3_pa,
                         X1=X1,X2=X2,X3=X3,k.pe=k.pe, time.t0=time.t0, time.tau=time.tau, Ti=Ti)) 
# Define initial Values
inits1 <- dump.format(list(c0=-4, c=c(0.1,0.1,0.1,0.1,0.1,0.1), u.tau=1, u.tau1=1, u.tau2=1, u.tau3=1, cp1=5, cp2.temp=9,
                           b0=0, b=c(0.25,-0.14,-0.38), a=1, w.tau=1, 
                           .RNG.name="base::Super-Duper", .RNG.seed=1))
inits2 <- dump.format(list(c0=-4.1, c=c(0.1,0.1,0.1,0.1,0.1,0.1)+0.01, u.tau=1, u.tau1=1, u.tau2=1, u.tau3=1, cp1=5.1, cp2.temp=9.1,
                           b0=0.1,b=c(0.25,-0.14,-0.38)+0.1, a=1.1,  w.tau=1, 
                           .RNG.name="base::Super-Duper", .RNG.seed=2))
# Run MCMC
res <- run.jags(model=modelfixcp, burnin=4000, sample=6000, 
                 monitor=c("B1", "B2","B3","c0", "c", "cp1", "cp2","u","u1","u2","u3",
                           "u.tau","u.tau.inv","u.tau1","u.tau1.inv","u.tau2","u.tau2.inv","u.tau3","u.tau3.inv",
                           "cp1.mu","cp1.tau", "cp2.temp",
                           "b0","b", "a","v","ga","ga1","ga2","ga3","w","w.tau","w.tau.inv","ll.a","ll.e","dev.a","dev.e","dic"), 
                 data=data, n.chains=2, inits=c(inits1,inits2), thin=10,  module='dic')

sum <- summary(res)
sum

# Save results
write.csv(sum,"Z:/EJCStudents/ShiJ/EPIC-CF/Result/JM1_PA_PEX_2FCP_RS.csv")

sum(sum[,11]<1.1)
sum.df <- round(as.data.frame(sum),2)

res_jm <- res$mcmc
dimnames(res_jm[[1]])
vars<-mcmc.list(res_jm[[1]][,c(1:9)],res_jm[[2]][,c(1:9)])
str(vars)
plot(vars[,1])

summary(vars)
traplot(vars)

## tests of convergence
gew<-geweke.diag(vars)

# Compute WAIC and DIC for model comparison
# Calculate WAIC
waic.jm.a <- waicf(res_jm,which(rownames(sum.df)=="ll.a[1]"),which(rownames(sum.df)=="ll.a[1734]"))
print(waic.jm.a)
##$waic
##[1] 28292.02

##$p_waic
##[1] 1183.228

##$lppd
##[1] -12962.78

##$p_waic_1
##[1] 813.2396
waic.a <- waic.jm.a$waic

waic.jm.e <- waicf(res_jm,which(rownames(sum.df)=="ll.e[1]"),which(rownames(sum.df)=="ll.e[1734]"))
print(waic.jm.e)
##$waic
##[1] 10844.37

##$p_waic
##[1] 605.1683

##$lppd
##[1] -4817.015

##$p_waic_1
##[1] 399.4678
waic.e <- waic.jm.e$waic

# Mean deviance
dicf <- function(mcmc.obj,from,to){
  log.like1 <- mcmc.obj[[1]][,c(from:to)]
  log.like2 <- mcmc.obj[[2]][,c(from:to)]
  log.like <- rbind(log.like1,log.like2)
  
  mean.dev <- -2*sum (colMeans(log.like)) 
  return(mean.dev) 
}
md.jm.a <- dicf(res_jm,which(rownames(sum.df)=="ll.a[1]"),which(rownames(sum.df)=="ll.a[1734]")) 
print(md.jm.a) # PA:26738.8

md.jm.e <- dicf(res_jm,which(rownames(sum.df)=="ll.e[1]"),which(rownames(sum.df)=="ll.e[1734]")) 
print(md.jm.e) # PE:10033.5

## calculate DIC = mean D + pd = 2*mean D - D mean
sum.df <- as.data.frame(sum)
c0 <- sum.df[4,4]
c1 <- sum.df[5,4]
c2 <- sum.df[6,4]
c3 <- sum.df[7,4] 
c4 <- sum.df[8,4] 
c5 <- sum.df[9,4] 
c6 <- sum.df[10,4] 
cp1 <- sum.df[11,4] 
cp2 <- sum.df[12,4] 
u <- sum.df[13:1746,4]
u1 <- sum.df[1747:3480,4]
u2 <- sum.df[3481:5214,4]
u3 <- sum.df[5215:6948,4]

b0 <- sum.df[6960,4]
b1 <- sum.df[6961,4]
b2 <- sum.df[6962,4]
b3 <- sum.df[6963,4] 
b <- c(b1,b2,b3)
a <- sum.df[6964,4]
v <- sum.df[6965:8698,4]

ga <- sum.df[8699,4]
ga1 <- sum.df[8700,4]
ga2 <- sum.df[8701,4]
ga3 <- sum.df[8702,4] 
w <- sum.df[8703:10436,4]

p2 <- matrix(NA, M, 126, byrow=TRUE)
I1<-matrix(NA, M, 126, byrow=TRUE)
I2<-matrix(NA, M, 126, byrow=TRUE)
lambda0 <- matrix(NA, M, 43, byrow=TRUE)
lambda <- matrix(NA, M, 43, byrow=TRUE)
L.a <- rep(NA,M)
L.e <- rep(NA,M)  
ll.a <- rep(NA,M)
ll.e <- rep(NA,M)  

for(i in 1:M){
  for(j in 1:k.pa[i]){
    I1[i,j]<-ifelse(X[i,j] < cp1,-1,1)
    I2[i,j]<-ifelse(X[i,j] < cp2,-1,1)
    p2[i,j]=exp(c0+(c1+u1[i])*(X[i,j]-cp1)+(c2+u2[i])*(X[i,j]-cp1)*I1[i,j]+(c3+u3[i])*(X[i,j]-cp2)*I2[i,j]+c4*X1_pa[i]+c5*X2_pa[i]+c6*X3_pa[i]+u[i])/(1+exp(c0+(c1+u1[i])*(X[i,j]-cp1)+(c2+u2[i])*(X[i,j]-cp1)*I1[i,j]+(c3+u3[i])*(X[i,j]-cp2)*I2[i,j]+c4*X1_pa[i]+c5*X2_pa[i]+c6*X3_pa[i]+u[i]))
  }
  for(j in 1:k.pe[i]){
    lambda0[i,j] <- a*(Ti[i,j])^(a-1)
    lambda[i,j] <- lambda0[i,j]*v[i]*exp(b0+b[1]*X1[i]+b[2]*X2[i]+b[3]*X3[i])
  }
  L.a[i] <- prod(((p2[i,1:k.pa[i]])^(Y[i,1:k.pa[i]]))*((1-p2[i,1:k.pa[i]])^(1-Y[i,1:k.pa[i]])))
  ll.a[i] <- log(L.a[i])
  L.e[i] <- ifelse(Ti[i,1]!=0, prod(lambda[i,1:k.pe[i]]) * exp(v[i]*exp(b0+b[1]*X1[i]+b[2]*X2[i]+b[3]*X3[i])*(time.t0[i]^a-time.tau[i]^a)), exp(v[i]*exp(b0+b[1]*X1[i]+b[2]*X2[i]+b[3]*X3[i])*(time.t0[i]^a-time.tau[i]^a)))
  ll.e[i] <- log(L.e[i])
}

log_lik0.a <- sum(ll.a[]) 
dev_hat.a <- -2*log_lik0.a

log_lik0.e <- sum(ll.e[]) 
dev_hat.e <- -2*log_lik0.e

dic.a <- md.jm.a*2-dev_hat.a ##DIC:28404.49
dic.a
dic.e <- md.jm.e*2-dev_hat.e ##DIC:10850.61
dic.e

total.dic <- dic.a+dic.e
total.dic
total.waic <- waic.a+waic.e
total.waic

cbind(round(dic.a),round(waic.a),round(dic.e),round(waic.e), round(total.dic),round(total.waic))

#########################################################################################
# End of Script
#########################################################################################
