########################################################################
#####################################################################
#####################################################################
## Name: Joint model for PA+ and PEX (Left truncation)
## Author: Jiayuan Shi
## Date: 
##
## Description:  

##			  
##
## Components:
## 		- 
##		  
## Input:  Data-PA-cohort,
##    Data-multiple-final-cohort-EPICstart-Bt-LengthPEx.csv
##		
## Output:
##		- 
## Notes:
## 	
## 
rm(list=ls())
library(coda)
library(rjags)
library(runjags)
library(mcmcplots)
library(tidyverse)

#########################################################################################
#########################################################################################
##   Functions to calculate WAIC
## 
#########################################################################################
##### calculating the variance of each individual contribution to the likelihood
colVars <- function (a){
  ##a <- log_lik
  diff <- a - matrix (colMeans(a), nrow(a), ncol(a), byrow=TRUE)
  vars <- apply(diff^2,2,sum)/(nrow(a)-1)
  return (vars)
}

##### Calculation of Waic. Returns lppd, p_waic_1, p_waic_2, and waic, which we define
## as 2*(lppd - p_waic_2), as recommmended in BDA
waic <- function (log_lik){
  lppd <- sum (log (colMeans(exp(log_lik))))
  p_waic_1 <- 2*sum (log(colMeans(exp(log_lik))) - colMeans(log_lik))
  p_waic_2 <- sum (colVars(log_lik))
  waic_2 <- -2*lppd + 2*p_waic_2
  return (list (waic=waic_2, p_waic=p_waic_2, lppd=lppd, p_waic_1=p_waic_1))
}

##### function to read the CODA object and calculate the WAIC
## input: CODA object; from = index where the iterations for the likelihood start; 
##         to = index where the iterations for the likelihood end
## output: WAIC
waicf <- function(mcmc.obj,from,to){
  log.like1 <- mcmc.obj[[1]][,c(from:to)]
  log.like2 <- mcmc.obj[[2]][,c(from:to)]
  log.like <- rbind(log.like1,log.like2)
  
  waico <- waic(log.like)
  return(waico)
}


#############################################################################
##################use matrix data #################
##################read in data for PA+ #################
##################Reading in main dataset##########################
#################################################################

dirg <- "Z:/EJCStudents/ShiJ/EPIC-CF/Data/Source/"
dat.pa <- read.csv(file=paste(dirg,"/Data-PA-cohort.csv",sep=""))
dat.pa <- dat.pa[,-1]
head(dat.pa, n=10)

dat.pa_f1 <-aggregate(dat.pa$VisitAge, by=list(dat.pa$cffidno),
                      FUN=min, na.rm=TRUE)
names(dat.pa_f1) <- c('cffidno','min.aage')
head(dat.pa_f1)

dat.pa_f2 <-aggregate(dat.pa$VisitAge, by=list(dat.pa$cffidno),
                      FUN=max, na.rm=TRUE)
names(dat.pa_f2) <- c('cffidno','max.aage')
head(dat.pa_f2)

dat.pa_f <- merge(dat.pa_f1, dat.pa_f2, by=c('cffidno'))
dat.pa_f$fage.a <-  dat.pa_f$max.aage - dat.pa_f$min.aage


#############################################################################
##################use matrix data #################
##################read in data for PEX #################
##################Reading in main dataset##########################
#################################################################

dat.pe0 <- read.csv(file=paste(dirg,"/Data-multiple-final-cohort-EPICstart-Bt-LengthPEx.csv",sep=""))
dat.pe0 <- dat.pe0[with(dat.pe0, order(cffidno, tstart)), ]
length(unique(dat.pe0$cffidno)) 

dat.pe <- dat.pe0[, c(1,5:7,21,24,25)]
summary(dat.pe)
head(dat.pe, n=10)

dat.pe_f1 <-aggregate(dat.pe$tstart, by=list(dat.pe$cffidno),
                      FUN=min, na.rm=TRUE)
names(dat.pe_f1) <- c('cffidno','min.eage')
head(dat.pe_f1)
nrow(dat.pe_f1)
dat.pe_f2 <-aggregate(dat.pe$tend, by=list(dat.pe$cffidno),
                      FUN=max, na.rm=TRUE)
names(dat.pe_f2) <- c('cffidno','max.eage')
head(dat.pe_f2)

dat.pe_f <- merge(dat.pe_f1, dat.pe_f2, by=c('cffidno'))
dat.pe_f$fage.e <-  dat.pe_f$max.eage - dat.pe_f$min.eage


dat_f <- merge(dat.pa_f, dat.pe_f, by='cffidno')
dat_f$min.age <- pmax(dat_f$min.aage, dat_f$min.eage, na.rm = TRUE)
dat_f$max.age <- pmin(dat_f$max.aage, dat_f$max.eage, na.rm = TRUE)
dat_f$fage <-  dat_f$max.age - dat_f$min.age

dat_f1 <- dat_f[,c('cffidno',"max.age","min.age","fage")]
dat.pa1 <- merge(dat.pa,dat_f1,by="cffidno")
dat.pa2 <- subset(dat.pa1, VisitAge>=min.age & VisitAge<=max.age)

dat.pe1 <- merge(dat.pe,dat_f1,by="cffidno")
dat.pe2 <- subset(dat.pe1, tstart>=min.age & tend<=max.age)

# minimum fu at least 5
dat.pa_sub <- subset(dat.pa2, fage>=5) 
length(unique(dat.pa_sub$cffidno)) #1415

dat.pe_sub <- subset(dat.pe2, fage>=5) 
length(unique(dat.pe_sub$cffidno)) #1415


##########################PA datast##########################
dat.pa <- dat.pa_sub
M <- length(unique(dat.pa$cffidno))
length(dat.pa$cffidno)

count.pa <- dat.pa %>% count(cffidno)
max(count.pa$n)
##########################Assigning unique number to each subject##########################
dat.pa1 <- dat.pa[order(dat.pa$cffidno) , ]
dat.pa2 <- dat.pa1 %>% group_by(cffidno) %>% mutate(time = c(1:length(cffidno)))

Yd.temp.pa <- data.frame(cffidno = rep(unique(dat.pa$cffidno),each=116), time = 1:116)
head(Yd.temp.pa,n=12)

Y.epic.pa <- merge(dat.pa2,Yd.temp.pa,by=c('cffidno','time'),all.y=TRUE)
nrow(dat.pa)
nrow(Y.epic.pa)
str(Y.epic.pa)
head(Y.epic.pa[,1:6],150 )

#################Readingin data for X vectors#############################
X.dat.pa <- dat.pa1[!duplicated(dat.pa1$cffidno,dat.pa1$sexf,dat.pa1$mut), ]
X.dat.pa1 <- model.matrix(cltpa ~ sexf + as.factor(mut), data = X.dat.pa)
X1_pa <- as.numeric(X.dat.pa1[,2]) ## sexf: female
X2_pa <- as.numeric(X.dat.pa1[,3]) ## mut2: 1 alleles D F508
X3_pa <- as.numeric(X.dat.pa1[,4]) ## mut3: 0 alleles D F508

#################Readingin data for VisitAge matrix#############################

X <- matrix(Y.epic.pa$VisitAge, M, 116, byrow=TRUE)

#################Readingin data for Pa+ matrix#############################

Y <- matrix(Y.epic.pa$cltpa, M, 116, byrow=TRUE)
Y[1:2,]##need to remove first row of column headers
dim(Y)

#################input variables for simulation#####################
#### checking for how many individuals we have NAs in the middle of followup

k.pa=rep(NA,M)

ids.pa <- unique(Y.epic.pa$cffidno) ## 103104 103125 103129 103145 103147
for (i in 1:M){
  na.indices.pa <- which(Y.epic.pa$VisitAge[Y.epic.pa$cffidno==ids.pa[i]] %in% NA)
  if (length(na.indices.pa)==0){
    k.pa[i] <- 116} else{
      k.pa[i] <- min(na.indices.pa)-1}
}


##########################PE datast##########################
dat.pe <- dat.pe_sub
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
time1 <- time[,c(1,11)]  
dat.pe1 <- merge(timeS,timeE)
dat.pe2 <- merge(dat.pe1,time1,all=TRUE)
dat.pe2$t[which(is.na(dat.pe2$t))] <- 0

length(which(dat.pe2$t== 0))

M <- length(unique(dat.pe2$cffidno)) 
count.pe <- dat.pe2 %>% count(cffidno)
max(count.pe$n)
##########################Assigning unique number to each subject##########################
dat.pe3 <- dat.pe2 %>% group_by(cffidno) %>% mutate(time = c(1:length(cffidno)))

Yd.temp.pe <- data.frame(cffidno = rep(unique(dat.pe$cffidno),each=42), time = 1:42) 
head(Yd.temp.pe,n=12)

Y.epic.pe <- merge(dat.pe3,Yd.temp.pe,by=c('cffidno','time'),all.y=TRUE)
nrow(dat.pe)
nrow(Y.epic.pe)
str(Y.epic.pe)

#################Readingin data for time matrix#############################
Ti <- matrix(Y.epic.pe$t, M, 42, byrow=TRUE)
Ti[1:2,]##need to remove first row of column headers
dim(Ti)
Ti[1,]

#################Readingin data for X, time.t0, time.tau vectors#############################
X1 <- as.numeric(X.dat.pe1[,2]) ## sexf: female
X2 <- as.numeric(X.dat.pe1[,3]) ## mut2: 1 alleles D F508
X3 <- as.numeric(X.dat.pe1[,4]) ## mut3: 0 alleles D F508
##X4 <- as.numeric(X.dat.pe1[,5]) ## rcaucasian: white

time.t0 <- timeS$time.t0
time.tau <- timeE$time.tau

#################input variables for simulation#####################
#### checking for how many individuals we have NAs in the middle of followup

k.pe=rep(NA,M)

ids.pe <- unique(Y.epic.pe$cffidno) ## 103104 103125 103129 103145 103147
for (i in 1:M){
  na.indices.pe <- which(Y.epic.pe$t[Y.epic.pe$cffidno==ids.pe[i]] %in% NA)
  if (length(na.indices.pe)==0){
    k.pe[i] <- 42} else{
      k.pe[i] <- min(na.indices.pe)-1}
}


#################################################################
####### Joint model finding 2 fixed change point for PA+
####### rec event model for PEX 
##---------------------------------------------------------------         
#################################################################
#### 
modelrancp <- "
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
        logit(p2[i,j]) <- c0 + (c[1]+u1[i]) * (X[i,j]-cp1) + (c[2]+u2[i]) * (X[i,j]-cp1) * (2*step(X[i,j]-cp1)-1) + (c[3]+u3[i]) * (X[i,j]-cp2) * (2*step(X[i,j]-cp2)-1)  + c[4] * X1_pa[i] + c[5] * X2_pa[i] + c[6] * X3_pa[i] + u[i]
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
  c0 ~ dnorm(0,0.0001)
	for (k in 1:6){
	      c[k] ~ dnorm(0,0.0001)	
	}
  ## prior distributions
	u.tau ~ dgamma(0.001,0.001)
	u.tau1 ~ dgamma(0.001,0.001)
	u.tau2 ~ dgamma(0.001,0.001)
	u.tau3 ~ dgamma(0.001,0.001)
	cp1 ~ dunif(0,max)
	cp2 ~ dunif(cp1,max)
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


####Observed DATA 
data <- dump.format(list(X=X, Y=Y, M=M, max=max(dat.pa$VisitAge), k.pa=k.pa, X1_pa=X1_pa,X2_pa=X2_pa,X3_pa=X3_pa,
                         X1=X1,X2=X2,X3=X3,k.pe=k.pe, time.t0=time.t0, time.tau=time.tau, Ti=Ti)) 
###initial Values
inits1 <- dump.format(list(c0=-4, c=c(0.1,0.1,0.1,0.1,-0.1,-0.4), u.tau=1, u.tau1=1, u.tau2=1, u.tau3=1, cp1=5, cp2=15,
                           b0=-4, b=c(0.2,-0.2,-0.4), a=1.7, w.tau=1, 
                           .RNG.name="base::Super-Duper", .RNG.seed=1))
inits2 <- dump.format(list(c0=-4.1, c=c(0.1,0.1,0.1,0.1,-0.1,-0.4)+0.01, u.tau=1, u.tau1=1, u.tau2=1, u.tau3=1, cp1=5.1, cp2=15.1,
                           b0=-4.1,b=c(0.2,-0.2,-0.4)+0.1, a=1.8,  w.tau=1, 
                           .RNG.name="base::Super-Duper", .RNG.seed=2))
#### Run the model and produce plots
res <- run.jags(model=modelrancp, burnin=1000, sample=5000, 
                 monitor=c("B1", "B2","B3","cp1", "cp2","c0", "c", 
                           "u.tau.inv","u.tau1.inv","u.tau2.inv","u.tau3.inv",
                           "b0","b", "a","ga","ga1","ga2","ga3","w.tau.inv",
                           "u.tau","u.tau1","u.tau2","u.tau3","w.tau",
                           "u","u1","u2","u3", "v","w","ll.a","ll.e","dev.a","dev.e"), 
                 data=data, n.chains=2, method = "parallel",inits=c(inits1,inits2), thin=8)

sum <- summary(res)
sum

write.csv(sum,"Z:/EJCStudents/ShiJ/EPIC-CF/Result/Fixed_CP/restricted/JM1_PA_PEX_2FCP_RS_restricted.csv")
save(res, file="Z:/EJCStudents/ShiJ/EPIC-CF/Result/Fixed_CP/restricted/JM1_PA_PEX_2FCP_RS_restricted.RData") 

load("Z:/EJCStudents/ShiJ/EPIC-CF/Result/Fixed_CP/restricted/JM1_PA_PEX_2FCP_RS_restricted.RData")
sum <- read.csv("Z:/EJCStudents/ShiJ/EPIC-CF/Result/Fixed_CP/restricted/JM1_PA_PEX_2FCP_RS_restricted.csv")

sum(sum[,12]<1.1)
rownames(sum) <- sum$X
sum.df <- round(as.data.frame(sum[,2:12]),2)

library(dplyr)
sum.df <- sum.df %>%
  mutate(Mean_CI = paste0(Mean, " (", Lower95, ", ", Upper95, ")")) %>%
  select(Mean_CI, everything()) 
write.csv(sum.df,"Z:/EJCStudents/ShiJ/EPIC-CF/Result/Fixed_CP/JM1_PA_PEX_2FCP_RS_restricted_R.csv")


res_jm <- res$mcmc
dimnames(res_jm[[1]])
vars<-mcmc.list(res_jm[[1]][,c(1:29)],res_jm[[2]][,c(1:29)])
str(vars)
plot(vars[,1])

summary(vars)
traplot(vars)


## tests of convergence
gew<-geweke.diag(vars)

## calculate WAIC
waic.jm.a <- waicf(res_jm,which(rownames(sum.df)=="ll.a[1]"),which(rownames(sum.df)=="ll.a[1415]"))
print(waic.jm.a)
waic.a <- waic.jm.a$waic

waic.jm.e <- waicf(res_jm,which(rownames(sum.df)=="ll.e[1]"),which(rownames(sum.df)=="ll.e[1415]"))
print(waic.jm.e)
waic.e <- waic.jm.e$waic

## mean deviance
dicf <- function(mcmc.obj,from,to){
  log.like1 <- mcmc.obj[[1]][,c(from:to)]
  log.like2 <- mcmc.obj[[2]][,c(from:to)]
  log.like <- rbind(log.like1,log.like2)
  
  mean.dev <- -2*sum (colMeans(log.like)) 
  return(mean.dev) 
}
md.jm.a <- dicf(res_jm,which(rownames(sum.df)=="ll.a[1]"),which(rownames(sum.df)=="ll.a[1415]")) 
print(md.jm.a)

md.jm.e <- dicf(res_jm,which(rownames(sum.df)=="ll.e[1]"),which(rownames(sum.df)=="ll.e[1415]")) 
print(md.jm.e)


## calculate DIC = mean D + pd = 2*mean D - D mean
sum.df <- as.data.frame(sum[,2:12])
cp1 <- sum.df[4,4] 
cp2 <- sum.df[5,4] 
c0 <- sum.df[6,4]
c1 <- sum.df[7,4]
c2 <- sum.df[8,4]
c3 <- sum.df[9,4] 
c4 <- sum.df[10,4] 
c5 <- sum.df[11,4] 
c6 <- sum.df[12,4] 

b0 <- sum.df[17,4]
b1 <- sum.df[18,4]
b2 <- sum.df[19,4]
b3 <- sum.df[20,4] 
b <- c(b1,b2,b3)
a <- sum.df[21,4]
ga <- sum.df[22,4]
ga1 <- sum.df[23,4]
ga2 <- sum.df[24,4]
ga3 <- sum.df[25,4] 

u <- sum.df[32:1446,4]
u1 <- sum.df[1447:2861,4]
u2 <- sum.df[2862:4276,4]
u3 <- sum.df[4277:5691,4]
v <- sum.df[5692:7106,4]
w <- sum.df[7107:8521,4]

p2 <- matrix(NA, M, 116, byrow=TRUE)
I1<-matrix(NA, M, 116, byrow=TRUE)
I2<-matrix(NA, M, 116, byrow=TRUE)
lambda0 <- matrix(NA, M, 42, byrow=TRUE)
lambda <- matrix(NA, M, 42, byrow=TRUE)
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


dic.a <- md.jm.a*2-dev_hat.a 
dic.a
dic.e <- md.jm.e*2-dev_hat.e 
dic.e

total.dic <- dic.a+dic.e
total.dic
total.waic <- waic.a+waic.e
total.waic

cbind(round(dic.a),round(waic.a),round(dic.e),round(waic.e), round(total.dic),round(total.waic))
#      [,1]  [,2] [,3] [,4]  [,5]  [,6]
#[1,] 25359 25245 9889 9876 35248 35121
