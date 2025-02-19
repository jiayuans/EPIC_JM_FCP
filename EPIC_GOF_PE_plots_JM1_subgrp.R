########################################################################
#####################################################################
#####################################################################
## Name: Check GOF of JM1 for PE - Population level
## Author: Jiayuan Shi
## Date: Oct 2024
##
## Output: GOF_PE_plots_JM1_population_bysubgrp_comb_n1000.png
#####################################################################
#####################################################################
##################################################################### 

library(tidyverse)
library(survival)

#################################################################
## --------------------------------------------------------------
####### PA dataset
## --------------------------------------------------------------
dirg <- "Z:/EJCStudents/ShiJ/EPIC-CF/Data/Source/"
dat_pa <- read.csv(file=paste(dirg,"/Data-PA-cohort.csv",sep=""))
dat_pa <- dat_pa[with(dat_pa, order(cffidno, VisitAge)),]
head(dat_pa)
summary(dat_pa)
nrow(dat_pa)

length(unique(dat_pa$cffidno)) # 1734 individuals
X.dat.pa <- dat_pa[!duplicated(dat_pa$cffidno,dat_pa$sexf,dat_pa$mut), c(2,6,8)]

#################################################################
## --------------------------------------------------------------
####### PE dataset
## --------------------------------------------------------------
dat_pe0 <- read.csv(file=paste(dirg,"/Data-multiple-final-cohort-EPICstart-Bt-LengthPEx.csv",sep=""))
head(dat_pe0)
nrow(dat_pe0) 
#M<- length(unique(dat_pe0$cffidno)) # 1734
dat_pe0<-dat_pe0[,c(1:3,5,6,7)]

dat_pe <- merge(X.dat.pa,dat_pe0,by="cffidno")

#################################################################
## --------------------------------------------------------------
####### Posterior results
## --------------------------------------------------------------
result <- read.csv(file="Z:/EJCStudents/ShiJ/EPIC-CF/Result/Joint_PA_PEX_model_addrandom_rec_model2_covariate_lt.csv")
c0 <- round(result[4,5],2)
c1 <- round(result[5,5],2)
c2 <- round(result[6,5],2)
c3 <- round(result[7,5],2) 
c4 <- round(result[8,5],2) 
c5 <- round(result[9,5],2) 
c6 <- round(result[10,5],2) 
cp1 <- round(result[11,5],2) 
cp2 <- round(result[12,5],2) 
u <- round(result[13:1746,5],2)
u1 <- round(result[1747:3480,5],2)
u2 <- round(result[3481:5214,5],2)
u3 <- round(result[5215:6948,5],2)

u.tau.inv <- round(result[6950,5],2)
u1.tau.inv <- round(result[6952,5],2)
u2.tau.inv <- round(result[6954,5],2)
u3.tau.inv <- round(result[6956,5],3)

c=c(c1,c2,c3,c4,c5,c6)
cp=c(cp1,cp2)

b0 <- round(result[6960,5],2)
b1 <- round(result[6961,5],2)
b2 <- round(result[6962,5],2)
b3 <- round(result[6963,5],2) 
b <- c(b1,b2,b3)
a <- round(result[6964,5],2)
v <- round(result[6965:8698,5],2)

ga <- round(result[8699,5],2)
ga1 <- round(result[8700,5],2)
ga2 <- round(result[8701,5],2)
ga3 <- round(result[8702,5],2) 

w <- round(result[8703:10436,5],2) 
w.tau.inv <- round(result[10438,5],2)

#################################################################
## --------------------------------------------------------------
####### PE plots - observed
## --------------------------------------------------------------
####### Plotting the NElson-Aalen estimator
# N-A Cumulative Number of Events plot
cox <- coxph(Surv(tstart,tend,status) ~ 1, data=dat_pe)
na.fit <- survfit(cox,type="aalen")

h.aalen <- (-log(c(1,na.fit$surv)))
aalen.est <- cbind(time=c(0,na.fit$time),h.aalen)
aalen.est <- as.data.frame(aalen.est)
plot(aalen.est$time, aalen.est$h.aalen, type="s", xlab="Age (Years)",
     ylab="Cumulative Number of Events")
mtext("JM1\n")


# N-A Cumulative Number of Events plot for females and males
cox.f <- coxph(Surv(tstart,tend,status) ~ 1, data=dat_pe, subset=(sexf==1))
na.fit.f <- survfit(cox.f,type="aalen")
h.aalen.f <- (-log(c(1,na.fit.f$surv)))
aalen.est.f <- cbind(time=c(0,na.fit.f$time), h.aalen.f)
aalen.est.f <- as.data.frame(aalen.est.f)

cox.m <- coxph(Surv(tstart,tend,status) ~ 1, data=dat_pe, subset=(sexf==0))
na.fit.m <- survfit(cox.m,type="aalen")
h.aalen.m <- (-log(c(1,na.fit.m$surv)))
aalen.est.m <- cbind(time=c(0,na.fit.m$time), h.aalen.m)
aalen.est.m <- as.data.frame(aalen.est.m)

plot(aalen.est.f$time, aalen.est.f$h.aalen, col="indianred4", type="s", xlab="Age (Years)",
     ylab="Cumulative Number of Events")
mtext("JM1 by Gender\n")
lines(aalen.est.m$time, aalen.est.m$h.aalen,type="s",col="royalblue2")

legend(0,9, c("Female","Male"), lwd=c(2,2),col=c("indianred4","royalblue2"))


# N-A Cumulative Number of Events plot for mutation types
cox.2alleles <- coxph(Surv(tstart,tend,status) ~ 1, data=dat_pe, subset=(mut==1))
na.fit.2alleles <- survfit(cox.2alleles,type="aalen")
h.aalen.2alleles <- (-log(c(1,na.fit.2alleles$surv)))
aalen.est.2alleles <- cbind(time=c(0,na.fit.2alleles$time), h.aalen.2alleles)
aalen.est.2alleles <- as.data.frame(aalen.est.2alleles)

cox.1alleles  <- coxph(Surv(tstart,tend,status) ~ 1, data=dat_pe, subset=(mut==2))
na.fit.1alleles <- survfit(cox.1alleles,type="aalen")
h.aalen.1alleles <- (-log(c(1,na.fit.1alleles$surv)))
aalen.est.1alleles <- cbind(time=c(0,na.fit.1alleles$time), h.aalen.1alleles)
aalen.est.1alleles <- as.data.frame(aalen.est.1alleles)

cox.0alleles  <- coxph(Surv(tstart,tend,status) ~ 1, data=dat_pe, subset=(mut==3))
na.fit.0alleles <- survfit(cox.0alleles,type="aalen")
h.aalen.0alleles <- (-log(c(1,na.fit.0alleles$surv)))
aalen.est.0alleles <- cbind(time=c(0,na.fit.0alleles$time), h.aalen.0alleles)
aalen.est.0alleles <- as.data.frame(aalen.est.0alleles)

plot(aalen.est.2alleles$time, aalen.est.2alleles$h.aalen, col="indianred4", type="s", xlab="Age (Years)",
     ylab="Cumulative Number of Events")
mtext("JM1 for F508del groups\n")
lines(aalen.est.1alleles$time, aalen.est.1alleles$h.aalen,type="s",col="royalblue2")
lines(aalen.est.0alleles$time, aalen.est.0alleles$h.aalen,type="s",col="lightgreen")

legend(0,9, c("Two alleles F508del","One allele F508del","Other or unknown F508del"), lwd=c(2,2,2),col=c("indianred4","royalblue2","lightgreen"))



#################################################################
## --------------------------------------------------------------
####### PE plots - simulated
## --------------------------------------------------------------
#################Reading in data for covariates#############################
X.dat_pe <- dat_pe[!duplicated(dat_pe$cffidno,dat_pe$sexf,dat_pe$mut), ]
X.dat_pe1 <- model.matrix(status ~ sexf + as.factor(mut), data = X.dat_pe)
X1_pe <- as.numeric(X.dat_pe1[,2]) ## sexf: female
X2_pe <- as.numeric(X.dat_pe1[,3]) ## mut2: one allele F508del
X3_pe <- as.numeric(X.dat_pe1[,4]) ## mut3: other or unknown F508del

## extract time of first visit and last visit for each individual
dat_pe1 <-aggregate(dat_pe$tend, by=list(dat_pe$cffidno),
                    FUN=max, na.rm=TRUE)
names(dat_pe1) <- c('cffidno','age.max')

dat_pe2 <-aggregate(dat_pe$tstart, by=list(dat_pe$cffidno),
                    FUN=min, na.rm=TRUE)
names(dat_pe2) <- c('cffidno','age.min')

t<-dat_pe2$age.min
tt<-dat_pe1$age.max

N<-length(tt)

#########################################################################
# Functions that create an event times dataset for a poisson process (continuous data)
# -------------- Simulating the event times -----
set.seed(123)
NHPP<-function(a,b,T){
  mu <- b*T^a  # Mean of the Poisson process up to time T
  n <-rpois(1, mu)  #  number of events poisson
  if (n!=0) {
    u <- runif(n,0,1) # n uniforms
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

# -------------- Building the simulated poisson data -----
poisson.d <- function(alpha,beta1,beta2,beta3,beta0,x1,x2,x3,TTei){
  c_0i <- rnorm(N,0,sqrt(w.tau.inv))
  b_0i <- rnorm(N,0,sqrt(u.tau.inv))
  b_1i <- rnorm(N,0,sqrt(u1.tau.inv))
  b_2i <- rnorm(N,0,sqrt(u2.tau.inv))
  b_3i <- rnorm(N,0,sqrt(u3.tau.inv))
  vi <- exp(ga*b_0i+c_0i+ga1*b_1i+ga2*b_2i+ga3*b_3i)
  times <- NHPP(b=vi[1]*exp(beta1*x1[1]+beta2*x2[1]+beta3*x3[1])*exp(beta0),a=alpha,T=TTei[1])
  start <-  times[,1]
  stop <- times[,2]
  status <- times[,3]
  n.rec <- times[,4]
  id <- rep(1,length(stop))
  x1i <- rep(x1[1],length(stop))
  x2i <- rep(x2[1],length(stop))
  x3i <- rep(x3[1],length(stop))
  Tei <- rep(TTei[1],length(stop))
  for (i in 2:length(x1)){
    times2 <- NHPP(b=vi[i]*exp(beta0+beta1*x1[i]+beta2*x2[i]+beta3*x3[i]),a=alpha,T=TTei[i]) 
    start2 <-  times2[,1]
    stop2 <- times2[,2]
    status2 <- times2[,3]
    n.rec2 <- times2[,4]
    id <- c(id,rep(i,length(stop2)))
    x1i <- c(x1i,rep(x1[i],length(stop2)))
    x2i <- c(x2i,rep(x2[i],length(stop2)))
    x3i <- c(x3i,rep(x3[i],length(stop2)))
    Tei <- c(Tei,rep(TTei[i],length(stop2)))
    
    start <- c(start,start2)
    stop <- c(stop,stop2)
    status <- c(status,status2)
    n.rec <- c(n.rec,n.rec2)
  }
  return(data.frame(id,x1i,x2i,x3i,Tei,n.rec,start,stop,status))
}

simdat.pe_list <- list()
I=1000
num <- rep(NA,I)

for (r in 1:I){
## simulated recurrent event data  
simdat.pe00 <- poisson.d(alpha=a,beta1=b1,beta2=b2,beta3=b3,beta0=b0,x1=X1_pe,x2=X2_pe,x3=X3_pe,TTei=tt)

id <- c(1:N)
timeS <- as.data.frame(cbind(id,t)) ## left truncation time
timeE <- as.data.frame(cbind(id,tt))

## generate the simulated data with left truncation
simdat.pe0 <- merge(simdat.pe00, timeS,all=TRUE)
simdat.pe <- subset(simdat.pe0, stop >= t) ## not include those event times before the truncation time
simdat.pe1 <- simdat.pe
simdat.pe1$start <- ifelse(simdat.pe1$start < simdat.pe1$t, simdat.pe1$t, simdat.pe1$start) ## for those start times before the truncation time, set them to be the the truncation time
simdat.pe_list[[r]] <- simdat.pe1
num[r] <- length(unique(simdat.pe1$id)) ## still have 1734 individuals
} 

png(file = "Z:/EJCStudents/ShiJ/EPIC-CF/Result/figure_PE_GOF/GOF_PE_plots_JM1_population_bysubgrp_comb_n1000.png", width = 700, height = 700) 

par(mfrow=c(3,2))

## females and two alleles F508del
cox.f.2 <- coxph(Surv(tstart,tend,status) ~ 1, data=dat_pe, subset=(sexf==1 & mut==1))
na.fit.f.2 <- survfit(cox.f.2,type="aalen")
h.aalen.f.2 <- (-log(c(1,na.fit.f.2$surv)))
aalen.est.f.2 <- cbind(time=c(0,na.fit.f.2$time), h.aalen.f.2)
aalen.est.f.2 <- as.data.frame(aalen.est.f.2)

plot(aalen.est.f.2$time, aalen.est.f.2$h.aalen, type="s", xlab="Age (Years)",
     ylab="Cumulative Number of Events")
mtext("JM1 for females and two alleles F508del")


for (j in 1:length(simdat.pe_list)) {
  cox <- coxph(Surv(start,stop,status) ~ 1, data=simdat.pe_list[[j]], subset=(x1i==1 & x2i==0 & x3i==0), timefix = FALSE )
  na.fit <- survfit(cox,type="aalen")
  h.aalen <- (-log(c(1,na.fit$surv)))
  aalen.est <- cbind(time=c(0,na.fit$time),h.aalen)
  aalen.est <- as.data.frame(aalen.est)
  lines(aalen.est$time, aalen.est$h.aalen, type="s", col="pink",lwd = 3)
}
lines(aalen.est.f.2$time, aalen.est.f.2$h.aalen, type="s")

legend("topleft", c("observed","simulated"), lty=c(1,1), col = c("black","black","pink"))

## females and one allele F508del
cox.f.1 <- coxph(Surv(tstart,tend,status) ~ 1, data=dat_pe, subset=(sexf==1 & mut==2))
na.fit.f.1 <- survfit(cox.f.1,type="aalen")
h.aalen.f.1 <- (-log(c(1,na.fit.f.1$surv)))
aalen.est.f.1 <- cbind(time=c(0,na.fit.f.1$time), h.aalen.f.1)
aalen.est.f.1 <- as.data.frame(aalen.est.f.1)

plot(aalen.est.f.1$time, aalen.est.f.1$h.aalen, type="s", xlab="Age (Years)",
     ylab="Cumulative Number of Events")
mtext("JM1 for females and one allele F508del")


for (j in 1:length(simdat.pe_list)) {
  cox <- coxph(Surv(start,stop,status) ~ 1, data=simdat.pe_list[[j]], subset=(x1i==1 & x2i==1 & x3i==0),timefix = FALSE  )
  na.fit <- survfit(cox,type="aalen")
  h.aalen <- (-log(c(1,na.fit$surv)))
  aalen.est <- cbind(time=c(0,na.fit$time),h.aalen)
  aalen.est <- as.data.frame(aalen.est)
  lines(aalen.est$time, aalen.est$h.aalen, type="s", col="pink",lwd = 3)
}
lines(aalen.est.f.1$time, aalen.est.f.1$h.aalen, type="s")

legend("topleft", c("observed","simulated"), lty=c(1,1), col = c("black","black","pink"))


## females and other or unknown F508del
cox.f.0 <- coxph(Surv(tstart,tend,status) ~ 1, data=dat_pe, subset=(sexf==1 & mut==3),timefix = FALSE )
na.fit.f.0 <- survfit(cox.f.0,type="aalen")
h.aalen.f.0 <- (-log(c(1,na.fit.f.0$surv)))
aalen.est.f.0 <- cbind(time=c(0,na.fit.f.0$time), h.aalen.f.0)
aalen.est.f.0 <- as.data.frame(aalen.est.f.0)

plot(aalen.est.f.0$time, aalen.est.f.0$h.aalen, type="s", xlab="Age (Years)",
     ylab="Cumulative Number of Events")
mtext("JM1 for females and other or unknown F508del")


for (j in 1:length(simdat.pe_list)) {
  cox <- coxph(Surv(start,stop,status) ~ 1, data=simdat.pe_list[[j]], subset=(x1i==1 & x2i==0 & x3i==1),timefix = FALSE )
  na.fit <- survfit(cox,type="aalen")
  h.aalen <- (-log(c(1,na.fit$surv)))
  aalen.est <- cbind(time=c(0,na.fit$time),h.aalen)
  aalen.est <- as.data.frame(aalen.est)
  lines(aalen.est$time, aalen.est$h.aalen, type="s", col="pink",lwd = 3)
}
lines(aalen.est.f.0$time, aalen.est.f.0$h.aalen, type="s")

legend("topleft", c("observed","simulated"), lty=c(1,1), col = c("black","black","pink"))


## males and two alleles F508del
cox.m.0 <- coxph(Surv(tstart,tend,status) ~ 1, data=dat_pe, subset=(sexf==0 & mut==1))
na.fit.m.0 <- survfit(cox.m.0,type="aalen")
h.aalen.m.0 <- (-log(c(1,na.fit.m.0$surv)))
aalen.est.m.0 <- cbind(time=c(0,na.fit.m.0$time), h.aalen.m.0)
aalen.est.m.0 <- as.data.frame(aalen.est.m.0)

plot(aalen.est.m.0$time, aalen.est.m.0$h.aalen, type="s", xlab="Age (Years)",
     ylab="Cumulative Number of Events")
mtext("JM1 for males and two alleles F508del")

for (j in 1:length(simdat.pe_list)) {
  cox <- coxph(Surv(start,stop,status) ~ 1, data=simdat.pe_list[[j]], subset=(x1i==0 & x2i==0 & x3i==0) ,timefix = FALSE )
  na.fit <- survfit(cox,type="aalen")
  h.aalen <- (-log(c(1,na.fit$surv)))
  aalen.est <- cbind(time=c(0,na.fit$time),h.aalen)
  aalen.est <- as.data.frame(aalen.est)
  lines(aalen.est$time, aalen.est$h.aalen, type="s", col="pink",lwd = 3)
}
lines(aalen.est.m.0$time, aalen.est.m.0$h.aalen, type="s")

legend("topleft", c("observed","simulated"), lty=c(1,1), col = c("black","pink"))



## males and one allele F508del
cox.m.0 <- coxph(Surv(tstart,tend,status) ~ 1, data=dat_pe, subset=(sexf==0 & mut==2))
na.fit.m.0 <- survfit(cox.m.0,type="aalen")
h.aalen.m.0 <- (-log(c(1,na.fit.m.0$surv)))
aalen.est.m.0 <- cbind(time=c(0,na.fit.m.0$time), h.aalen.m.0)
aalen.est.m.0 <- as.data.frame(aalen.est.m.0)

plot(aalen.est.m.0$time, aalen.est.m.0$h.aalen, type="s", xlab="Age (Years)",
     ylab="Cumulative Number of Events")
mtext("JM1 for males and one allele F508del")

for (j in 1:length(simdat.pe_list)) {
  cox <- coxph(Surv(start,stop,status) ~ 1, data=simdat.pe_list[[j]], subset=(x1i==0 & x2i==1 & x3i==0),timefix = FALSE  )
  na.fit <- survfit(cox,type="aalen")
  h.aalen <- (-log(c(1,na.fit$surv)))
  aalen.est <- cbind(time=c(0,na.fit$time),h.aalen)
  aalen.est <- as.data.frame(aalen.est)
  lines(aalen.est$time, aalen.est$h.aalen, type="s", col="pink",lwd = 3)
}
lines(aalen.est.m.0$time, aalen.est.m.0$h.aalen, type="s")

legend("topleft", c("observed","simulated"), lty=c(1,1), col = c("black","pink"))


## males and other or unknown F508del
cox.m.0 <- coxph(Surv(tstart,tend,status) ~ 1, data=dat_pe, subset=(sexf==0 & mut==3))
na.fit.m.0 <- survfit(cox.m.0,type="aalen")
h.aalen.m.0 <- (-log(c(1,na.fit.m.0$surv)))
aalen.est.m.0 <- cbind(time=c(0,na.fit.m.0$time), h.aalen.m.0)
aalen.est.m.0 <- as.data.frame(aalen.est.m.0)

plot(aalen.est.m.0$time, aalen.est.m.0$h.aalen, type="s", xlab="Age (Years)",
     ylab="Cumulative Number of Events")
mtext("JM1 for males and other or unknown F508del")

for (j in 1:length(simdat.pe_list)) {
  cox <- coxph(Surv(start,stop,status) ~ 1, data=simdat.pe_list[[j]], subset=(x1i==0 & x2i==0 & x3i==1),timefix = FALSE  )
  na.fit <- survfit(cox,type="aalen")
  h.aalen <- (-log(c(1,na.fit$surv)))
  aalen.est <- cbind(time=c(0,na.fit$time),h.aalen)
  aalen.est <- as.data.frame(aalen.est)
  lines(aalen.est$time, aalen.est$h.aalen, type="s", col="pink",lwd = 3)
}
lines(aalen.est.m.0$time, aalen.est.m.0$h.aalen, type="s")

legend("topleft", c("observed","simulated"), lty=c(1,1), col = c("black","pink"))



dev.off()
