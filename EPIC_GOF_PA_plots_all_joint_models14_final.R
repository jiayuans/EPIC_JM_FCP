########################################################################
#####################################################################
#####################################################################
## Name: Compare GOF of all JMs (JM1-JM4) for PA - Overall
## Author: Jiayuan Shi
## Date: Oct 2024
##
## Output: GOF_PA_plots_overall_JMs14_n1000.png
#####################################################################
#####################################################################
##################################################################### 

library(tidyverse)
library(survival)
setwd("C:/UCHealth/RA/Project/EPIC-CF/Analysis_Jiayuan")
getwd()

#################################################################
## --------------------------------------------------------------
####### PA dataset
## --------------------------------------------------------------
dat_pa <- read.csv( file="C:/UCHealth/RA/Project/EPIC-CF/Data/Source/Data-PA-cohort.csv")
dat_pa <- dat_pa[with(dat_pa, order(cffidno, VisitAge)),]
head(dat_pa)
summary(dat_pa)
nrow(dat_pa)

length(unique(dat_pa$cffidno)) # 1734 individuals

#################################################################
## --------------------------------------------------------------
####### percentages are calculated by # PA in an interval / total visits in that interval
## --------------------------------------------------------------

#################################################################
## --------------------------------------------------------------
####### PA plots - observed
## --------------------------------------------------------------
dat_pa1 <- dat_pa[with(dat_pa, order(cffidno, VisitAge)),2:5]
#dat_pa1$age.bin <- ceiling(dat_pa1$VisitAge) ## one bin per year
#dat_pa1$age_yr <- dat_pa1$age.bin-0.5  ## mid point of the year interval
dat_pa1$age.bin <- ceiling(dat_pa1$VisitAge*12/2) ## one bin per 2 month
dat_pa1$age.mid <- dat_pa1$age.bin*2-1  ## mid point of the month interval
dat_pa1$age_yr <- dat_pa1$age.mid/12 ## from month to year

## year
dat_pa.n <-aggregate(dat_pa1$cltpa, by=list(dat_pa1$age_yr),
                     FUN=sum, na.rm=TRUE)
names(dat_pa.n) <- c('age_yr','n_pa')

dat_pa.tot <- dat_pa1 %>%
  group_by(age_yr) %>%
  summarise(count = n())
names(dat_pa.tot) <- c('age_yr','n_visit')

dat_pa.cal0 <- merge(dat_pa.tot, dat_pa.n, by=c('age_yr'))
dat_pa.cal <- subset(dat_pa.cal0, age_yr<21)
dat_pa.cal$perc <- dat_pa.cal$n_pa/dat_pa.cal$n_visit


#################################################################
## --------------------------------------------------------------
####### PA plots - simulated - JM1
## --------------------------------------------------------------

#################################################################
## --------------------------------------------------------------
####### Posterior results
## --------------------------------------------------------------
result <- read.csv(file="C:/UCHealth/RA/Project/EPIC-CF/Analysis_Jiayuan/Result/Joint_PA_PEX_model_addrandom_rec_model2_covariate_lt.csv")
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


dat_pa2 <- dat_pa[,-1]
head(dat_pa2, n=10)

M <- length(unique(dat_pa2$cffidno))
count <- dat_pa2 %>% count(cffidno)
max.count <- max(count$n)

##########################Assigning unique number to each subject##########################
dat_pa3 <- dat_pa2 %>% group_by(cffidno) %>% mutate(time = c(1:length(cffidno)))

Yd.temp <- data.frame(cffidno = rep(unique(dat_pa2$cffidno),each=max.count), time = 1:max.count)
head(Yd.temp,n=12)

Y.epic <- merge(dat_pa3,Yd.temp,by=c('cffidno','time'),all.y=TRUE)

#################Readingin data for X vectors#############################
X.dat.pa <- dat_pa2[!duplicated(dat_pa2$cffidno,dat_pa2$sexf,dat_pa2$mut), ]
X.dat.pa1 <- model.matrix(cltpa ~ sexf + as.factor(mut), data = X.dat.pa)
X1_pa <- as.numeric(X.dat.pa1[,2]) ## sexf: female
X2_pa <- as.numeric(X.dat.pa1[,3]) ## mut2: 1 alleles D F508
X3_pa <- as.numeric(X.dat.pa1[,4]) ## mut3: 0 alleles D F508

#################Readingin data for VisitAge matrix#############################
X <- matrix(Y.epic$VisitAge, M, max.count, byrow=TRUE)
X[1,]##need to remove first row of column headers
dim(X)
X[1,]

#################input variables for simulation#####################
#### checking for how many individuals we have NAs in the middle of followup
sum.na <- rep(NA,M)
k=rep(NA,M)

ids <- unique(Y.epic$cffidno) ## 103104 103125 103129 103145 103147
for (i in 1:M){
  na.indices <- which(Y.epic$VisitAge[Y.epic$cffidno==ids[i]] %in% NA)
  if (length(na.indices)==0){
    k[i] <- max.count} else{
      k[i] <- min(na.indices)-1}
}
kk=max(k)

###set number of iterations#################################
I=1000

Y_list <- list()
set.seed(123)
#######################################################
for (r in 1:I){
  u<-rnorm(M,0,sqrt(u.tau.inv))
  u1<-rnorm(M,0,sqrt(u1.tau.inv))
  u2<-rnorm(M,0,sqrt(u2.tau.inv))
  u3<-rnorm(M,0,sqrt(u3.tau.inv))
  
  I1<-matrix(NA, nrow=M, ncol=kk, byrow=TRUE)
  I2<-matrix(NA, nrow=M, ncol=kk, byrow=TRUE)
  p2<-matrix(NA, nrow=M, ncol=kk, byrow=TRUE)
  Y<-matrix(NA, nrow=M, ncol=kk, byrow=TRUE)
  
  for (i in 1:M){
    for (j in 1:k[i]){
      I1[i,j]<-ifelse(X[i,j]< cp[1],-1,1)
      I2[i,j]<-ifelse(X[i,j]< cp[2],-1,1)
      p2[i,j]<-exp(c0+(c[1]+u1[i])*(X[i,j]-cp[1])+(c[2]+u2[i])*(X[i,j]-cp[1])*I1[i,j]+(c[3]+u3[i])*(X[i,j]-cp[2])*I2[i,j]+c[4]*X1_pa[i]+c[5]*X2_pa[i]+c[6]*X3_pa[i]+u[i])/(1+exp(c0+(c[1]+u1[i])*(X[i,j]-cp[1])+(c[2]+u2[i])*(X[i,j]-cp[1])*I1[i,j]+(c[3]+u3[i])*(X[i,j]-cp[2])*I2[i,j]+c[4]*X1_pa[i]+c[5]*X2_pa[i]+c[6]*X3_pa[i]+u[i]))
      Y[i,j]=rbinom(1, 1, p2[i,j])
    }
  }
  Y_list[[r]] <- Y
} 

print(Y_list[[1]])

###################################################################
dat_pa.cal.sim_list <- list()

for (i in 1:length(Y_list)) {
  Y <- Y_list[[i]]
  X_df <- as.data.frame(X)
  Y_df <- as.data.frame(Y)
  X_df$cffidno <- unique(dat_pa$cffidno)
  Y_df$cffidno <- unique(dat_pa$cffidno)
  
  X_df_long = X_df %>%
    pivot_longer(V1:V126, names_to = "visit", values_to = "VisitAge") 
  
  Y_df_long = Y_df %>%
    pivot_longer(V1:V126, names_to = "visit", values_to = "cltpa.sim") 
  
  dat.sim <- cbind(X_df_long,Y_df_long[,3])
  dat.sim1 <- subset(dat.sim, is.na(VisitAge) != T)
  dat.sim2 <- dat.sim1[with(dat.sim1, order(cffidno, VisitAge)),]
  #dat.sim2$age.bin <- ifelse(dat.sim2$VisitAge!=0, ceiling(dat.sim2$VisitAge),1)
  #dat.sim2$age_yr <- ifelse(dat.sim2$VisitAge!=0, dat.sim2$age.bin-0.5,0)
  dat.sim2$age.bin <- ifelse(dat.sim2$VisitAge!=0, ceiling(dat.sim2$VisitAge*12/2),1)
  dat.sim2$age <- ifelse(dat.sim2$VisitAge!=0, dat.sim2$age.bin*2-1,0)
  dat.sim2$age_yr <- dat.sim2$age/12
  
  ## year
  dat_pa.n.sim <-aggregate(dat.sim2$cltpa.sim, by=list(dat.sim2$age_yr),
                           FUN=sum, na.rm=TRUE)
  names(dat_pa.n.sim) <- c('age_yr','n_pa')
  
  dat_pa.tot.sim <- dat.sim2 %>%
    group_by(age_yr) %>%
    summarise(count = n())
  names(dat_pa.tot.sim) <- c('age_yr','n_visit')
  
  dat_pa.cal.sim0 <- merge(dat_pa.tot.sim, dat_pa.n.sim, by=c('age_yr'))
  dat_pa.cal.sim <- subset(dat_pa.cal.sim0, age_yr<21)
  dat_pa.cal.sim$perc <- dat_pa.cal.sim$n_pa/dat_pa.cal.sim$n_visit
  
  dat_pa.cal.sim_list[[i]] <- dat_pa.cal.sim
}



#################################################################
## --------------------------------------------------------------
####### PA plots - observed + simulated
## --------------------------------------------------------------

png(file = "C:/UCHealth/RA/Project/EPIC-CF/Analysis_Jiayuan/Result/figure_PA_GOF/GOF_PA_plots_overall_JMs14_n1000.png", width = 700, height = 700) 
par(mfrow=c(2,2))
plot(dat_pa.cal$age_yr,dat_pa.cal$perc, type = "l", xlab = "Age (Year)", ylab = "Proportion of PA", main="JM1", ylim = c(0,1)) ## renamed the JM, this is the one with random slope + intercept
for (i in 1:length(Y_list)) {
  lines(dat_pa.cal.sim_list[[i]]$age_yr,dat_pa.cal.sim_list[[i]]$perc, type = "l",col="pink",lwd = 3)
}
lines(dat_pa.cal$age_yr,dat_pa.cal$perc, type = "l")

legend(0,1, c("observed","simulated"), lty=c(1,1),col=c("black","pink"))





#################################################################
## --------------------------------------------------------------
####### PA plots - simulated - JM2
## --------------------------------------------------------------

#################################################################
## --------------------------------------------------------------
####### Posterior results
## --------------------------------------------------------------
result <- read.csv(file="C:/UCHealth/RA/Project/EPIC-CF/Analysis_Jiayuan/Result/Joint_PA_PEX_model_addrandom_rec_model1_covariate_lt.csv")
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

u.tau.inv <- round(result[1747,5],2)

c=c(c1,c2,c3,c4,c5,c6)
cp=c(cp1,cp2)


dat_pa2 <- dat_pa[,-1]
head(dat_pa2, n=10)

M <- length(unique(dat_pa2$cffidno))
count <- dat_pa2 %>% count(cffidno)
max.count <- max(count$n)

##########################Assigning unique number to each subject##########################
dat_pa3 <- dat_pa2 %>% group_by(cffidno) %>% mutate(time = c(1:length(cffidno)))

Yd.temp <- data.frame(cffidno = rep(unique(dat_pa2$cffidno),each=max.count), time = 1:max.count)
head(Yd.temp,n=12)

Y.epic <- merge(dat_pa3,Yd.temp,by=c('cffidno','time'),all.y=TRUE)

#################Readingin data for X vectors#############################
X.dat.pa <- dat_pa2[!duplicated(dat_pa2$cffidno,dat_pa2$sexf,dat_pa2$mut), ]
X.dat.pa1 <- model.matrix(cltpa ~ sexf + as.factor(mut), data = X.dat.pa)
X1_pa <- as.numeric(X.dat.pa1[,2]) ## sexf: female
X2_pa <- as.numeric(X.dat.pa1[,3]) ## mut2: 1 alleles D F508
X3_pa <- as.numeric(X.dat.pa1[,4]) ## mut3: 0 alleles D F508

#################Readingin data for VisitAge matrix#############################
X <- matrix(Y.epic$VisitAge, M, max.count, byrow=TRUE)
X[1,]##need to remove first row of column headers
dim(X)
X[1,]

#################input variables for simulation#####################
#### checking for how many individuals we have NAs in the middle of followup
sum.na <- rep(NA,M)
k=rep(NA,M)

ids <- unique(Y.epic$cffidno) ## 103104 103125 103129 103145 103147
for (i in 1:M){
  na.indices <- which(Y.epic$VisitAge[Y.epic$cffidno==ids[i]] %in% NA)
  if (length(na.indices)==0){
    k[i] <- max.count} else{
      k[i] <- min(na.indices)-1}
}
kk=max(k)

###set number of iterations#################################
I=1000

Y_list <- list()
set.seed(123)
#######################################################
for (r in 1:I){
  u<-rnorm(M,0,sqrt(u.tau.inv))
  
  I1<-matrix(NA, nrow=M, ncol=kk, byrow=TRUE)
  I2<-matrix(NA, nrow=M, ncol=kk, byrow=TRUE)
  p2<-matrix(NA, nrow=M, ncol=kk, byrow=TRUE)
  Y<-matrix(NA, nrow=M, ncol=kk, byrow=TRUE)
  
  for (i in 1:M){
    for (j in 1:k[i]){
      I1[i,j]<-ifelse(X[i,j]< cp[1],-1,1)
      I2[i,j]<-ifelse(X[i,j]< cp[2],-1,1)
      p2[i,j]<-exp(c0+c[1]*(X[i,j]-cp[1])+c[2]*(X[i,j]-cp[1])*I1[i,j]+c[3]*(X[i,j]-cp[2])*I2[i,j]+c[4]*X1_pa[i]+c[5]*X2_pa[i]+c[6]*X3_pa[i]+u[i])/(1+exp(c0+c[1]*(X[i,j]-cp[1])+c[2]*(X[i,j]-cp[1])*I1[i,j]+c[3]*(X[i,j]-cp[2])*I2[i,j]+c[4]*X1_pa[i]+c[5]*X2_pa[i]+c[6]*X3_pa[i]+u[i]))
      Y[i,j]=rbinom(1, 1, p2[i,j])
    }
  }
  Y_list[[r]] <- Y
} 

print(Y_list[[1]])

###################################################################
dat_pa.cal.sim_list <- list()

for (i in 1:length(Y_list)) {
  Y <- Y_list[[i]]
  X_df <- as.data.frame(X)
  Y_df <- as.data.frame(Y)
  X_df$cffidno <- unique(dat_pa$cffidno)
  Y_df$cffidno <- unique(dat_pa$cffidno)
  
  X_df_long = X_df %>%
    pivot_longer(V1:V126, names_to = "visit", values_to = "VisitAge") 
  
  Y_df_long = Y_df %>%
    pivot_longer(V1:V126, names_to = "visit", values_to = "cltpa.sim") 
  
  dat.sim <- cbind(X_df_long,Y_df_long[,3])
  dat.sim1 <- subset(dat.sim, is.na(VisitAge) != T)
  dat.sim2 <- dat.sim1[with(dat.sim1, order(cffidno, VisitAge)),]
  #dat.sim2$age.bin <- ifelse(dat.sim2$VisitAge!=0, ceiling(dat.sim2$VisitAge),1)
  #dat.sim2$age_yr <- ifelse(dat.sim2$VisitAge!=0, dat.sim2$age.bin-0.5,0)
  dat.sim2$age.bin <- ifelse(dat.sim2$VisitAge!=0, ceiling(dat.sim2$VisitAge*12/2),1)
  dat.sim2$age <- ifelse(dat.sim2$VisitAge!=0, dat.sim2$age.bin*2-1,0)
  dat.sim2$age_yr <- dat.sim2$age/12
  
  ## year
  dat_pa.n.sim <-aggregate(dat.sim2$cltpa.sim, by=list(dat.sim2$age_yr),
                           FUN=sum, na.rm=TRUE)
  names(dat_pa.n.sim) <- c('age_yr','n_pa')
  
  dat_pa.tot.sim <- dat.sim2 %>%
    group_by(age_yr) %>%
    summarise(count = n())
  names(dat_pa.tot.sim) <- c('age_yr','n_visit')
  
  dat_pa.cal.sim0 <- merge(dat_pa.tot.sim, dat_pa.n.sim, by=c('age_yr'))
  dat_pa.cal.sim <- subset(dat_pa.cal.sim0, age_yr<21)
  dat_pa.cal.sim$perc <- dat_pa.cal.sim$n_pa/dat_pa.cal.sim$n_visit
  
  dat_pa.cal.sim_list[[i]] <- dat_pa.cal.sim
}



#################################################################
## --------------------------------------------------------------
####### PA plots - observed + simulated
## --------------------------------------------------------------
plot(dat_pa.cal$age_yr,dat_pa.cal$perc, type = "l", xlab = "Age (Year)", ylab = "Proportion of PA", main="JM2", ylim = c(0,1)) 
for (i in 1:length(Y_list)) {
  lines(dat_pa.cal.sim_list[[i]]$age_yr,dat_pa.cal.sim_list[[i]]$perc, type = "l",col="pink",lwd = 3)
}
lines(dat_pa.cal$age_yr,dat_pa.cal$perc, type = "l")

legend(0,1, c("observed","simulated"), lty=c(1,1),col=c("black","pink"))



#################################################################
## --------------------------------------------------------------
####### PA plots - simulated - JM3
## --------------------------------------------------------------

#################################################################
## --------------------------------------------------------------
####### Posterior results
## --------------------------------------------------------------
result <- read.csv(file="C:/UCHealth/RA/Project/EPIC-CF/Analysis_Jiayuan/Result/Joint_PA_PEX_model_addrandom_rec_model0b_covariate_lt.csv")
c0 <- round(result[1,5],2)
c1 <- round(result[2,5],2)
c2 <- round(result[3,5],2)
c3 <- round(result[4,5],2)
c4 <- round(result[5,5],3) 

u <- round(result[6:1739,5],2)
u1 <- round(result[1740:3473,5],2)

u.tau.inv <- round(result[3475,5],2)
u1.tau.inv <- round(result[3477,5],2)

c=c(c1,c2,c3,c4)

dat_pa2 <- dat_pa[,-1]
head(dat_pa2, n=10)

M <- length(unique(dat_pa2$cffidno))
count <- dat_pa2 %>% count(cffidno)
max.count <- max(count$n)

##########################Assigning unique number to each subject##########################
dat_pa3 <- dat_pa2 %>% group_by(cffidno) %>% mutate(time = c(1:length(cffidno)))

Yd.temp <- data.frame(cffidno = rep(unique(dat_pa2$cffidno),each=max.count), time = 1:max.count)
head(Yd.temp,n=12)

Y.epic <- merge(dat_pa3,Yd.temp,by=c('cffidno','time'),all.y=TRUE)

#################Readingin data for X vectors#############################
X.dat.pa <- dat_pa2[!duplicated(dat_pa2$cffidno,dat_pa2$sexf,dat_pa2$mut), ]
X.dat.pa1 <- model.matrix(cltpa ~ sexf + as.factor(mut), data = X.dat.pa)
X1_pa <- as.numeric(X.dat.pa1[,2]) ## sexf: female
X2_pa <- as.numeric(X.dat.pa1[,3]) ## mut2: 1 alleles D F508
X3_pa <- as.numeric(X.dat.pa1[,4]) ## mut3: 0 alleles D F508

#################Readingin data for VisitAge matrix#############################
X <- matrix(Y.epic$VisitAge, M, max.count, byrow=TRUE)
X[1,]##need to remove first row of column headers
dim(X)
X[1,]

#################input variables for simulation#####################
#### checking for how many individuals we have NAs in the middle of followup
sum.na <- rep(NA,M)
k=rep(NA,M)

ids <- unique(Y.epic$cffidno) ## 103104 103125 103129 103145 103147
for (i in 1:M){
  na.indices <- which(Y.epic$VisitAge[Y.epic$cffidno==ids[i]] %in% NA)
  if (length(na.indices)==0){
    k[i] <- max.count} else{
      k[i] <- min(na.indices)-1}
}
kk=max(k)

###set number of iterations#################################
I=1000

Y_list <- list()
set.seed(123)
#######################################################
for (r in 1:I){
  u<-rnorm(M,0,sqrt(u.tau.inv))
  u1<-rnorm(M,0,sqrt(u1.tau.inv))
  
  I1<-matrix(NA, nrow=M, ncol=kk, byrow=TRUE)
  I2<-matrix(NA, nrow=M, ncol=kk, byrow=TRUE)
  p2<-matrix(NA, nrow=M, ncol=kk, byrow=TRUE)
  Y<-matrix(NA, nrow=M, ncol=kk, byrow=TRUE)
  
  for (i in 1:M){
    for (j in 1:k[i]){
      I1[i,j]<-ifelse(X[i,j]< cp[1],-1,1)
      I2[i,j]<-ifelse(X[i,j]< cp[2],-1,1)
      p2[i,j]=exp(c0+(c1+u1[i])*(X[i,j])+c2*X1_pa[i]+c3*X2_pa[i]+c4*X3_pa[i]+u[i])/(1+exp(c0+(c1+u1[i])*(X[i,j])+c2*X1_pa[i]+c3*X2_pa[i]+c4*X3_pa[i]+u[i]))
      Y[i,j]=rbinom(1, 1, p2[i,j])
    }
  }
  Y_list[[r]] <- Y
} 

print(Y_list[[1]])

###################################################################
dat_pa.cal.sim_list <- list()

for (i in 1:length(Y_list)) {
  Y <- Y_list[[i]]
  X_df <- as.data.frame(X)
  Y_df <- as.data.frame(Y)
  X_df$cffidno <- unique(dat_pa$cffidno)
  Y_df$cffidno <- unique(dat_pa$cffidno)
  
  X_df_long = X_df %>%
    pivot_longer(V1:V126, names_to = "visit", values_to = "VisitAge") 
  
  Y_df_long = Y_df %>%
    pivot_longer(V1:V126, names_to = "visit", values_to = "cltpa.sim") 
  
  dat.sim <- cbind(X_df_long,Y_df_long[,3])
  dat.sim1 <- subset(dat.sim, is.na(VisitAge) != T)
  dat.sim2 <- dat.sim1[with(dat.sim1, order(cffidno, VisitAge)),]
  #dat.sim2$age.bin <- ifelse(dat.sim2$VisitAge!=0, ceiling(dat.sim2$VisitAge),1)
  #dat.sim2$age_yr <- ifelse(dat.sim2$VisitAge!=0, dat.sim2$age.bin-0.5,0)
  dat.sim2$age.bin <- ifelse(dat.sim2$VisitAge!=0, ceiling(dat.sim2$VisitAge*12/2),1)
  dat.sim2$age <- ifelse(dat.sim2$VisitAge!=0, dat.sim2$age.bin*2-1,0)
  dat.sim2$age_yr <- dat.sim2$age/12
  
  ## year
  dat_pa.n.sim <-aggregate(dat.sim2$cltpa.sim, by=list(dat.sim2$age_yr),
                           FUN=sum, na.rm=TRUE)
  names(dat_pa.n.sim) <- c('age_yr','n_pa')
  
  dat_pa.tot.sim <- dat.sim2 %>%
    group_by(age_yr) %>%
    summarise(count = n())
  names(dat_pa.tot.sim) <- c('age_yr','n_visit')
  
  dat_pa.cal.sim0 <- merge(dat_pa.tot.sim, dat_pa.n.sim, by=c('age_yr'))
  dat_pa.cal.sim <- subset(dat_pa.cal.sim0, age_yr<21)
  dat_pa.cal.sim$perc <- dat_pa.cal.sim$n_pa/dat_pa.cal.sim$n_visit
  
  dat_pa.cal.sim_list[[i]] <- dat_pa.cal.sim
}



#################################################################
## --------------------------------------------------------------
####### PA plots - observed + simulated
## --------------------------------------------------------------
plot(dat_pa.cal$age_yr,dat_pa.cal$perc, type = "l", xlab = "Age (Year)", ylab = "Proportion of PA", main="JM3", ylim = c(0,1)) 
for (i in 1:length(Y_list)) {
  lines(dat_pa.cal.sim_list[[i]]$age_yr,dat_pa.cal.sim_list[[i]]$perc, type = "l",col="pink",lwd = 3)
}
lines(dat_pa.cal$age_yr,dat_pa.cal$perc, type = "l")

legend(0,1, c("observed","simulated"), lty=c(1,1),col=c("black","pink"))




#################################################################
## --------------------------------------------------------------
####### PA plots - simulated - JM4
## --------------------------------------------------------------

#################################################################
## --------------------------------------------------------------
####### Posterior results
## --------------------------------------------------------------
result <- read.csv(file="C:/UCHealth/RA/Project/EPIC-CF/Analysis_Jiayuan/Result/Joint_PA_PEX_model_addrandom_rec_model0a_covariate_lt.csv")
c0 <- round(result[1,5],2)
c1 <- round(result[2,5],2)
c2 <- round(result[3,5],2)
c3 <- round(result[4,5],2)
c4 <- round(result[5,5],3) 

u <- round(result[6:1739,5],2)
u.tau.inv <- round(result[1740,5],2)

c=c(c1,c2,c3,c4)

dat_pa2 <- dat_pa[,-1]
head(dat_pa2, n=10)

M <- length(unique(dat_pa2$cffidno))
count <- dat_pa2 %>% count(cffidno)
max.count <- max(count$n)

##########################Assigning unique number to each subject##########################
dat_pa3 <- dat_pa2 %>% group_by(cffidno) %>% mutate(time = c(1:length(cffidno)))

Yd.temp <- data.frame(cffidno = rep(unique(dat_pa2$cffidno),each=max.count), time = 1:max.count)
head(Yd.temp,n=12)

Y.epic <- merge(dat_pa3,Yd.temp,by=c('cffidno','time'),all.y=TRUE)

#################Readingin data for X vectors#############################
X.dat.pa <- dat_pa2[!duplicated(dat_pa2$cffidno,dat_pa2$sexf,dat_pa2$mut), ]
X.dat.pa1 <- model.matrix(cltpa ~ sexf + as.factor(mut), data = X.dat.pa)
X1_pa <- as.numeric(X.dat.pa1[,2]) ## sexf: female
X2_pa <- as.numeric(X.dat.pa1[,3]) ## mut2: 1 alleles D F508
X3_pa <- as.numeric(X.dat.pa1[,4]) ## mut3: 0 alleles D F508

#################Readingin data for VisitAge matrix#############################
X <- matrix(Y.epic$VisitAge, M, max.count, byrow=TRUE)
X[1,]##need to remove first row of column headers
dim(X)
X[1,]

#################input variables for simulation#####################
#### checking for how many individuals we have NAs in the middle of followup
sum.na <- rep(NA,M)
k=rep(NA,M)

ids <- unique(Y.epic$cffidno) ## 103104 103125 103129 103145 103147
for (i in 1:M){
  na.indices <- which(Y.epic$VisitAge[Y.epic$cffidno==ids[i]] %in% NA)
  if (length(na.indices)==0){
    k[i] <- max.count} else{
      k[i] <- min(na.indices)-1}
}
kk=max(k)

###set number of iterations#################################
I=1000

Y_list <- list()
set.seed(123)
#######################################################
for (r in 1:I){
  u<-rnorm(M,0,sqrt(u.tau.inv))
  
  I1<-matrix(NA, nrow=M, ncol=kk, byrow=TRUE)
  I2<-matrix(NA, nrow=M, ncol=kk, byrow=TRUE)
  p2<-matrix(NA, nrow=M, ncol=kk, byrow=TRUE)
  Y<-matrix(NA, nrow=M, ncol=kk, byrow=TRUE)
  
  for (i in 1:M){
    for (j in 1:k[i]){
      I1[i,j]<-ifelse(X[i,j]< cp[1],-1,1)
      I2[i,j]<-ifelse(X[i,j]< cp[2],-1,1)
      p2[i,j]=exp(c0+c1*(X[i,j])+c2*X1_pa[i]+c3*X2_pa[i]+c4*X3_pa[i]+u[i])/(1+exp(c0+c1*(X[i,j])+c2*X1_pa[i]+c3*X2_pa[i]+c4*X3_pa[i]+u[i]))
      Y[i,j]=rbinom(1, 1, p2[i,j])
    }
  }
  Y_list[[r]] <- Y
} 

print(Y_list[[1]])

###################################################################
dat_pa.cal.sim_list <- list()

for (i in 1:length(Y_list)) {
  Y <- Y_list[[i]]
  X_df <- as.data.frame(X)
  Y_df <- as.data.frame(Y)
  X_df$cffidno <- unique(dat_pa$cffidno)
  Y_df$cffidno <- unique(dat_pa$cffidno)
  
  X_df_long = X_df %>%
    pivot_longer(V1:V126, names_to = "visit", values_to = "VisitAge") 
  
  Y_df_long = Y_df %>%
    pivot_longer(V1:V126, names_to = "visit", values_to = "cltpa.sim") 
  
  dat.sim <- cbind(X_df_long,Y_df_long[,3])
  dat.sim1 <- subset(dat.sim, is.na(VisitAge) != T)
  dat.sim2 <- dat.sim1[with(dat.sim1, order(cffidno, VisitAge)),]
  #dat.sim2$age.bin <- ifelse(dat.sim2$VisitAge!=0, ceiling(dat.sim2$VisitAge),1)
  #dat.sim2$age_yr <- ifelse(dat.sim2$VisitAge!=0, dat.sim2$age.bin-0.5,0)
  dat.sim2$age.bin <- ifelse(dat.sim2$VisitAge!=0, ceiling(dat.sim2$VisitAge*12/2),1)
  dat.sim2$age <- ifelse(dat.sim2$VisitAge!=0, dat.sim2$age.bin*2-1,0)
  dat.sim2$age_yr <- dat.sim2$age/12
  
  ## year
  dat_pa.n.sim <-aggregate(dat.sim2$cltpa.sim, by=list(dat.sim2$age_yr),
                           FUN=sum, na.rm=TRUE)
  names(dat_pa.n.sim) <- c('age_yr','n_pa')
  
  dat_pa.tot.sim <- dat.sim2 %>%
    group_by(age_yr) %>%
    summarise(count = n())
  names(dat_pa.tot.sim) <- c('age_yr','n_visit')
  
  dat_pa.cal.sim0 <- merge(dat_pa.tot.sim, dat_pa.n.sim, by=c('age_yr'))
  dat_pa.cal.sim <- subset(dat_pa.cal.sim0, age_yr<21)
  dat_pa.cal.sim$perc <- dat_pa.cal.sim$n_pa/dat_pa.cal.sim$n_visit
  
  dat_pa.cal.sim_list[[i]] <- dat_pa.cal.sim
}



#################################################################
## --------------------------------------------------------------
####### PA plots - observed + simulated
## --------------------------------------------------------------
plot(dat_pa.cal$age_yr,dat_pa.cal$perc, type = "l", xlab = "Age (Year)", ylab = "Proportion of PA", main="JM4", ylim = c(0,1)) 
for (i in 1:length(Y_list)) {
  lines(dat_pa.cal.sim_list[[i]]$age_yr,dat_pa.cal.sim_list[[i]]$perc, type = "l",col="pink",lwd = 3)
}
lines(dat_pa.cal$age_yr,dat_pa.cal$perc, type = "l")

legend(0,1, c("observed","simulated"), lty=c(1,1),col=c("black","pink"))

dev.off()
