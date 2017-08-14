#########################################################
# ICPSR "Advanced Maximum Likelihood" 2017 - Survival
#
# Day six materials.
#
########################################################

library(RCurl)
library(foreign)
library(gtools)
library(plyr)
library(survival)
library(flexsurv)
library(nnet)
library(mstate)
library(texreg)

options(scipen = 5) # bias against scientific notation
options(digits = 2) # show fewer decimal places

# SCOTUS data:

scotusURL<-"https://raw.githubusercontent.com/PrisonRodeo/ICPSR-2017-git/master/Data/scotus.csv"
temp<-getURL(scotusURL)
scotus<-read.csv(textConnection(temp))
rm(temp)

# Stratification:

set.seed=7222009
Z<-rnorm(200)
X0<-rep(0,times=200)
X1<-rep(1,times=200)
T0<-rweibull(200,shape=1,scale=1/exp(1+0.5*Z))
T1<-rweibull(200,shape=0.5,scale=1/exp(1+0.5*Z))
C<-rep(1,times=400)
X<-append(X0,X1)
T<-append(T0,T1)
data<-as.data.frame(cbind(T,C,X,rep(Z,times=2)))
colnames(data)<-c("T","C","X","Z")

S<-Surv(data$T,data$C)

pdf("StratSs.pdf",6,5)
par(mar=c(4,4,2,2))
plot(survfit(S[data$X==0]~1),conf.int=F,lwd=2,
     xlab="Time",ylab="Survival")
lines(survfit(S[data$X==1]~1),conf.int=F,,lwd=2,col="red")
legend("topright",inset=0.1,bty="n",lwd=c(2,2),col=c("black","red"),
       c("p = 1.0","p = 0.5"),cex=1.2)
dev.off()

cox<-coxph(S~Z+X,data=data)
summary(cox)

cox.strata<-coxph(S~Z+strata(X),data=data)
summary(cox.strata)

pdf("FittedStratSs.pdf",6,5)
par(mar=c(4,4,2,2))
plot(survfit(cox),lwd=3,lty=1,col="black",conf.int=F,
     xlab="Time",ylab="Survival",xlim=c(0,3))
lines(survfit(cox.strata),lwd=2,col=c("black","red"),lty=2)
legend("topright",inset=0.05,bty="n",lwd=c(3,2,2),
       col=c("black","black","red"),lty=c(1,2,2),
       c("S(t), No Strata","S(t) | X=0","S(t) | X=1"),cex=1.2)
dev.off()

summary(survreg(S~Z+strata(X),data=data,dist="weibull"))

##################################
# Duration dependence...

# Unobserved heterogeneity plot:

pdf("UnobsHet.pdf",6,5)
par(mar=c(4,4,2,2))
plot(seq(0:20),rep(0.05,times=21),pch=NA,ylim=c(0.015,0.055),
     xlim=c(0,20),xlab="Time",ylab="Hazard")
abline(h=0.02,lwd=3,col="red")
abline(h=0.05,lwd=3,col="black")
abline(0.045,-0.001,lwd=3,lty=2,col="blue")
text(18,0.017,label=c("h(t) | Z=1"),col="red")
text(18,0.053,label=c("h(t) | Z=0"),col="black")
text(15.5,0.035,label=c("Estimated Hazard"),col="blue")
dev.off()

# Heterogeneity sim:

set.seed(7222009)
W<-rnorm(500)
X<-rnorm(500)
Z<-rnorm(500)
T<-rexp(500,rate=(exp(0+0.5*W+0.5*X-0.6*Z)))
C<-rep(1,times=500)
S<-Surv(T,C)

summary(survreg(S~W,dist="weibull"))
summary(survreg(S~W+X,dist="weibull"))
summary(survreg(S~W+X+Z,dist="weibull"))

M1<-survreg(S~W,dist="weibull")
M2<-survreg(S~W+X,dist="weibull")
M3<-survreg(S~W+X+Z,dist="weibull")

# Plot:

t<-cbind(1:60,1:60,1:60)
P<-c(1/(M1$scale),1/(M2$scale),1/(M3$scale))
DurDepHs<-t(apply(t,1,function(t) 0.02*P*((0.02*t)^(P-1))))

pdf("DurDepHs.pdf",6,5)
par(mar=c(4,4,2,2))
plot(t[,1],DurDepHs[,1],t="l",lwd=3,lty=1,col="green",
     xlab="Time",ylab="Weibull Hazard",ylim=c(0.01,0.04))
lines(t[,2],DurDepHs[,2],t="l",lwd=3,lty=2,col="blue")
lines(t[,3],DurDepHs[,3],t="l",lwd=3,lty=3,col="red")
abline(h=0.02,lty=4,lwd=2)
legend("topright",inset=.02,
       c("One Covariate","Two Covariates","Correct Specification","True Hazard"),
       lty=c(1,2,3,4),lwd=c(3,3,3,2),col=c("green","blue","red","black"),
       cex=1.2,bty="n")
dev.off()

# Parameterized duration dependence:

ct.weib<-flexsurvreg(scotus.S~age+pension+pagree,
                     data=scotus,dist="weibull")
ct.weib.DD<-flexsurvreg(scotus.S~age+pension+pagree+shape(age),
                        data=scotus,dist="weibull")

# Plots (this code is a hotter mess than usual, and could use
# a solid smack with a few Hadleyverse tools):

t<-1:max(scotus$service)

Age<-seq(min(scotus$age),max(scotus$age),by=1)
P.vary<-exp(ct.weib.DD$coefficients[1]+(ct.weib.DD$coefficients[6]*Age))

pdf("PbyAge.pdf",6,5)
par(mar=c(4,4,2,2))
plot(Age,P.vary,t="l",lwd=3,lty=1,ylab="Estimate of p")
abline(h=1,lty=2,lwd=2)
dev.off()            

age55<-55
age75<-75

P<-c(exp(scotus.weib$coefficients[1]),
     exp(ct.weib$coefficients[1]),
     exp(ct.weib.DD$coefficients[1]+(ct.weib.DD$coefficients[6]*age55)),
     exp(ct.weib.DD$coefficients[1]+(ct.weib.DD$coefficients[6]*age75)))
XB55<-ct.weib$coefficients[2]+(ct.weib$coefficients[3]*age55)+
  (ct.weib$coefficients[4]*median(scotus$pension))+
  (ct.weib$coefficients[5]*median(scotus$pagree)) 
XB75<-ct.weib$coefficients[2]+(ct.weib$coefficients[3]*age75)+
  (ct.weib$coefficients[4]*median(scotus$pension))+
  (ct.weib$coefficients[5]*median(scotus$pagree)) 
XB55DD<-ct.weib.DD$coefficients[2]+(ct.weib.DD$coefficients[3]*age55)+
  (ct.weib.DD$coefficients[4]*median(scotus$pension))+
  (ct.weib.DD$coefficients[5]*median(scotus$pagree)) 
XB75DD<-ct.weib.DD$coefficients[2]+(ct.weib.DD$coefficients[3]*age75)+
  (ct.weib.DD$coefficients[4]*median(scotus$pension))+
  (ct.weib.DD$coefficients[5]*median(scotus$pagree)) 
XB<-c(XB55,XB75,XB55DD,XB75DD)

h55<-dweibull(t,P[1],scale=XB[1]) / (1 - pweibull(t,P[1],scale=XB[1]))
h75<-dweibull(t,P[2],scale=XB[2]) / (1 - pweibull(t,P[2],scale=XB[2]))
h55DD<-dweibull(t,P[3],scale=XB[3]) / (1 - pweibull(t,P[3],scale=XB[3]))
h75DD<-dweibull(t,P[4],scale=XB[4]) / (1 - pweibull(t,P[4],scale=XB[4]))

pdf("DDHazards.pdf",6,5)
par(mar=c(4,4,2,2))
plot(t,h75DD,t="l",lwd=3,lty=1,col="red",ylim=c(0.2,0.8),
     xlab="Time (in years)",ylab="Hazard")
lines(t,h55DD,lwd=3,lty=1,col="blue")
lines(t,h75,lwd=3,lty=2,col="red")
lines(t,h55,lwd=3,lty=2,col="blue")
legend("topleft",inset=0.02,
       c("Age 55 (p varying)","Age 75 (p varying)","Age 55 (p fixed)",
         "Age 75 (p fixed)"),lty=c(1,1,2,2),lwd=c(3,3,3,3),
       col=c("blue","red","blue","red"),bty="n")
dev.off()

# SCOTUS data...

scotusURL<-"https://raw.githubusercontent.com/PrisonRodeo/ICPSR-2017-git/master/Data/scotus2.csv"
temp<-getURL(scotusURL)
scotus<-read.csv(textConnection(temp))
rm(temp)

scotus.S<-Surv(scotus$svcstart,scotus$service,scotus$retire)


# Illustrating competing risks: 
  
scotus.SR<-Surv(scotus$svcstart,scotus$service,scotus$retire)
scotus.SD<-Surv(scotus$svcstart,scotus$service,scotus$death)

pdf("CR-KMs.pdf",6,5)
plot(survfit(scotus.SR~1),mark.time=F,lwd=c(3,1,1),
     xlab="Time (in years)",ylab="Survival")
par(new=TRUE)
plot(survfit(scotus.SD~1),mark.time=F,lwd=c(4,2,2),
     lty=c(3,3,3),col=c("red","red","red"),
     xlab="Time (in years)",ylab="Survival")
legend("topright",bty="n",inset=0.02,lty=c(1,2),lwd=c(3,4),
       col=c("black","red"),c("Retirement","Death"))
dev.off()

scotus$C<-scotus$retire+scotus$death
scotus.SC<-Surv(scotus$svcstart,scotus$service,scotus$C)

scotus.C<-coxph(scotus.SC~age+chief+south+pension+pagree,
                data=scotus,method="efron")
scotus.R<-coxph(scotus.SR~age+chief+south+pension+pagree,
                data=scotus,method="efron")
scotus.D<-coxph(scotus.SD~age+chief+south+pension+pagree,
                data=scotus,method="efron")

# # Pretty table:
# scotus.C.out<-extract(scotus.C,include.rsquared=F,include.maxrs=F,
#                       include.zph=F,include.nobs=F,
#                       include.missings=F)
# scotus.R.out<-extract(scotus.R,include.rsquared=F,include.maxrs=F,
#                       include.zph=F,include.nobs=F,
#                       include.missings=F)
# scotus.D.out<-extract(scotus.D,include.rsquared=F,include.maxrs=F,
#                       include.zph=F,include.nobs=F,
#                       include.missings=F)
# texreg(l=list(scotus.C.out,scotus.R.out,scotus.D.out),
#        file="scotusCRTable.tex",
#        custom.model.names=c("Combined","Retirement","Death"),
#        custom.coef.names=c("Age","Chief","South",
#                            "Pension Eligibility","Party Agreement"),
#        stars=numeric(0),caption="",label="",table=F)

# MNL:

scotus$lnT<-log(scotus$service)
scotus.MNL<-multinom(threecat~chief+south+age+pension+pagree+lnT,
                     data=scotus,na.action=na.omit)

# # Pretty table:
# scotus.MNL.out<-extract(scotus.MNL,include.deviance=F,include.maxrs=F,
#                       include.zph=F,include.nobs=F,
#                       include.missings=F)
# texreg(l=scotus.MNL.out,
#        file="scotusMNLTable.tex",       
#        custom.model.names=c("Retirement","Death"),
#        custom.coef.names=c("Intercept","Age","Chief","South",
#                            "Pension Eligibility","Party Agreement",
#                            "log(Time)"),
#        stars=numeric(0),caption="",label="",table=F)


# Multiple Events...

ORURL<-"https://raw.githubusercontent.com/PrisonRodeo/ICPSR-2017-git/master/Data/OR.csv"
temp<-getURL(ORURL)
OR<-read.csv(textConnection(temp))
rm(temp)

OR<-OR[order(OR$dyadid,OR$year),]
OR$one<-rep(1,times=nrow(OR))
OR<-ddply(OR,"dyadid",mutate,eventno=cumsum(dispute)-dispute+1,
          altstart=cumsum(one)-1,altstop=cumsum(one))

listvars<-c("dyadid","year","start","stop",
       "altstart","altstop","dispute","eventno")
ORsm<-OR[listvars]
print(ORsm[which(ORsm$dyadid==2130 & ORsm$year<1966),])

# First events:

OR1st<-OR[OR$eventno==1,]
OR.1st<-Surv(OR1st$altstart,OR1st$altstop,OR1st$dispute)
OR.Cox.1st<-coxph(OR.1st~allies+contig+capratio+growth+democracy+
                   trade+cluster(dyadid),data=OR1st,method="efron")

# AG:

OR.AGS<-Surv(OR$altstart,OR$altstop,OR$dispute)
OR.Cox.AG<-coxph(OR.AGS~allies+contig+capratio+growth+democracy+
                   trade+cluster(dyadid),data=OR,method="efron")

# PWP Elapsed:

OR.PWPES<-Surv(OR$altstart,OR$altstop,OR$dispute)
OR.Cox.PWPE<-coxph(OR.PWPES~allies+contig+capratio+growth+democracy+
                   trade+strata(eventno)+cluster(dyadid),data=OR,
                   method="efron")

# PWP Gap:

OR.PWPGS<-Surv(OR$start,OR$stop,OR$dispute)
OR.Cox.PWPG<-coxph(OR.PWPGS~allies+contig+capratio+growth+democracy+
                     trade+strata(eventno)+cluster(dyadid),data=OR,
                     method="efron")

# WLW:

OR.expand<-OR[rep(1:nrow(OR),each=max(OR$eventno)),]
OR.expand<-ddply(OR.expand,c("dyadid","year"),mutate,
                 eventrisk=cumsum(one))
OR.expand$dispute<-ifelse(OR.expand$eventno==OR.expand$eventrisk 
                          & OR.expand$dispute==1,1,0)

OR.expand.S<-Surv(OR.expand$altstart,OR.expand$altstop,
                  OR.expand$dispute)
OR.Cox.WLW<-coxph(OR.expand.S~allies+contig+capratio+growth+
                    democracy+trade+strata(eventno)+
                    cluster(dyadid),data=OR.expand,
                    method="efron")

# # Pretty table:
# 
# OR.1st.out<-extract(OR.Cox.1st,include.rsquared=F,include.maxrs=F,
#                       include.zph=F,include.nobs=F,
#                       include.missings=F)
# OR.AG.out<-extract(OR.Cox.AG,include.rsquared=F,include.maxrs=F,
#                     include.zph=F,include.nobs=F,
#                     include.missings=F)
# OR.PWPE.out<-extract(OR.Cox.PWPE,include.rsquared=F,include.maxrs=F,
#                     include.zph=F,include.nobs=F,
#                     include.missings=F)
# OR.PWPG.out<-extract(OR.Cox.PWPG,include.rsquared=F,include.maxrs=F,
#                     include.zph=F,include.nobs=F,
#                     include.missings=F)
# OR.WLW.out<-extract(OR.Cox.WLW,include.rsquared=F,include.maxrs=F,
#                     include.zph=F,include.nobs=F,
#                     include.missings=F)
# 
# texreg(l=list(OR.1st.out,OR.AG.out,OR.PWPE.out,OR.PWPG.out,OR.WLW.out),
#        file="RepeatedEventsTable.tex",
#        custom.model.names=c("First","AG","PWP-E","PWP-G","WLW"),
#        custom.coef.names=c("Allies","Contiguity","Capability Ratio",
#                            "Growth","Democracy","Trade"),
#        stars=numeric(0),caption="",label="",table=F)

# Parameter change

OR$capXevent<-OR$capratio*OR$eventno
OR.Cox.BVary<-coxph(OR.PWPGS~allies+contig+growth+democracy+
                     trade+capratio+capXevent+strata(eventno)+
                     cluster(dyadid),data=OR,
                     method="efron")

