# ICPSR 2017 - "Advanced MLE"
#
# A quick primer on how one can transform "non-
# time varying" data into "time varying" data (e.g.,
# in order to use discrete-time approaches).

require(survival)

# Generate some simulated non-time-varying data
# (exponential draws, N=400, 30% censored 
# observations, one binary covariate with B=1)...

set.seed(1337)
N<-400
C<-rbinom(N,1,0.3)
X<-rep(c(0,1),times=N/2)
T<-rexp(N,exp(-1+1*X))
T<-round(T,0)+1
ID<-seq(1:N)

# Data:

NTV<-data.frame(ID=ID,T=T,C=C,X=X)

# Survival object:

Sntv<-Surv(NTV$T,NTV$C)

# Now fit a Cox model...

Cox.ntv<-coxph(Sntv~X,NTV,method="efron")
summary(Cox.ntv)

# Now, expand this into time-varying data:

TV<-data.frame(ID=rep(NTV$ID,NTV$T),
               Tmax=rep(NTV$T,NTV$T),
               X=rep(NTV$X,NTV$T))
TV$T<-with(TV,ave(Tmax,ID,FUN=seq_along)) # time-varying "T"
TV$C<-ifelse(TV$T==TV$Tmax,1,0) # time-varying "C"

# Time-varying survival object:
Stv<-Surv(TV$T-1,TV$T,TV$C)

# Cox model on time-varying data:

Cox.tv<-coxph(Stv~X,TV,method="efron")
summary(Cox.tv)

# Compare (e.g.) the estimated betas and confidence intervals
# from the two models:

pd<-data.frame(B = c(Cox.ntv$coefficients,Cox.tv$coefficients),
               LB = c(confint(Cox.ntv)[1],confint(Cox.tv)[1]),
               UB = c(confint(Cox.ntv)[2],confint(Cox.tv)[2]),
               label = as.factor(c("Non-Time Varying","Time-Varying")))

with(pd, plot(0,0, xlim=c(0.5,2.5),ylim=c(0.3,1.4),type="n",
              xaxt="n",xlab="",ylab="Estimate"))
with(pd, points(label,B,pch=19))
with(pd, segments(1,LB[1],1,UB[1]))
with(pd, segments(2,LB[2],2,UB[2]))
axis(1,at=c(1,2),labels=pd$label)


# IOW, the structure of the data has no effect on the 
# estimates, as long as the underlying information in
# the data is the same.
