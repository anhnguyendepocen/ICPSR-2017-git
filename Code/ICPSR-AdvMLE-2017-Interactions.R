############################################
# ICPSR Advanced MLE -- 2016
#
# Interaction terms in survival models...
############################################

library(gtools)
library(plyr)
library(reshape2)
library(RColorBrewer)
library(lattice)
library(gplots)
library(ggplot2)
library(texreg)
library(survival)

options(scipen = 6) # bias against scientific notation
options(digits = 2) # show fewer decimal places

# Interactions: Simulation

set.seed(1122334)
N<-400
X<-rnorm(N)
W<- -X + rchisq(N,5) # X and W are correlated
Z<-rbinom(N,1,0.5)
XB<-(-0.2 + 0.3*W + (0.5*X) - (0.5*Z) - (1*X*Z))
T<-abs(rexp(N,rate=exp(XB)))
C<-rbinom(N,1,0.8)
df<-data.frame(X=X,W=W,Z=Z,T=T,C=C)
S<-Surv(df$T,df$C)

p <- ggplot(df, aes(x=X,y=log(T),color=as.factor(Z),
                    (shape=20)))
p <- p + theme_classic() + geom_point()+
         geom_smooth(method=loess,se=FALSE,size=1) +
         scale_colour_hue(l=50) +
         ggsave("IntxScatter.pdf",width=5,height=4)
print(p)

# Cox:

coxFit<-coxph(S~W+X+Z+X*Z)
summary(coxFit)

# Weibull:
# WFit<-survreg(S~W+X+Z+X*Z,dist="weibull")
# summary(WFit)

library(car)
linearHypothesis(coxFit,"Z+X:Z")

# Hazard ratios:

VCV<-vcov(coxFit)
coxOR<-data.frame(Z=c(0,1),
        B=c(coxFit$coefficients[2],coxFit$coefficients[2]+coxFit$coefficients[4]),
        se=c(sqrt(VCV[2,2]),sqrt(VCV[2,2]+VCV[4,4]+2*(VCV[2,4]))))
z<-qnorm(0.975)
coxOR$OR<-exp(coxOR$B)
coxOR$UB<-exp(coxOR$B + z*coxOR$se)
coxOR$LB<-exp(coxOR$B - z*coxOR$se)

pdf("InterORs.pdf",5,5)
with(coxOR, plot(Z,OR,ylim=c(0,2.5),xaxt="n",pch=19,
     main="Estimated Odds Ratios for a One-Unit\nChange in X at Z=0 and Z=1"))
with(coxOR, lines(c(0,0),c(UB[1],LB[1]),lwd=2))
with(coxOR, lines(c(1,1),c(UB[2],LB[2]),lwd=2))
abline(h=1,lwd=1,lty=2)
axis(1,at=c(0,1))
dev.off()

# Predicted survival curves:

curves<-5
dfZ0<-data.frame(W=rep(mean(df$W),times=curves),
                 X=quantile(X,probs=seq(0,1,length=curves)),
                 Z=rep(0,times=curves))
SZ0<-survfit(coxFit,dfZ0)

pdf("SZ0.pdf",4,5)
plot(SZ0,mark.time=FALSE,lty=seq(1:5),lwd=2,
     main="Predicted Survival for\nQuintiles of X | Z=0",
     xlab="Time",ylab="Predicted Survival")
legend("topright",lty=seq(1:5),bty="n",lwd=2,
       legend=c("Minimum","25th","Median","75th","Maximum"))
dev.off()

dfZ1<-data.frame(W=rep(mean(df$W),times=curves),
                 X=quantile(X,probs=seq(0,1,length=curves)),
                 Z=rep(1,times=curves))
SZ1<-survfit(coxFit,dfZ1)

pdf("SZ1.pdf",4,5)
plot(SZ1,mark.time=FALSE,lty=seq(1:5),lwd=2,
     main="Predicted Survival for\nQuintiles of X | Z=1",
     xlab="Time",ylab="Predicted Survival")
legend("topright",lty=seq(1:5),bty="n",lwd=2,
       legend=c("Minimum","25th","Median","75th","Maximum"))
dev.off()

##########
# Two Continuous Covariates:

set.seed(2719)
N<-400
X1<-rnorm(N)
W<- -X1 + rchisq(N,5) # X1 and W are correlated
X2<-rpois(N,20)
XB<-(-0.2 + 0.3*W + (0.5*X1) - (0.05*X2) - (0.05*X1*X2))
T<-abs(rexp(N,rate=exp(XB)))
C<-rbinom(N,1,0.8)
df2<-data.frame(W=W,X1=X1,X2=X2,T=T,C=C)
S2<-Surv(df2$T,df2$C)

# Cox:

cox2<-coxph(S2~W+X1+X2+X1*X2,df2)
summary(cox2)

# Heatmap / contour plots...

# Create colorblind-friendly palette:

RdYlGn<-brewer.pal(11,"RdYlGn")
RdBu<-brewer.pal(11,"RdBu")
YlGn<-brewer.pal(9,"YlGn")
Greens<-brewer.pal(9,"Greens")
RWG<-c("#990000","#FFFFFF","#006633") # Red-White-Green
CBFriend<-c("#ca0020","#f4a582","#92c5de","#0571b0")

map<-11
C2sim<-expand.grid(x=quantile(df2$X1,probs=seq(0,1,length=map)),
                   y=quantile(df2$X2,probs=seq(0,1,length=map)))
C2sim$X1<-C2sim$x
C2sim$X2<-C2sim$y
C2sim$W<-mean(df2$W)
C2hats<-survfit(cox2,C2sim)

# Relative risk heatmap:

C2hats2<-predict(cox2,C2sim,type="risk")
HMData<-data.frame(X1=C2sim$X1,X2=C2sim$X2,RR=C2hats2)
HMData2<-melt(HMData,id=c("X1","X2"))
HMData2<-dcast(HMData2,X1~X2)
row.names(HMData2)<-round(HMData2$X1,2)
HMData2$X1<-NULL
colnames(HMData2)<-c("Min","10th","20th","30th","40th","Median",
                     "60th","70th","80th","90th","Max")


pdf("Cox-Heatmap.pdf",10,9)
heatmap.2(as.matrix(HMData2),
          lwid=c(0.05,0.95),
          lhei=c(0.15,0.85),
          Rowv=FALSE,Colv=FALSE,
          dendrogram="none",
          col=colorRampPalette(CBFriend)(40),
          cexRow=1.0,cexCol=1.0,
          trace=c("none"),
          xlab="X2",ylab="X1",
          scale=c("column"),
          srtCol=45,
          key=TRUE)
dev.off()

FLlevel<-

pdf("CoxLevelPlot.pdf",5,4)
with(HMData,levelplot(log(RR)~X1*X2,
            col.regions=colorRampPalette(CBFriend)(100)))
dev.off()


# # Read some data:
# 
# JudgesURL<-"http://www.fjc.gov/history/export/jb.txt"
# temp<-getURL(JudgesURL)
# Judges<-read.csv(textConnection(temp))
# rm(temp)





