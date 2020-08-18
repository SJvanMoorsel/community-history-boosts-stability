####Revision Nov. 2019
#------------------------------------------------------------
#Stability analysis  May 2019
rm(list=ls())
ls()
library(pascal)
library(asreml) #https://www.vsni.co.uk/software/asreml
library(ggplot2)
##### LOADING DATA
# Modify to your directory
d=read.csv("YourDirectory/Data/StabilityMetrics.csv",sep=";")
dim(d)
str(d)
names(d)
head(d)
#--------------------------------
# ----- This plot should be removed!
# Plot!="B1A12"|Soil!="Org"
d=d[- which(d$Plot=="B1A12"&d$SH=="Org"),] #(already removed in dataset available on Pangaea)
d=d[d$Plot!="B2A05",] #excluding fes pra monoculture plots (already removed in dataset available on Pangaea)
# ----- Remove extreme outlier Ono.vic mono:
d=d[d$Plot != "B2A15",] 
dim(d)
# transforming variables, making factors and contrasts:
d$logstab=log(d$Stab);d$logstab2=log(d$noFloodStab)
d$popstab=1/d$popCV;d$popstab2=1/d$noFlood.popCV
d$logpopstab=log(1/d$popCV);d$logpopstab2=log(1/d$noFlood.popCV)
d$async=1-d$Sync;d$async2=1-d$noFlood.Sync
d$logSR=log2(d$SR);d$fSR=factor(d$SR)
d$leg=factor(as.numeric(d$nLegume>0),labels=c("no_leg","legumes"));d$propL=d$nLegume/d$SR
d$gr=factor(as.numeric(d$nGrass>0),labels=c("no_gr","grass"));d$propG=d$nGrass/d$SR
d$th=factor(as.numeric(d$nTherb>0),labels=c("no_th","tall"));d$propT=d$nTherb/d$SR
d$sh=factor(as.numeric(d$nSherb>0),labels=c("no_sh","short"));d$propS=d$nSherb/d$SR
#--------------------------------
# ANOVAs for Table 1  (six harvests, excluding flood harvest)
#--------------------------------
m1=asreml(logstab2~1+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH,data=d)
test.asreml(m1)
m2=asreml(logpopstab2~1+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH,data=d)
test.asreml(m2)
m3=asreml.nvc(async2~1+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH,data=d)
test.asreml(m3)
#--------------------------------
# % sum of squares for Table 1 (six harvests, excluding flood harvest)
#--------------------------------
lm1=lm(terms(logstab2~1+logSR+SH+PH+SH:logSR+PH:logSR+Plot+Plot:SH,keep.order=T),data=d)
anova(lm1)
lm2=lm(terms(logpopstab2~1+logSR+SH+PH+SH:logSR+PH:logSR+Plot+Plot:SH,keep.order=T),data=d)
anova(lm2)
lm3=lm(terms(async2~1+logSR+SH+PH+SH:logSR+PH:logSR+Plot+Plot:SH,keep.order=T),data=d)
anova(lm3)
#--------------------------------
# ANOVAs for Table 2 
#--------------------------------
m1=asreml(resistRaw~1+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH ,data=d)
test.asreml(m1)
m1=asreml.nvc(recovRaw~1+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH ,data=d)
test.asreml(m1)
m1=asreml(resilRaw~1+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH ,data=d)
test.asreml(m1)
#--------------------------------
# % sum of squares for Table 2
#--------------------------------
lm1=lm(terms(resistRaw~1+logSR+SH+PH+SH:logSR+PH:logSR+Plot+Plot:SH,keep.order=T),data=d)
anova(lm1)
lm2=lm(terms(recovRaw~1+logSR+SH+PH+SH:logSR+PH:logSR+Plot+Plot:SH,keep.order=T),data=d)
anova(lm2)
lm3=lm(terms(resilRaw~1+logSR+SH+PH+SH:logSR+PH:logSR+Plot+Plot:SH,keep.order=T),data=d)
anova(lm3)
#--------------------------------
# ANOVAs for Table 3 
#--------------------------------
dnoNA=d[!(is.na(d$stabBefore)),]
dnoNA$logstabBefore=log(dnoNA$stabBefore)
mm=asreml(logstabBefore~1+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH,data=dnoNA)
test.asreml(mm)
dnoNA=d[!(is.na(d$stabAfter)),]
dnoNA$logstabAfter=log(dnoNA$stabAfter)
mm=asreml.nvc(logstabAfter~1+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH,data=dnoNA)
test.asreml(mm)
#--------------------------------
# ANOVAs for Table S1  (all seven harvests)
#--------------------------------
m1=asreml(logstab~1+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH,data=d)
test.asreml(m1)
m2=asreml(logpopstab~1+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH,data=d)
test.asreml(m2)
m3=asreml.nvc(async~1+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH,data=d)
test.asreml(m3)
#--------------------------------
# ANOVAs for Table S2 
#--------------------------------
m1=asreml.nvc(turnBefAft~1+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH,data=d)
test.asreml(m1)
#--------------------------------
# ANOVAs for Table S3 
#--------------------------------
m1=asreml(resistRaw~1+BMbefore+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH ,data=d)
test.asreml(m1)
m1=asreml.nvc(recovRaw~1+BMbefore+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH ,data=d)
test.asreml(m1)
m1=asreml(resilRaw~1+BMbefore+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH ,data=d)
test.asreml(m1)
#--------------------------------
# Figure 2A - Diversity stability by plant history (6 harvests)
#--------------------------------
par(mfrow=c(2,2))
par(mar=c(5,4,4,2)+0.1,oma=c(1,2,1,1)+0.1)
mm=asreml(logstab2~1+logSR+PH+SH+logSR:PH+logSR:SH
          ,random=~Plot/SH,data=d)
test.asreml(mm)
mtemp=lm(terms(logstab2~Plot:SH,keep.order=T),data=d)
d$res_logstab2=residuals(mtemp)
for (i in c(1,2,4,8)) {
  d1=subset(d,d$SR==i)
  d$adj_logstab2[d$SR==i]=d1$res_logstab2+mean(d1$logstab2)
}
x.old=tapply(d[d$PH=="old",]$logSR,d[d$PH=="old",]$Plot,mean,na.rm=T)
y.old=tapply(d[d$PH=="old",]$noFloodStab,d[d$PH=="old",]$Plot,mean,na.rm=T)
x.new=tapply(d[d$PH=="new",]$logSR,d[d$PH=="new",]$Plot,mean,na.rm=T)
y.new=tapply(d[d$PH=="new",]$noFloodStab,d[d$PH=="new",]$Plot,mean,na.rm=T)
xi=aggregate(d$logSR,list(Plot=d$Plot,PH=d$PH),mean)
yi=aggregate(exp(d$adj_logstab2),list(Plot=d$Plot,PH=d$PH),mean)
plot(x=jitter(xi[xi$PH=="old","x"]-0.05,factor=0.2),y=yi[yi$PH=="old","x"]
     ,xaxt="n",bty="n",tck=0.02,mgp=c(1.5,0.1,0)
     ,ylab=bquote(paste("Stability (", italic("CV"["com"])^-1, ")"))
     ,xlab="Planted richness"
     ,pch=21,col=alpha("black",0.6),bg=alpha("skyblue",0.6),lwd=1.5
     ,cex=0.9,tck=0.02,cex.axis=1.2,cex.lab=1.5
     ,ylim=c(0.5,4),xlim=c(-0.1,log2(8)+0.1),log="y")
points(x=jitter(xi[xi$PH=="new","x"]+0.05,factor=0.2),yi[yi$PH=="new","x"]
       ,pch=21,col=alpha("black",0.6),bg=alpha("firebrick",0.6),lwd=1.5,cex=0.9)
axis(side=1,at=log2(c(1,2,4,8)),labels=c(1,2,4,8),tck=0.02,mgp=c(1.5,0.1,0),cex.axis=1.2)
pr.old=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.old,na.rm=T),max(x.old,na.rm=T),l=30)
                                           ,PH="old"))$prediction$pvals
pr.new=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.new,na.rm=T),max(x.new,na.rm=T),l=30)
                                           ,PH="new"))$prediction$pvals
polygon(c(pr.new$logSR,rev(pr.new$logSR))
        ,exp(c(pr.new[,3]-pr.new[,4],rev(pr.new[,3]+pr.new[,4])))
        ,border=NA,col=alpha("firebrick",0.2))
polygon(c(pr.old$logSR,rev(pr.old$logSR))
        ,exp(c(pr.old[,3]-pr.old[,4],rev(pr.old[,3]+pr.old[,4])))
        ,border=NA,col=alpha("skyblue",0.2))
lines(pr.old$logSR,exp(pr.old[,3]),lwd=3,col="skyblue")
lines(pr.new$logSR,exp(pr.new[,3]),lwd=3,col="firebrick")
#--------------------------------
# Figure 2B  - Diversity pop. stability by plant history (6 harvests)
#--------------------------------
mm=asreml(logpopstab2~1+logSR+PH+SH+logSR:PH+logSR:SH
          ,random=~Plot/SH,data=d)
test.asreml(mm)
mtemp=lm(terms(logpopstab2~Plot:SH,keep.order=T),data=d)
d$res_logpopstab2=residuals(mtemp)
for (i in c(1,2,4,8)) {
  d1=subset(d,d$SR==i)
  d$adj_logpopstab2[d$SR==i]=d1$res_logpopstab2+mean(d1$logpopstab2)
}
x.old=tapply(d[d$PH=="old",]$logSR,d[d$PH=="old",]$Plot,mean,na.rm=T)
y.old=tapply(d[d$PH=="old",]$popstab2,d[d$PH=="old",]$Plot,mean,na.rm=T)
x.new=tapply(d[d$PH=="new",]$logSR,d[d$PH=="new",]$Plot,mean,na.rm=T)
y.new=tapply(d[d$PH=="new",]$popstab2,d[d$PH=="new",]$Plot,mean,na.rm=T)
xi=aggregate(d$logSR,list(Plot=d$Plot,PH=d$PH),mean)
yi=aggregate(exp(d$adj_logpopstab2),list(Plot=d$Plot,PH=d$PH),mean)
#par(mfrow=c(1,1))
plot(x=jitter(xi[xi$PH=="old","x"]-0.05,factor=0.2),y=yi[yi$PH=="old","x"]
     ,xaxt="n",bty="n",tck=0.02,mgp=c(1.5,0.1,0)
     ,ylab=bquote(paste("Population stability (", italic("CV"["pop"])^-1, ")"))
     ,xlab="Planted richness"
     ,pch=21,col=alpha("black",0.6),bg=alpha("skyblue",0.6),lwd=1.5
     ,cex=0.9,tck=0.02,cex.axis=1.2,cex.lab=1.5
     ,ylim=c(0.6,2),xlim=c(-0.1,log2(8)+0.1),log="y")
points(x=jitter(xi[xi$PH=="new","x"]+0.05,factor=0.2),yi[yi$PH=="new","x"]
       ,pch=21,col=alpha("black",0.6),bg=alpha("firebrick",0.6),lwd=1.5,cex=0.9)
axis(side=1,at=log2(c(1,2,4,8)),labels=c(1,2,4,8),tck=0.02,mgp=c(1.5,0.1,0),cex.axis=1.2)
pr.old=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.old,na.rm=T),max(x.old,na.rm=T),l=30)
                                           ,PH="old"))$prediction$pvals
pr.new=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.new,na.rm=T),max(x.new,na.rm=T),l=30)
                                           ,PH="new"))$prediction$pvals
polygon(c(pr.new$logSR,rev(pr.new$logSR))
        ,exp(c(pr.new[,3]-pr.new[,4],rev(pr.new[,3]+pr.new[,4])))
        ,border=NA,col=alpha("firebrick",0.2))
polygon(c(pr.old$logSR,rev(pr.old$logSR))
        ,exp(c(pr.old[,3]-pr.old[,4],rev(pr.old[,3]+pr.old[,4])))
        ,border=NA,col=alpha("skyblue",0.2))
lines(pr.old$logSR,exp(pr.old[,3]),lwd=3,col="skyblue")
lines(pr.new$logSR,exp(pr.new[,3]),lwd=3,col="firebrick")
#--------------------------------
# Figure 2C  - Diversity asynchrony by plant history (6 harvests)
#--------------------------------
mm=asreml(async2~1+logSR+PH+SH+logSR:PH+logSR:SH
          ,random=~Plot/SH,data=d)
test.asreml(mm)
mtemp=lm(terms(async2~Plot:SH,keep.order=T),data=d)
d$res_async2=residuals(mtemp)
for (i in c(1,2,4,8)) {
  d1=subset(d,d$SR==i)
  d$adj_async2[d$SR==i]=d1$res_async2+mean(d1$async2)
}
x.old=tapply(d[d$PH=="old",]$logSR,d[d$PH=="old",]$Plot,mean,na.rm=T)
y.old=tapply(d[d$PH=="old",]$async2,d[d$PH=="old",]$Plot,mean,na.rm=T)
x.new=tapply(d[d$PH=="new",]$logSR,d[d$PH=="new",]$Plot,mean,na.rm=T)
y.new=tapply(d[d$PH=="new",]$async2,d[d$PH=="new",]$Plot,mean,na.rm=T)
xi=aggregate(d$logSR,list(Plot=d$Plot,PH=d$PH),mean)
yi=aggregate(d$adj_async2,list(Plot=d$Plot,PH=d$PH),mean)
#par(mfrow=c(1,1))
plot(x=jitter(xi[xi$PH=="old","x"]-0.05,factor=0.2),y=yi[yi$PH=="old","x"]
     ,xaxt="n",bty="n",tck=0.02,mgp=c(1.5,0.1,0)
     ,ylab=bquote(paste("Asynchrony (", 1-theta, ")"))
     ,xlab="Planted richness"
     ,pch=21,col=alpha("black",0.6),bg=alpha("skyblue",0.6),lwd=1.5
     ,cex=0.9,tck=0.02,cex.axis=1.2,cex.lab=1.5
     ,ylim=c(0,1),xlim=c(-0.1,log2(8)+0.1))
points(x=jitter(xi[xi$PH=="new","x"]+0.05,factor=0.2),yi[yi$PH=="new","x"]
       ,pch=21,col=alpha("black",0.6),bg=alpha("firebrick",0.6),lwd=1.5,cex=0.9)
axis(side=1,at=log2(c(1,2,4,8)),labels=c(1,2,4,8),tck=0.02,mgp=c(1.5,0.1,0),cex.axis=1.2)
pr.old=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.old,na.rm=T),max(x.old,na.rm=T),l=30)
                                           ,PH="old"))$prediction$pvals
pr.new=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.new,na.rm=T),max(x.new,na.rm=T),l=30)
                                           ,PH="new"))$prediction$pvals
polygon(c(pr.new$logSR,rev(pr.new$logSR))
        ,c(pr.new[,3]-pr.new[,4],rev(pr.new[,3]+pr.new[,4]))
        ,border=NA,col=alpha("firebrick",0.2))
polygon(c(pr.old$logSR,rev(pr.old$logSR))
        ,c(pr.old[,3]-pr.old[,4],rev(pr.old[,3]+pr.old[,4]))
        ,border=NA,col=alpha("skyblue",0.2))
lines(pr.old$logSR,pr.old[,3],lwd=3,col="skyblue")
lines(pr.new$logSR,pr.new[,3],lwd=3,col="firebrick")
#--------------------------------
# Figure 2D  - Residuals Stability vs. residuals asynchrony
#--------------------------------
d$noFloodAsync=1-d$noFlood.Sync
m2=asreml.nvc(logstab2~1+noFloodAsync+logSR+PH+noFloodAsync:logSR+noFloodAsync:PH
              ,random=~Plot/SH+noFloodAsync:Plot,data=d)
test.asreml(m2)
# correlation with residuals:
mtemp=lm(terms(logstab2~Plot:SH,keep.order=T),data=d)
d$res_logstab2=residuals(mtemp)
mtemp=lm(terms(noFloodAsync~Plot:SH,keep.order=T),data=d)
d$res_async2=residuals(mtemp)
#par(mfrow=c(1,1))
#plot(x=d$res_async2,y=d$res_logstab2)
lm=lm(res_logstab2~res_async2+PH,data=d) # res_async2:PH is zero because of adjustment
anova(lm) #both significant (P<0.001 and P<0.01)
#but residual df in above should be corrected from 271 to 135:
lm=lm(res_logstab2~res_async2+Plot/SH+PH,data=d) # res_async2:PH is zero because of adjustment
anova(lm)
pf(8.1272,1,135,lower.tail=F)
mm=asreml(res_logstab2~res_async2+PH,data=d) # res_async2:PH is zero because of adjustment
test.asreml(mm)
x.old=tapply(d[d$PH=="old",]$res_async2,list(d[d$PH=="old",]$Plot,d[d$PH=="old",]$SH),mean,na.rm=T)
y.old=tapply(d[d$PH=="old",]$res_logstab2,list(d[d$PH=="old",]$Plot,d[d$PH=="old",]$SH),mean,na.rm=T)
x.new=tapply(d[d$PH=="new",]$res_async2,list(d[d$PH=="new",]$Plot,d[d$PH=="new",]$SH),mean,na.rm=T)
y.new=tapply(d[d$PH=="new",]$res_logstab2,list(d[d$PH=="new",]$Plot,d[d$PH=="new",]$SH),mean,na.rm=T)
xi=aggregate(d$res_async2,list(Plot=d$Plot,SH=d$SH,PH=d$PH),mean)
yi=aggregate(d$res_logstab2,list(Plot=d$Plot,SH=d$SH,PH=d$PH),mean)
#par(mfrow=c(1,1))
plot(x=xi[xi$PH=="old","x"],y=yi[yi$PH=="old","x"]
     ,xaxt="n",bty="n",tck=0.02,mgp=c(1.5,0.1,0)
     ,ylab="Stability residuals"
     ,xlab="Asynchrony residuals"
     ,pch=21,col=alpha("black",0.6),bg=alpha("skyblue",0.6),lwd=1.5
     ,cex=0.9,tck=0.02,cex.axis=1.2,cex.lab=1.5
     ,ylim=c(-0.7,0.7),xlim=c(-0.3,0.3))
points(x=xi[xi$PH=="new","x"],yi[yi$PH=="new","x"]
       ,pch=21,col=alpha("black",0.6),bg=alpha("firebrick",0.6),lwd=1.5,cex=0.9)
axis(side=1,at=c(-0.3,-0.15,0,0.15,0.3),labels=c(-0.3,-0.15,0,0.15,0.3),tck=0.02,mgp=c(1.5,0.1,0),cex.axis=1.2)
pr.old=predict(mm,classify="PH:res_async2",list(res_async2=seq(min(x.old,na.rm=T),max(x.old,na.rm=T),l=30)
                                           ,PH="old"))$prediction$pvals
pr.new=predict(mm,classify="PH:res_async2",list(res_async2=seq(min(x.new,na.rm=T),max(x.new,na.rm=T),l=30)
                                           ,PH="new"))$prediction$pvals
polygon(c(pr.new$res_async2,rev(pr.new$res_async2))
        ,c(pr.new[,3]-pr.new[,4],rev(pr.new[,3]+pr.new[,4]))
        ,border=NA,col=alpha("firebrick",0.2))
polygon(c(pr.old$res_async2,rev(pr.old$res_async2))
        ,c(pr.old[,3]-pr.old[,4],rev(pr.old[,3]+pr.old[,4]))
        ,border=NA,col=alpha("skyblue",0.2))
lines(pr.old$res_async2,pr.old[,3],lwd=3,col="skyblue")
lines(pr.new$res_async2,pr.new[,3],lwd=3,col="firebrick")
#--------------------------------
# Figure 4A  - Diversity resistance by plant history
#--------------------------------
par(mfrow=c(1,3))
par(mar=c(5,6,4,2)+0.1,oma=c(1,1,1,1)+0.1)
dnoNA=d[!(is.na(d$resistRaw)),]
mm=asreml(resistRaw~1+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH ,data=dnoNA)
test.asreml(mm)
x.old=tapply(dnoNA[dnoNA$PH=="old",]$logSR,dnoNA[dnoNA$PH=="old",]$Plot,mean,na.rm=T)
y.old=tapply(dnoNA[dnoNA$PH=="old",]$resistRaw,dnoNA[dnoNA$PH=="old",]$Plot,mean,na.rm=T)
x.new=tapply(dnoNA[dnoNA$PH=="new",]$logSR,dnoNA[dnoNA$PH=="new",]$Plot,mean,na.rm=T)
y.new=tapply(dnoNA[dnoNA$PH=="new",]$resistRaw,dnoNA[dnoNA$PH=="new",]$Plot,mean,na.rm=T)
mtemp=lm(terms(resistRaw~Plot:SH,keep.order=T),data=dnoNA)
dnoNA$res_resistRaw=residuals(mtemp)
dim(dnoNA)
for (i in c(1,2,4,8)) {
  d1=subset(dnoNA,dnoNA$SR==i)
  dnoNA$adj_resistRaw[dnoNA$SR==i]=d1$res_resistRaw+mean(d1$resistRaw)
}
xi=aggregate(dnoNA$logSR,list(Plot=dnoNA$Plot,PH=dnoNA$PH),mean)
yi=aggregate(dnoNA$adj_resistRaw,list(Plot=dnoNA$Plot,PH=dnoNA$PH),mean)
plot(jitter(x=xi[xi$PH=="old","x"]-0.05,factor=0.2),y=yi[yi$PH=="old","x"]
     ,xaxt="n",bty="n",tck=0.03,mgp=c(3,0.5,0) #default 3,1,0 (1.5,0.1,0)
     ,ylab=bquote(paste("Resistance (g/m"^2, ")"))
     ,xlab="Planted richness",xlim=c(-0.1,log2(8)+0.1)
     ,pch=21,col=alpha("black",0.6),bg=alpha("skyblue",0.6),lwd=1.5
     ,cex=1.5,tck=0.03,cex.axis=1.7,cex.lab=2
     ,ylim=c(-240,60))
points(x=jitter(xi[xi$PH=="new","x"]+0.05,factor=0.2),y=yi[yi$PH=="new","x"]
       ,pch=21,col=alpha("black",0.6),bg=alpha("firebrick",0.6),lwd=1.5,cex=1.5)
axis(side=1,at=log2(c(1,2,4,8)),labels=c(1,2,4,8),tck=0.03,mgp=c(3,0.5,0),cex.axis=1.7)
pr.old=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.old,na.rm=T),max(x.old,na.rm=T),l=30)
                                           ,PH="old"))$prediction$pvals
pr.new=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.new,na.rm=T),max(x.new,na.rm=T),l=30)
                                           ,PH="new"))$prediction$pvals
polygon(c(pr.new$logSR,rev(pr.new$logSR))
        ,c(pr.new[,3]-pr.new[,4],rev(pr.new[,3]+pr.new[,4]))
        ,border=NA,col=alpha("firebrick",0.2))
polygon(c(pr.old$logSR,rev(pr.old$logSR))
        ,c(pr.old[,3]-pr.old[,4],rev(pr.old[,3]+pr.old[,4]))
        ,border=NA,col=alpha("skyblue",0.2))
lines(pr.old$logSR,pr.old[,3],lwd=3,col="skyblue")
lines(pr.new$logSR,pr.new[,3],lwd=3,col="firebrick")
lines(c(-1,10), c(0,0), lty=3, col="black")
#--------------------------------
# Figure 4B  - Diversity recovery by plant history
#--------------------------------
dnoNA=d[!(is.na(d$recovRaw)),]
mm=asreml(recovRaw~1+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH ,data=dnoNA)
test.asreml(mm)
x.old=tapply(dnoNA[dnoNA$PH=="old",]$logSR,dnoNA[dnoNA$PH=="old",]$Plot,mean,na.rm=T)
y.old=tapply(dnoNA[dnoNA$PH=="old",]$recovRaw,dnoNA[dnoNA$PH=="old",]$Plot,mean,na.rm=T)
x.new=tapply(dnoNA[dnoNA$PH=="new",]$logSR,dnoNA[dnoNA$PH=="new",]$Plot,mean,na.rm=T)
y.new=tapply(dnoNA[dnoNA$PH=="new",]$recovRaw,dnoNA[dnoNA$PH=="new",]$Plot,mean,na.rm=T)
mtemp=lm(terms(recovRaw~Plot:SH,keep.order=T),data=dnoNA)
dnoNA$res_recovRaw=residuals(mtemp)
dim(dnoNA)
for (i in c(1,2,4,8)) {
  d1=subset(dnoNA,dnoNA$SR==i)
  dnoNA$adj_recovRaw[dnoNA$SR==i]=d1$res_recovRaw+mean(d1$recovRaw)
}
xi=aggregate(dnoNA$logSR,list(Plot=dnoNA$Plot,PH=dnoNA$PH),mean)
yi=aggregate(dnoNA$adj_recovRaw,list(Plot=dnoNA$Plot,PH=dnoNA$PH),mean)
plot(jitter(x=xi[xi$PH=="old","x"]-0.05,factor=0.2),y=yi[yi$PH=="old","x"]
     ,xaxt="n",bty="n",tck=0.03,mgp=c(3,0.5,0)
     ,ylab=bquote(paste("Recovery (g/m"^2, ")"))
     ,xlab="Planted richness",xlim=c(-0.1,log2(8)+0.1)
     ,pch=21,col=alpha("black",0.6),bg=alpha("skyblue",0.6),lwd=1.5
     ,cex=1.5,tck=0.03,cex.axis=1.7,cex.lab=2
     ,ylim=c(-50,250))
points(x=jitter(xi[xi$PH=="new","x"]+0.05,factor=0.2),y=yi[yi$PH=="new","x"]
       ,pch=21,col=alpha("black",0.6),bg=alpha("firebrick",0.6),lwd=1.5,cex=1.5)
axis(side=1,at=log2(c(1,2,4,8)),labels=c(1,2,4,8),tck=0.03,mgp=c(3,0.5,0),cex.axis=1.7)
pr.old=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.old,na.rm=T),max(x.old,na.rm=T),l=30)
                                           ,PH="old"))$prediction$pvals
pr.new=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.new,na.rm=T),max(x.new,na.rm=T),l=30)
                                           ,PH="new"))$prediction$pvals
polygon(c(pr.new$logSR,rev(pr.new$logSR))
        ,c(pr.new[,3]-pr.new[,4],rev(pr.new[,3]+pr.new[,4]))
        ,border=NA,col=alpha("firebrick",0.2))
polygon(c(pr.old$logSR,rev(pr.old$logSR))
        ,c(pr.old[,3]-pr.old[,4],rev(pr.old[,3]+pr.old[,4]))
        ,border=NA,col=alpha("skyblue",0.2))
lines(pr.old$logSR,pr.old[,3],lwd=3,col="skyblue")
lines(pr.new$logSR,pr.new[,3],lwd=3,col="firebrick")
lines(c(-1,10), c(0,0), lty=3, col="black")
#--------------------------------
# Figure 4C  - Diversity resilience by plant history
#--------------------------------
dnoNA=d[!(is.na(d$resilRaw)),]
mm=asreml(resilRaw~1+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH ,data=dnoNA)
test.asreml(mm)
x.old=tapply(dnoNA[dnoNA$PH=="old",]$logSR,dnoNA[dnoNA$PH=="old",]$Plot,mean,na.rm=T)
y.old=tapply(dnoNA[dnoNA$PH=="old",]$resilRaw,dnoNA[dnoNA$PH=="old",]$Plot,mean,na.rm=T)
x.new=tapply(dnoNA[dnoNA$PH=="new",]$logSR,dnoNA[dnoNA$PH=="new",]$Plot,mean,na.rm=T)
y.new=tapply(dnoNA[dnoNA$PH=="new",]$resilRaw,dnoNA[dnoNA$PH=="new",]$Plot,mean,na.rm=T)
mtemp=lm(terms(resilRaw~Plot:SH,keep.order=T),data=dnoNA)
dnoNA$res_resilRaw=residuals(mtemp)
dim(dnoNA)
for (i in c(1,2,4,8)) {
  d1=subset(dnoNA,dnoNA$SR==i)
  dnoNA$adj_resilRaw[dnoNA$SR==i]=d1$res_resilRaw+mean(d1$resilRaw)
}
xi=aggregate(dnoNA$logSR,list(Plot=dnoNA$Plot,PH=dnoNA$PH),mean)
yi=aggregate(dnoNA$adj_resilRaw,list(Plot=dnoNA$Plot,PH=dnoNA$PH),mean)
plot(jitter(x=xi[xi$PH=="old","x"]-0.05,factor=0.2),y=yi[yi$PH=="old","x"]
     ,xaxt="n",bty="n",tck=0.03,mgp=c(3,0.5,0)
     ,ylab=bquote(paste("Resilience (g/m"^2, ")"))
     ,xlab="Planted richness",xlim=c(-0.1,log2(8)+0.1)
     ,pch=21,col=alpha("black",0.6),bg=alpha("skyblue",0.6),lwd=1.5
     ,cex=1.5,tck=0.03,cex.axis=1.7,cex.lab=2
     ,ylim=c(-220,220))
points(x=jitter(xi[xi$PH=="new","x"]+0.05,factor=0.2),y=yi[yi$PH=="new","x"]
       ,pch=21,col=alpha("black",0.6),bg=alpha("firebrick",0.6),lwd=1.5,cex=1.5)
axis(side=1,at=log2(c(1,2,4,8)),labels=c(1,2,4,8),tck=0.02,mgp=c(3,0.5,0),cex.axis=1.7)
pr.old=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.old,na.rm=T),max(x.old,na.rm=T),l=30)
                                           ,PH="old"))$prediction$pvals
pr.new=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.new,na.rm=T),max(x.new,na.rm=T),l=30)
                                           ,PH="new"))$prediction$pvals
polygon(c(pr.new$logSR,rev(pr.new$logSR))
        ,c(pr.new[,3]-pr.new[,4],rev(pr.new[,3]+pr.new[,4]))
        ,border=NA,col=alpha("firebrick",0.2))
polygon(c(pr.old$logSR,rev(pr.old$logSR))
        ,c(pr.old[,3]-pr.old[,4],rev(pr.old[,3]+pr.old[,4]))
        ,border=NA,col=alpha("skyblue",0.2))
lines(pr.old$logSR,pr.old[,3],lwd=3,col="skyblue")
lines(pr.new$logSR,pr.new[,3],lwd=3,col="firebrick")
lines(c(-1,10), c(0,0), lty=3, col="black")
#--------------------------------
# Figure 5A  - Diversity pre-flood stability by plant history
#--------------------------------
par(mfrow=c(1,2))
dnoNA=d[!(is.na(d$stabBefore)),]
dnoNA$logstabBefore=log(dnoNA$stabBefore)
mm=asreml(logstabBefore~1+logSR+PH+SH+PH:logSR+SH:logSR,random=~Plot/SH,data=dnoNA)
test.asreml(mm)
x.old=tapply(dnoNA[dnoNA$PH=="old",]$logSR,dnoNA[dnoNA$PH=="old",]$Plot,mean,na.rm=T)
y.old=tapply(dnoNA[dnoNA$PH=="old",]$stabBefore,dnoNA[dnoNA$PH=="old",]$Plot,mean,na.rm=T)
x.new=tapply(dnoNA[dnoNA$PH=="new",]$logSR,dnoNA[dnoNA$PH=="new",]$Plot,mean,na.rm=T)
y.new=tapply(dnoNA[dnoNA$PH=="new",]$stabBefore,dnoNA[dnoNA$PH=="new",]$Plot,mean,na.rm=T)
mtemp=lm(terms(logstabBefore~Plot:SH,keep.order=T),data=dnoNA)
dnoNA$res_logstabBefore=residuals(mtemp)
dim(dnoNA)
for (i in c(1,2,4,8)) {
  d1=subset(dnoNA,dnoNA$SR==i)
  dnoNA$adj_logstabBefore[dnoNA$SR==i]=d1$res_logstabBefore+mean(d1$logstabBefore)
}
xi=aggregate(dnoNA$logSR,list(Plot=dnoNA$Plot,PH=dnoNA$PH),mean)
yi=aggregate(exp(dnoNA$adj_logstabBefore),list(Plot=dnoNA$Plot,PH=dnoNA$PH),mean)
plot(jitter(x=xi[xi$PH=="old","x"]-0.05,factor=0.2),y=yi[yi$PH=="old","x"]
     ,xaxt="n",bty="n",tck=0.02,mgp=c(1.5,0.1,0)
     ,ylab=bquote(paste("Pre-flood stability (", italic("CV"["com"])^-1, ")"))
     ,xlab="Planted richness",xlim=c(-0.1,log2(8)+0.1)
     ,pch=21,col=alpha("black",0.6),bg=alpha("skyblue",0.6),lwd=1.5
     ,cex=0.9,tck=0.02,cex.axis=1.2,cex.lab=1.5
     ,ylim=c(1,6),log="y")
points(x=jitter(xi[xi$PH=="new","x"]+0.05,factor=0.2),y=yi[yi$PH=="new","x"]
       ,pch=21,col=alpha("black",0.6),bg=alpha("firebrick",0.6),lwd=1.5,cex=0.9)
axis(side=1,at=log2(c(1,2,4,8)),labels=c(1,2,4,8),tck=0.02,mgp=c(1.5,0.1,0),cex.axis=1.2)
pr.old=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.old,na.rm=T),max(x.old,na.rm=T),l=30)
                                           ,PH="old"))$prediction$pvals
pr.new=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.new,na.rm=T),max(x.new,na.rm=T),l=30)
                                           ,PH="new"))$prediction$pvals
polygon(c(pr.new$logSR,rev(pr.new$logSR))
        ,exp(c(pr.new[,3]-pr.new[,4],rev(pr.new[,3]+pr.new[,4])))
        ,border=NA,col=alpha("firebrick",0.2))
polygon(c(pr.old$logSR,rev(pr.old$logSR))
        ,exp(c(pr.old[,3]-pr.old[,4],rev(pr.old[,3]+pr.old[,4])))
        ,border=NA,col=alpha("skyblue",0.2))
lines(pr.old$logSR,exp(pr.old[,3]),lwd=3,col="skyblue")
lines(pr.new$logSR,exp(pr.new[,3]),lwd=3,col="firebrick")
#--------------------------------
# Figure 5B - Diversity post-flood stability by plant history
#--------------------------------
dnoNA=d[!(is.na(d$stabAfter)),]
dnoNA$logstabAfter=log(dnoNA$stabAfter)
mm=asreml.nvc(logstabAfter~1+logSR+PH+SH+PH:logSR+SH:logSR,random=~Plot/SH,data=dnoNA)
test.asreml(mm)
x.old=tapply(dnoNA[dnoNA$PH=="old",]$logSR,dnoNA[dnoNA$PH=="old",]$Plot,mean,na.rm=T)
y.old=tapply(dnoNA[dnoNA$PH=="old",]$stabAfter,dnoNA[dnoNA$PH=="old",]$Plot,mean,na.rm=T)
x.new=tapply(dnoNA[dnoNA$PH=="new",]$logSR,dnoNA[dnoNA$PH=="new",]$Plot,mean,na.rm=T)
y.new=tapply(dnoNA[dnoNA$PH=="new",]$stabAfter,dnoNA[dnoNA$PH=="new",]$Plot,mean,na.rm=T)
mtemp=lm(terms(logstabAfter~Plot:SH,keep.order=T),data=dnoNA)
dnoNA$res_logstabAfter=residuals(mtemp)
dim(dnoNA)
for (i in c(1,2,4,8)) {
  d1=subset(dnoNA,dnoNA$SR==i)
  dnoNA$adj_logstabAfter[dnoNA$SR==i]=d1$res_logstabAfter+mean(d1$logstabAfter)
}
xi=aggregate(dnoNA$logSR,list(Plot=dnoNA$Plot,PH=dnoNA$PH),mean)
yi=aggregate(exp(dnoNA$adj_logstabAfter),list(Plot=dnoNA$Plot,PH=dnoNA$PH),mean)
plot(jitter(x=xi[xi$PH=="old","x"]-0.05,factor=0.2),y=yi[yi$PH=="old","x"]
     ,xaxt="n",bty="n",tck=0.02,mgp=c(1.5,0.1,0)
     ,ylab=bquote(paste("Post-flood stability (", italic("CV"["com"])^-1, ")"))
     ,xlab="Planted richness",xlim=c(-0.1,log2(8)+0.1)
     ,pch=21,col=alpha("black",0.6),bg=alpha("skyblue",0.6),lwd=1.5
     ,cex=0.9,tck=0.02,cex.axis=1.2,cex.lab=1.5
     ,ylim=c(1,6),log="y")
points(x=jitter(xi[xi$PH=="new","x"]+0.05,factor=0.2),y=yi[yi$PH=="new","x"]
       ,pch=21,col=alpha("black",0.6),bg=alpha("firebrick",0.6),lwd=1.5,cex=0.9)
axis(side=1,at=log2(c(1,2,4,8)),labels=c(1,2,4,8),tck=0.02,mgp=c(1.5,0.1,0),cex.axis=1.2)
pr.old=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.old,na.rm=T),max(x.old,na.rm=T),l=30)
                                           ,PH="old"))$prediction$pvals
pr.new=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.new,na.rm=T),max(x.new,na.rm=T),l=30)
                                           ,PH="new"))$prediction$pvals
polygon(c(pr.new$logSR,rev(pr.new$logSR))
        ,exp(c(pr.new[,3]-pr.new[,4],rev(pr.new[,3]+pr.new[,4])))
        ,border=NA,col=alpha("firebrick",0.2))
polygon(c(pr.old$logSR,rev(pr.old$logSR))
        ,exp(c(pr.old[,3]-pr.old[,4],rev(pr.old[,3]+pr.old[,4])))
        ,border=NA,col=alpha("skyblue",0.2))
lines(pr.old$logSR,exp(pr.old[,3]),lwd=3,col="skyblue")
lines(pr.new$logSR,exp(pr.new[,3]),lwd=3,col="firebrick")
#--------------------------------
# Figure S3A  - Diversity biomass-corrected resistance by plant history
#--------------------------------
par(mfrow=c(1,3))
par(mar=c(5,6,4,2)+0.1,oma=c(1,1,1,1)+0.1)
dnoNA=d[!(is.na(d$resistRaw)),]
#lm1=lm(terms(resistRaw~BMbefore+logSR+SH+SH:logSR+PH+PH:logSR,keep.order=T),data=dnoNA)
#anova(lm1)
#lm2=lm(terms(resistRaw~BMbefore,keep.order=T),data=dnoNA)
#anova(lm2)
#dnoNA$res_resistRaw=residuals(lm2)
#dnoNA$res2_resistRaw=dnoNA$resistRaw-fitted(lm2)
#plot(dnoNA$res_resistRaw,dnoNA$res2_resistRaw)
#dnoNA$adj_resistRaw=dnoNA$res_resistRaw+mean(dnoNA$resistRaw)
#lm3=lm(terms(adj_resistRaw~logSR+SH+SH:logSR+PH+PH:logSR,keep.order=T),data=dnoNA)
#anova(lm3)
#mm1=asreml(resistRaw~BMbefore+logSR+SH+SH:logSR+PH+PH:logSR,keep.order=T,data=dnoNA)
#test.asreml(mm1)
#mm2=asreml(resistRaw~BMbefore,data=dnoNA)
#test.asreml(mm2)
#dnoNA$res_resistRaw=residuals(mm2)
#dnoNA$res2_resistRaw=dnoNA$resistRaw-fitted(mm2)
#plot(dnoNA$res_resistRaw,dnoNA$res2_resistRaw)
#dnoNA$adj_resistRaw=dnoNA$res_resistRaw+mean(dnoNA$resistRaw)
#mm3=asreml(adj_resistRaw~logSR+SH+SH:logSR+PH+PH:logSR,keep.order=T,data=dnoNA)
#test.asreml(mm3)
mm=asreml(resistRaw~1+BMbefore+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH ,data=dnoNA)
test.asreml(mm)
#mmadj=asreml(adj_resistRaw~1+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH ,data=dnoNA)
#test.asreml(mmadj)
x.old=tapply(dnoNA[dnoNA$PH=="old",]$logSR,dnoNA[dnoNA$PH=="old",]$Plot,mean,na.rm=T)
y.old=tapply(dnoNA[dnoNA$PH=="old",]$resistRaw,dnoNA[dnoNA$PH=="old",]$Plot,mean,na.rm=T)
x.new=tapply(dnoNA[dnoNA$PH=="new",]$logSR,dnoNA[dnoNA$PH=="new",]$Plot,mean,na.rm=T)
y.new=tapply(dnoNA[dnoNA$PH=="new",]$resistRaw,dnoNA[dnoNA$PH=="new",]$Plot,mean,na.rm=T)
mtemp=lm(terms(resistRaw~BMbefore,keep.order=T),data=dnoNA)
dnoNA$res_resistRaw=residuals(mtemp)
dnoNA$adj_resistRaw=dnoNA$res_resistRaw+mean(dnoNA$resistRaw)
mtemp1=lm(terms(adj_resistRaw~Plot:SH,keep.order=T),data=dnoNA)
dnoNA$res_resistRaw1=residuals(mtemp1)
for (i in c(1,2,4,8)) {
  d1=subset(dnoNA,dnoNA$SR==i)
  dnoNA$adj_resistRaw1[dnoNA$SR==i]=d1$res_resistRaw1+mean(d1$adj_resistRaw)
}
xi=aggregate(dnoNA$logSR,list(Plot=dnoNA$Plot,PH=dnoNA$PH),mean)
#yi=aggregate(dnoNA$adj_resistRaw,list(Plot=dnoNA$Plot,PH=dnoNA$PH),mean)
yi=aggregate(dnoNA$adj_resistRaw1,list(Plot=dnoNA$Plot,PH=dnoNA$PH),mean)
plot(jitter(x=xi[xi$PH=="old","x"]-0.05,factor=0.2),y=yi[yi$PH=="old","x"]
     ,xaxt="n",bty="n",tck=0.03,mgp=c(3,0.5,0) #default 3,1,0 (1.5,0.1,0)
     ,ylab=bquote(paste("Resistance (g/m"^2, ")"))
     ,xlab="Planted richness",xlim=c(-0.1,log2(8)+0.1)
     ,pch=21,col=alpha("black",0.6),bg=alpha("skyblue",0.6),lwd=1.5
     ,cex=1.5,tck=0.03,cex.axis=1.7,cex.lab=2
#     ,ylim=c(-200,150))
     ,ylim=c(-200,50))
points(x=jitter(xi[xi$PH=="new","x"]+0.05,factor=0.2),y=yi[yi$PH=="new","x"]
       ,pch=21,col=alpha("black",0.6),bg=alpha("firebrick",0.6),lwd=1.5,cex=1.5)
axis(side=1,at=log2(c(1,2,4,8)),labels=c(1,2,4,8),tck=0.03,mgp=c(3,0.5,0),cex.axis=1.7)
pr.old=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.old,na.rm=T),max(x.old,na.rm=T),l=30)
                                           ,PH="old"))$prediction$pvals
pr.new=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.new,na.rm=T),max(x.new,na.rm=T),l=30)
                                           ,PH="new"))$prediction$pvals
#pr.oldadj=predict(mmadj,classify="logSR:PH",list(logSR=seq(min(x.old,na.rm=T),max(x.old,na.rm=T),l=30)
#                                           ,PH="old"))$prediction$pvals
#pr.newadj=predict(mmadj,classify="logSR:PH",list(logSR=seq(min(x.new,na.rm=T),max(x.new,na.rm=T),l=30)
#                                           ,PH="new"))$prediction$pvals
polygon(c(pr.new$logSR,rev(pr.new$logSR))
        ,c(pr.new[,3]-pr.new[,4],rev(pr.new[,3]+pr.new[,4]))
        ,border=NA,col=alpha("firebrick",0.2))
polygon(c(pr.old$logSR,rev(pr.old$logSR))
        ,c(pr.old[,3]-pr.old[,4],rev(pr.old[,3]+pr.old[,4]))
        ,border=NA,col=alpha("skyblue",0.2))
lines(pr.old$logSR,pr.old[,3],lwd=3,col="skyblue")
lines(pr.new$logSR,pr.new[,3],lwd=3,col="firebrick")
lines(c(-1,10), c(0,0), lty=3, col="black")
#--------------------------------
# Figure S3B  - Diversity biomass-corrected recovery by plant history
#--------------------------------
dnoNA=d[!(is.na(d$recovRaw)),]
mm=asreml(recovRaw~1+BMbefore+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH ,data=dnoNA)
test.asreml(mm)
x.old=tapply(dnoNA[dnoNA$PH=="old",]$logSR,dnoNA[dnoNA$PH=="old",]$Plot,mean,na.rm=T)
y.old=tapply(dnoNA[dnoNA$PH=="old",]$recovRaw,dnoNA[dnoNA$PH=="old",]$Plot,mean,na.rm=T)
x.new=tapply(dnoNA[dnoNA$PH=="new",]$logSR,dnoNA[dnoNA$PH=="new",]$Plot,mean,na.rm=T)
y.new=tapply(dnoNA[dnoNA$PH=="new",]$recovRaw,dnoNA[dnoNA$PH=="new",]$Plot,mean,na.rm=T)
mtemp=lm(terms(recovRaw~BMbefore,keep.order=T),data=dnoNA)
dnoNA$res_recovRaw=residuals(mtemp)
dnoNA$adj_recovRaw=dnoNA$res_recovRaw+mean(dnoNA$recovRaw)
mtemp1=lm(terms(adj_recovRaw~Plot:SH,keep.order=T),data=dnoNA)
dnoNA$res_recovRaw1=residuals(mtemp1)
for (i in c(1,2,4,8)) {
  d1=subset(dnoNA,dnoNA$SR==i)
  dnoNA$adj_recovRaw1[dnoNA$SR==i]=d1$res_recovRaw1+mean(d1$adj_recovRaw)
}
xi=aggregate(dnoNA$logSR,list(Plot=dnoNA$Plot,PH=dnoNA$PH),mean)
#yi=aggregate(dnoNA$adj_recovRaw,list(Plot=dnoNA$Plot,PH=dnoNA$PH),mean)
yi=aggregate(dnoNA$adj_recovRaw1,list(Plot=dnoNA$Plot,PH=dnoNA$PH),mean)
plot(jitter(x=xi[xi$PH=="old","x"]-0.05,factor=0.2),y=yi[yi$PH=="old","x"]
     ,xaxt="n",bty="n",tck=0.03,mgp=c(3,0.5,0)
     ,ylab=bquote(paste("Recovery (g/m"^2, ")"))
     ,xlab="Planted richness",xlim=c(-0.1,log2(8)+0.1)
     ,pch=21,col=alpha("black",0.6),bg=alpha("skyblue",0.6),lwd=1.5
     ,cex=1.5,tck=0.03,cex.axis=1.7,cex.lab=2
     ,ylim=c(-100,350))
points(x=jitter(xi[xi$PH=="new","x"]+0.05,factor=0.2),y=yi[yi$PH=="new","x"]
       ,pch=21,col=alpha("black",0.6),bg=alpha("firebrick",0.6),lwd=1.5,cex=1.5)
axis(side=1,at=log2(c(1,2,4,8)),labels=c(1,2,4,8),tck=0.03,mgp=c(3,0.5,0),cex.axis=1.7)
pr.old=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.old,na.rm=T),max(x.old,na.rm=T),l=30)
                                           ,PH="old"))$prediction$pvals
pr.new=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.new,na.rm=T),max(x.new,na.rm=T),l=30)
                                           ,PH="new"))$prediction$pvals
polygon(c(pr.new$logSR,rev(pr.new$logSR))
        ,c(pr.new[,3]-pr.new[,4],rev(pr.new[,3]+pr.new[,4]))
        ,border=NA,col=alpha("firebrick",0.2))
polygon(c(pr.old$logSR,rev(pr.old$logSR))
        ,c(pr.old[,3]-pr.old[,4],rev(pr.old[,3]+pr.old[,4]))
        ,border=NA,col=alpha("skyblue",0.2))
lines(pr.old$logSR,pr.old[,3],lwd=3,col="skyblue")
lines(pr.new$logSR,pr.new[,3],lwd=3,col="firebrick")
lines(c(-1,10), c(0,0), lty=3, col="black")
#--------------------------------
# Figure S3C  - Diversity biomass-corrected resilience by plant history
#--------------------------------
dnoNA=d[!(is.na(d$resilRaw)),]
mm=asreml(resilRaw~1+BMbefore+logSR+SH+PH+SH:logSR+PH:logSR,random=~Plot/SH ,data=dnoNA)
test.asreml(mm)
x.old=tapply(dnoNA[dnoNA$PH=="old",]$logSR,dnoNA[dnoNA$PH=="old",]$Plot,mean,na.rm=T)
y.old=tapply(dnoNA[dnoNA$PH=="old",]$resilRaw,dnoNA[dnoNA$PH=="old",]$Plot,mean,na.rm=T)
x.new=tapply(dnoNA[dnoNA$PH=="new",]$logSR,dnoNA[dnoNA$PH=="new",]$Plot,mean,na.rm=T)
y.new=tapply(dnoNA[dnoNA$PH=="new",]$resilRaw,dnoNA[dnoNA$PH=="new",]$Plot,mean,na.rm=T)
mtemp=lm(terms(resilRaw~BMbefore,keep.order=T),data=dnoNA)
dnoNA$res_resilRaw=residuals(mtemp)
dnoNA$adj_resilRaw=dnoNA$res_resilRaw+mean(dnoNA$resilRaw)
mtemp1=lm(terms(adj_resilRaw~Plot:SH,keep.order=T),data=dnoNA)
dnoNA$res_resilRaw1=residuals(mtemp1)
for (i in c(1,2,4,8)) {
  d1=subset(dnoNA,dnoNA$SR==i)
  dnoNA$adj_resilRaw1[dnoNA$SR==i]=d1$res_resilRaw1+mean(d1$adj_resilRaw)
}
xi=aggregate(dnoNA$logSR,list(Plot=dnoNA$Plot,PH=dnoNA$PH),mean)
#yi=aggregate(dnoNA$adj_resilRaw,list(Plot=dnoNA$Plot,PH=dnoNA$PH),mean)
yi=aggregate(dnoNA$adj_resilRaw1,list(Plot=dnoNA$Plot,PH=dnoNA$PH),mean)
plot(jitter(x=xi[xi$PH=="old","x"]-0.05,factor=0.2),y=yi[yi$PH=="old","x"]
     ,xaxt="n",bty="n",tck=0.03,mgp=c(3,0.5,0)
     ,ylab=bquote(paste("Resilience (g/m"^2, ")"))
     ,xlab="Planted richness",xlim=c(-0.1,log2(8)+0.1)
     ,pch=21,col=alpha("black",0.6),bg=alpha("skyblue",0.6),lwd=1.5
     ,cex=1.5,tck=0.03,cex.axis=1.7,cex.lab=2
#     ,ylim=c(-500,500))
     ,ylim=c(-400,300))
points(x=jitter(xi[xi$PH=="new","x"]+0.05,factor=0.2),y=yi[yi$PH=="new","x"]
       ,pch=21,col=alpha("black",0.6),bg=alpha("firebrick",0.6),lwd=1.5,cex=1.5)
axis(side=1,at=log2(c(1,2,4,8)),labels=c(1,2,4,8),tck=0.03,mgp=c(3,0.5,0),cex.axis=1.7)
pr.old=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.old,na.rm=T),max(x.old,na.rm=T),l=30)
                                           ,PH="old"))$prediction$pvals
pr.new=predict(mm,classify="logSR:PH",list(logSR=seq(min(x.new,na.rm=T),max(x.new,na.rm=T),l=30)
                                           ,PH="new"))$prediction$pvals
polygon(c(pr.new$logSR,rev(pr.new$logSR))
        ,c(pr.new[,3]-pr.new[,4],rev(pr.new[,3]+pr.new[,4]))
        ,border=NA,col=alpha("firebrick",0.2))
polygon(c(pr.old$logSR,rev(pr.old$logSR))
        ,c(pr.old[,3]-pr.old[,4],rev(pr.old[,3]+pr.old[,4]))
        ,border=NA,col=alpha("skyblue",0.2))
lines(pr.old$logSR,pr.old[,3],lwd=3,col="skyblue")
lines(pr.new$logSR,pr.new[,3],lwd=3,col="firebrick")
lines(c(-1,10), c(0,0), lty=3, col="black")

#--------------------------------
# Figure S5 - Effect Sizes Barplot
#--------------------------------
setwd("YourDirectory")

d<- read.csv("Data_Fig.S5.csv", sep=";")
str(d)
d$panel <- NA
#add a column to define panel A or B
d[d$response=="Stability" | d$response=="Population variance" | d$response=="Asynchrony",]$panel <- "A"
d[d$response=="Resistance" | d$response=="Recovery" | d$response=="Resilience",]$panel <- "B"

library(ggplot2)
library(cowplot)
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#000000", "#0072B2", "#D55E00", "#CC79A7")

d$Factor <- factor(d$Factor, levels=c("Log richness (Rlog)", "Soil treatment (SH)","Co-occurrence history (CH)", "SH x Rlog", "CH x Rlog"))

d$response <- factor(d$response, levels=c("Stability", "Population variance","Asynchrony", "Resistance", "Resilience", "Recovery"))

p1 <- ggplot() + geom_bar(aes(y = SS, x = response, fill = Factor), data = d[d$panel=="A",],
                          stat="identity") + 
  labs(x="Response variable", y="% SS") +
  theme(legend.position = "none")+
  scale_fill_manual(values=cbPalette) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


p2 <- ggplot() + geom_bar(aes(y = SS, x = response, fill = Factor), data = d[d$panel=="B",],
                          stat="identity")+ 
  labs(x="Response variable", y="% SS") +
  theme(legend.title = element_blank())+
  scale_fill_manual(values=cbPalette) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

quartz()
plot_grid(p1,p2, labels = c('A', 'B'), label_size = 10, ncol=2, nrow=1)
