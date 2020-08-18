#Resilience and resistance towards the flood: compositional changes

rm(list=ls())
ls()

Sys.setenv(ASREML_LICENSE_FILE="/Applications/asreml3/bin/asreml.lic")

library(vegan) # Package for community data
library(pascal) # for ASREML
library(asreml)
library(asremlPlus) # nice functions for extracting data
library(colorRamps) # for making nice colors
library(scales) # for plotting transparent colors

##### LOADING DATA

setwd("your directory")
d = read.csv("vanMoorsel_EtAl_Data_Stability metrics.csv", sep=",")


head(d, n=20)
d$SR
# making factors and contrasts:
d$Plot=factor(d$Plot);nlevels(d$Plot)
d$SH=factor(d$SH);nlevels(d$SH)
d$SH_c=factor(ifelse(d$SH=='Str',1,2)) # contrast
d$PH=factor(d$PH);nlevels(d$PH)
d$facSR=factor(d$SR);nlevels(d$facSR)
d$logSR=log2(d$SR) # log species richness linear
d$logSR=d$logSR-mean(d$logSR,na.rm=T) # center
d$SR8=factor(as.numeric(d$SR==8))
d$SR4=factor(as.numeric(d$SR==4))
d$SR2=factor(as.numeric(d$SR==2))
head(d)
d$wei=as.numeric(d$Plot!="B1A12"|d$SH!="Org"|d$PH!="old")
d=d[d$wei>0,] # excluding outlier half-quadrat

d$wei=as.numeric(d$Plot!="B1A12"|d$SH!="Org")
d=d[d$wei>0,] # excluding outlier quadrat (both plant histories!!)
dim(d)
d1=d[d$SH!="Org",] # excluding SH Org
dim(d1)

#--------------------------------
# ANOVA stability:
par(mfrow=c(1,2))
hist(d$Stab,n=100)
tapply(d$Stab,list(d$facSR,d$PH),mean)
tapply(d$Stab,list(d$facSR,d$SH),mean)
tapply(d$Stab,list(d$PH,d$facSR,d$SH),mean)
plot(d$Stab~d$logSR)
m=aov(terms(Stab~1
            +(logSR+facSR)+Plot
            +(SH_c+SH)
            +(SH_c+SH):(logSR+facSR)+Plot:SH
            +PH
            +PH:(logSR+facSR)
            +PH:(SH_c+SH)
            +PH:(SH_c+SH):(logSR+facSR)
            ,keep.order=T)
      ,data=d)
summary.aov(m,intercept=T)
quartz()
par(mfrow=c(2,2))
plot(m)
mm=asreml(Stab~1
          +facSR
          +SH
          +PH
          +SH:facSR
          +PH:facSR
          +PH:SH
          +facSR:SH:PH
          ,random =~Plot/SH 
          ,data=d)
test.asreml(mm)
plot(mm)
pvStab=predict(mm,classify="SH")
pvStab$predictions
pvStab=predict(mm,classify="facSR:PH")
pvStab$predictions
mm=asreml(Stab~1
          +facSR
          +SH
          +PH
          +SH:facSR
          +PH:facSR
          +PH:SH
          +facSR:SH:PH
          ,random =~Plot/SH 
          ,data=d1) #this is data excluding the Org soil history! --> effects much stronger
test.asreml(mm)
mm=asreml(Stab~1
          +logSR
          +SH
          +PH
          +SH:logSR
          +PH:logSR
          ,random =~Plot/SH 
          ,data=d) #d1)
test.asreml(mm)
pvStab=predict(mm,classify="logSR:PH",levels=list(logSR=0))
pvStab$predictions # Marginally selected communities are more stable.

mm=asreml(Stab~1
          +Sync
          +logSR
          +SH
          +PH
          +Sync:logSR
          +Sync:PH
          ,random =~Plot/SH 
          ,data=d)
test.asreml(mm)
mm=asreml(Stab~1
          +Sync
#          +logSR
#          +SH
          +PH
          +Sync:PH
          ,random =~Plot/SH 
          ,data=d) #d1)
test.asreml(mm) # ************************************
coefficients(mm)$fixed # Synchrony reduces stability; selected communities are more stable.
pvStab=predict(mm,classify="Sync:PH",levels=list(Sync=0))
pvStab$predictions # Synchrony reduces stability; selected communities are more stable.
pvStab=predict(mm,classify="Sync:PH",levels=list(Sync=1))
pvStab$predictions # Synchrony reduces stability; selected communities are more stable.
d$logpopCV=log(d$popCV)
mm=asreml(Stab~1
          +logpopCV
          +logSR
          +SH
          +PH
          +logpopCV:logSR
          +logpopCV:PH
          ,random =~Plot/SH 
          ,data=d) #d1)
test.asreml(mm)  # ************************************
coefficients(mm)$fixed # Population CV reduces(!) stability.
pvStab=predict(mm,classify="logpopCV:PH",levels=list(logpopCV=-0.5))
pvStab$predictions # At low population CV selected communities are less stable.
pvStab=predict(mm,classify="logpopCV:PH",levels=list(logpopCV=0.5))
pvStab$predictions # But at high population CV selected communities are more stable.

#--------------------------------
# ANOVA cv:
mm=asreml(log(CV)~1
          +logSR
          +SH
          +PH
          +SH:logSR
          +PH:logSR
          +PH:SH
          +logSR:SH:PH
          ,random =~Plot/SH 
          ,data=d)
test.asreml(mm)
plot(mm)
coefficients(mm)$fixed # Selected communities appear to have higher community-level CV.
pvCV=predict(mm,classify="SH")
pvCV$predictions
pvCV=predict(mm,classify="PH")
pvCV$predictions # Selected communities actually have lower community-level CV.
mm=asreml(log(CV)~1
          +logSR
          +SH
          +PH
          +SH:logSR
          +PH:logSR
          +PH:SH
          +logSR:SH:PH
          ,random =~Plot/SH 
          ,data=d1)
test.asreml(mm)

mm=asreml(log(CV)~1
          +Sync
          +logSR
          +SH
          +PH
          +Sync:logSR
          +Sync:PH
          ,random =~Plot/SH 
          ,data=d) #d1)
test.asreml(mm) # ************************************
pvCV=predict(mm,classify="Sync:PH",levels=list(Sync=-0))
pvCV$predictions # At low synchrony selected communities have similar community-level CV.
pvCV=predict(mm,classify="Sync:PH",levels=list(Sync=1))
pvCV$predictions # But at maximum synchrony they have lower community-level CVs than na?ve communities.
mm=asreml(log(CV)~1
          +Sync
          +PH
          +Sync:PH
          ,random =~Plot/SH 
          ,data=d) #d1)
test.asreml(mm) # ************************************ here only marginal:
pvCV=predict(mm,classify="Sync:PH",levels=list(Sync=-0))
pvCV$predictions # At low synchrony selected communities have similar community-level CV.
pvCV=predict(mm,classify="Sync:PH",levels=list(Sync=1))
pvCV$predictions # But at maximum synchrony they have lower community-level CVs than na?ve communities.

mm=asreml(log(CV)~1
          +logpopCV
          +PH
          +logpopCV:PH
          ,random =~Plot/SH 
          ,data=d) #d1)
test.asreml(mm)
mm=asreml(log(CV)~1
          +logpopCV
          +logSR
          +SH
          ,random =~Plot/SH 
          ,data=d) #d1)
test.asreml(mm)
coefficients(mm)$fixed

#--------------------------------
# ANOVA sync:
tapply(d$Sync,list(d$facSR,d$PH),mean)
tapply(d$Sync,list(d$facSR,d$SH),mean)
tapply(d$Sync,list(d$PH,d$facSR,d$SH),mean)
par(mfrow=c(1,2))
plot(d$Sync~d$logSR)
plot(log(d$CV)~d$Sync)
mm=asreml(Sync~1
#          +Stab
          +logSR
          +SH
          +PH
          +SH:logSR
          +PH:logSR
          +PH:SH
          +logSR:SH:PH
          ,random =~Plot/SH 
          ,data=d)
test.asreml(mm) # Plany history does not influence synchrony.
plot(mm)
coefficients(mm)$fixed
pvSync=predict(mm,classify="SH")
pvSync$predictions # But communities are much less synchronized on original soil.
pvSync=predict(mm,classify="logSR:PH:SH",levels=list(logSR=-1))
pvSync$predictions
pvSync=predict(mm,classify="logSR:PH:SH",levels=list(logSR=1))
pvSync$predictions

#--------------------------------
# ANOVA population cv:
mm=asreml(logpopCV~1
          +logSR
          +SH
          +PH
#          +SH:logSR
          +PH:logSR
          ,random =~Plot/SH 
          ,data=d)
test.asreml(mm) # ************************************
plot(mm)
coefficients(mm)$fixed # The population CV is generally lower in selected communities.
pvpopCV=predict(mm,classify="logSR:PH",levels=list(logSR=-1))
pvpopCV$predictions # This is due to a much reduced polulation CV in selected communities at low diversity.
pvpopCV=predict(mm,classify="logSR:PH",levels=list(logSR=1))
pvpopCV$predictions # At high diversity the trend is actually reversed.

mm=asreml(logpopCV~1
          +Sync
          +logSR
          +SH
          +PH
          +Sync:logSR
          +Sync:PH
          ,random =~Plot/SH 
          ,data=d) #d1) # (Much nicer with d1.)
test.asreml(mm)
mm=asreml(logpopCV~1
          +Sync
          +logSR
          +SH
          +PH
          +Sync:PH
          ,random =~Plot/SH 
          ,data=d) #d1)
test.asreml(mm)
coefficients(mm)$fixed
mm=asreml(logpopCV~1
          +Sync
          +PH
          +Sync:PH
          ,random =~Plot/SH 
          ,data=d) #d1)
test.asreml(mm) # ************************************
coefficients(mm)$fixed # The population CV is generally lower in selected communities.
pvpopCV=predict(mm,classify="Sync:PH",levels=list(Sync=0))
pvpopCV$predictions
pvpopCV=predict(mm,classify="Sync:PH",levels=list(Sync=1))
pvpopCV$predictions # This is especially true at maximum synchrony.
# In other words, synchrony is more strongly negatively related to population CV
# in selected than in na?ve communities:
quartz()
par(mfrow=c(1,3))
plot(d$logpopCV~d$Sync)
plot(d[d$PH=="old",]$logpopCV~d[d$PH=="old",]$Sync
     ,xlab="Synchrony"
     ,ylab=expression("log(population-level CV)"))
lm1=lm(d[d$PH=="old",]$logpopCV~d[d$PH=="old",]$Sync) 
abline(lm1)
plot(d[d$PH=="new",]$logpopCV~d[d$PH=="new",]$Sync
     ,xlab="Synchrony"
     ,ylab=expression("log(population-level CV)"))
lm2=lm(d[d$PH=="new",]$logpopCV~d[d$PH=="new",]$Sync) 
abline(lm2)



#Figure A describeds the relationship between synchrony and community-level variance (for old/new)
#Figure B describes the relationship between synchrony and population-level variance (for old/new)
#Figure C describes the relationship between population-level variance and diversity (for old/new)
#Figure D describes the relationship between community-level variance and diversity (for old/new)

#A

quartz()
plot(d[d$PH=="old",]$Sync, d[d$PH=="old",]$CV, pch=21, xlab = "Synchrony", ylab="Community-level CV")
points(d[d$PH=="new",]$Sync, d[d$PH=="new",]$CV, pch=19)
abline(lm(d[d$PH=="new",]$CV~d[d$PH=="new",]$Sync), lty=1, lwd=2)
abline(lm(d[d$PH=="old",]$CV~d[d$PH=="old",]$Sync), lty=2, lwd=2)
legend("top", bty="n", pch=c(19, 21), lty=c(1:2), c("Naïve communities", "Selected communities"))

#B

quartz()
plot(d[d$PH=="old",]$Sync, d[d$PH=="old",]$logpopCV, pch=21, xlab = "Synchrony", ylab="Log Population-level CV")
points(d[d$PH=="new",]$Sync, d[d$PH=="new",]$logpopCV, pch=19)
abline(lm(d[d$PH=="new",]$logpopCV~d[d$PH=="new",]$Sync), lty=1, lwd=2)
abline(lm(d[d$PH=="old",]$logpopCV~d[d$PH=="old",]$Sync), lty=2, lwd=2)
legend("top", bty="n", pch=c(19, 21), lty=c(1:2), c("Naïve communities", "Selected communities"))


#C
d$logSR
quartz()
plot(d[d$PH=="old",]$logSR, d[d$PH=="old",]$logpopCV, pch=21, xlab = "Log SR", ylab="Log Population-level CV")
points(d[d$PH=="new",]$logSR, d[d$PH=="new",]$logpopCV, pch=19)
abline(lm(d[d$PH=="new",]$logpopCV~d[d$PH=="new",]$logSR), lty=1, lwd=2)
abline(lm(d[d$PH=="old",]$logpopCV~d[d$PH=="old",]$logSR), lty=2, lwd=2)
legend("top", bty="n", pch=c(19, 21), lty=c(1:2), c("Naïve communities", "Selected communities"))

#D
d$logSR
quartz()
plot(d[d$PH=="old",]$logSR, d[d$PH=="old",]$CV, pch=21, xlab = "Log SR", ylab="Community-level CV")
points(d[d$PH=="new",]$logSR, d[d$PH=="new",]$CV, pch=19)
abline(lm(d[d$PH=="new",]$CV~d[d$PH=="new",]$logSR), lty=1, lwd=2)
abline(lm(d[d$PH=="old",]$CV~d[d$PH=="old",]$logSR), lty=2, lwd=2)
legend("top", bty="n", pch=c(19, 21), lty=c(1:2), c("Naïve communities", "Selected communities"))


#For the paper
#A: x-axis Synchrony, y-axis: Stability

quartz()
par(mfrow=c(1,4), oma=c(0,0,0.1,0.1))
plot(d[d$PH=="old",]$Sync, d[d$PH=="old",]$Stab, pch=19, xlab = "Synchrony", ylab="Stability", xlim=c(0,1))
points(d[d$PH=="new",]$Sync, d[d$PH=="new",]$Stab, pch=21)
abline(lm(d[d$PH=="new",]$Stab~d[d$PH=="new",]$Sync), lty=2, lwd=2)
abline(lm(d[d$PH=="old",]$Stab~d[d$PH=="old",]$Sync), lty=1, lwd=2)
legend("topright", bty="n", pch=c(19, 21), lty=c(1:2), c("Selected communities", "Naïve communities"))

#B: x-axis Synchrony, y-axis: CV

plot(d[d$PH=="old",]$Sync, d[d$PH=="old",]$CV, pch=19, xlab = "Synchrony", ylab="Community-level CV", xlim=c(0,1))
points(d[d$PH=="new",]$Sync, d[d$PH=="new",]$CV, pch=21)
abline(lm(d[d$PH=="new",]$CV~d[d$PH=="new",]$Sync), lty=2, lwd=2)
abline(lm(d[d$PH=="old",]$CV~d[d$PH=="old",]$Sync), lty=1, lwd=2)


#C: x-axis Synchrony, yaxis: popCV

plot(d[d$PH=="old",]$Sync, d[d$PH=="old",]$logpopCV, pch=19, xlab = "Synchrony", ylab="Log Population-level CV", xlim=c(0,1))
points(d[d$PH=="new",]$Sync, d[d$PH=="new",]$logpopCV, pch=21)
abline(lm(d[d$PH=="new",]$logpopCV~d[d$PH=="new",]$Sync), lty=2, lwd=2)
abline(lm(d[d$PH=="old",]$logpopCV~d[d$PH=="old",]$Sync), lty=1, lwd=2)

#D x.axis: popCV, yaxis: Stability

plot(d[d$PH=="old",]$logpopCV, d[d$PH=="old",]$Stab, pch=19, xlab = "Log Population-level CV", ylab="Stability")
points(d[d$PH=="new",]$logpopCV, d[d$PH=="new",]$Stab, pch=21)
abline(lm(d[d$PH=="new",]$Stab~d[d$PH=="new",]$logpopCV), lty=2, lwd=2)
abline(lm(d[d$PH=="old",]$Stab~d[d$PH=="old",]$logpopCV), lty=1, lwd=2)






#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ANOVA turnover before flood:
m1=asreml(turnBefore~1
          +logSR
         # +SH
          +PH
          #+SH:logSR
          +PH:logSR
         # +PH:SH
         # +logSR:SH:PH
          ,random =~Plot/SH 
          ,data=d)
test.asreml(m1)
plot(m1)

#--------------------------------
# ANOVA turnover after flood:
mm=asreml(turnAfter~1
          +logSR
          +SH
          +PH
          +SH:logSR
          +PH:logSR
          +PH:SH
          +logSR:SH:PH
          ,random =~Plot/SH 
          ,data=d)
test.asreml(mm)
plot(mm)
m2=asreml(turnAfter~1
          +logSR
          +PH
          +PH:logSR
          ,random =~Plot/SH 
          ,data=d)
test.asreml(m2)
coefficients(mm)$fixed # Selected communities have higher turnover after flood at high SR.

#--------------------------------
# ANOVA turnover before to after flood:
m3=asreml(turnBefAft~1
          +logSR
          #+SH
          +PH
          #+SH:logSR
          +PH:logSR
         # +PH:SH
          #+logSR:SH:PH
          ,random =~Plot/SH 
          ,data=d)
test.asreml(m3)
plot(m3)

#--------------------------------
head(d)
# mean of ANOVA turnover before and after flood:
d$meanTurn = ((d$turnAfter+d$turnBefore)/2)
str(d$meanTurn)
str(d$turnAfter)
m4=asreml(meanTurn ~1
          +logSR
          #+SH
          +PH
          #+SH:logSR
          +PH:logSR
          # +PH:SH
          #+logSR:SH:PH
          ,random =~Plot/SH 
          ,data=d)
test.asreml(m4)
plot(m4)
#-------------------------------------------------------------------------------------------------
#Barplot turnover figure
#raw data? or from predicted values?
quartz()
par(mfrow=c(2,6), las=1)
par(xpd=T)

m1=asreml(meanTurn~1
          +SR
          + PH
          +PH:SR
          ,random =~Plot/SH 
          ,data=d)
test.asreml(m1)
m2=asreml(turnBefore~1
          +SR
          + PH
          +PH:SR
          ,random =~Plot/SH 
          ,data=d)
test.asreml(m2)

m3=asreml(turnAfter~1
          +SR
          + PH
          +PH:SR
          ,random =~Plot/SH 
          ,data=d)
test.asreml(m3)


pred1 = predict(m1, classify="PH:SR")$predictions$pvals
pred2 = predict(m2, classify="PH:SR")$predictions$pvals
pred3 = predict(m3, classify="PH:SR")$predictions$pvals
#pred4 = predict(m4, classify="Plant:fSR")$predictions$pvals
#-------------------------------------------------------------------------------------------------#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------#-------------------------------------------------------------------------------------------------
head(d)
str(d)
aft <- aggr(d, factors=c("SR", "PH"), newcols=c("mean.aft=mean((turnAfter),na.rm=T)", "se.aft=se((turnAfter),na.rm=T)"))
bef <- aggr(d, factors=c("SR", "PH"), newcols=c("mean.bef=mean((turnBefore),na.rm=T)", "se.bef=se((turnBefore),na.rm=T)"))
befaft <- aggr(d, factors=c("SR", "PH"), newcols=c("mean.befaft=mean((turnBefAft),na.rm=T)", "se.befaft=se((turnBefAft),na.rm=T)"))

all <- aggr(d, factors=c("SR", "PH"), newcols=c("mean.all=mean((meanTurn),na.rm=T)", "se.all=se((meanTurn),na.rm=T)"))

quartz()
par(mfrow=c(1,3), las=1)
par(xpd=T)
warnings()


a3 <- barplot(bef$mean.bef, main="Pre-flood",
              xlab="", col=c("darkgrey","white"),
              ylab="Turnover (Bray-Curtis)",
              mgp=c(1,0.5,0),
              tck=-0.04,
              cex.axis=0.8,
              cex.main=0.8,
              space=c(0.2,0.2,0.5,0.2, 0.5, 0.2),
              ylim=c(0,1))

arrows(a3, bef$mean.bef - bef$se.bef, a3,
       bef$mean.bef + bef$se.bef , lwd = 1.5, angle = 90, code = 3, length = 0.02)
legend("topleft", c("Selected communities", "Naïve communities"), fill = c("darkgrey", "white"),
       border = "black", bty="n", cex=1)


a3 <- barplot(aft$mean.aft, main="Post-flood",
              xlab="", col=c("darkgrey","white"),
              ylab="Turnover (Bray-Curtis)",
              mgp=c(1,0.5,0),
              tck=-0.04,
              cex.axis=0.8,
              cex.main=0.8,
              space=c(0.2,0.2,0.5,0.2, 0.5, 0.2),
              ylim=c(0,1))

arrows(a3, aft$mean.aft - aft$se.aft, a3,
       aft$mean.aft + aft$se.aft , lwd = 1.5, angle = 90, code = 3, length = 0.02)



a1 <- barplot(befaft$mean.befaft, main="Post-flood vs. pre-flood",
              xlab="", col=c("darkgrey","white"),
              ylab="Turnover (Bray-Curtis)",
              mgp=c(1,0.5,0),
              tck=-0.04,
              cex.axis=0.8,
              cex.main=0.8,
              space=c(0.2,0.2,0.5,0.2, 0.5, 0.2),
              ylim=c(0,1))

arrows(a1, befaft$mean.befaft - befaft$se.befaft, a1,
       befaft$mean.befaft + befaft$se.befaft , lwd = 1.5, angle = 90, code = 3, length = 0.02)




#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ANOVA resistance to flood:
mm=asreml(resistRaw~1
          +logSR
          +SH
          +PH
          +SH:logSR
          +PH:logSR
          +PH:SH
          +logSR:SH:PH
          ,random =~Plot/SH 
          ,data=d)
test.asreml(mm)
plot(mm)
d$wei=as.numeric(d$Plot!="B1A12"|d$SH!="Org"|d$PH!="old")
mm=asreml(resistRaw~1
          +logSR
          +SH
          +PH
          +SH:logSR
          +PH:logSR
#          +PH:SH
#          +logSR:SH:PH
          ,random =~Plot/SH
          ,weights=wei
          ,data=d)
test.asreml(mm)
par(mfrow=c(1,1))
hist(residuals(mm),n=100)
coefficients(mm)$fixed # Selected communities have lower absolute resistance to flood.

mm=asreml(resistRel~1
          +logSR
          +SH
          +PH
          +SH:logSR
          +PH:logSR
          +PH:SH
          +logSR:SH:PH
          ,random =~Plot/SH 
          ,weights=wei
          ,data=d)
test.asreml(mm)
plot(mm)

#--------------------------------
# ANOVA recovery from flood:
mm=asreml(recovRaw~1
          +logSR
          +SH
          +PH
          +SH:logSR
          +PH:logSR
          +PH:SH
          +logSR:SH:PH
          ,random =~Plot/SH 
          ,data=d)
test.asreml(mm)
plot(mm)
mm=asreml(recovRaw~1
          +logSR
          +PH
          ,random =~Plot/SH
          ,weights=wei
          ,data=d)
test.asreml(mm)
par(mfrow=c(1,1))
hist(residuals(mm),n=100)
coefficients(mm)$fixed # Selected communities have higher absolute recovery after flood.

mm=asreml(recovRel~1
          +logSR
          +SH
          +PH
          +SH:logSR
          +PH:logSR
          +PH:SH
          +logSR:SH:PH
          ,random =~Plot/SH 
          ,data=d)
test.asreml(mm)
plot(mm)
mm=asreml(recovRel~1
          +logSR
          +PH
          ,random =~Plot/SH
          ,weights=wei
          ,data=d)
test.asreml(mm)
par(mfrow=c(1,1))
hist(residuals(mm),n=100)
coefficients(mm)$fixed # Selected communities have higher relative recovery after flood.

#--------------------------------
# ANOVA resilience during flood:
mm=asreml(resilRaw~1
          +logSR
          +SH
          +PH
          +SH:logSR
          +PH:logSR
          +PH:SH
          +logSR:SH:PH
          ,random =~Plot/SH 
          ,data=d)
test.asreml(mm)
plot(mm)
mm=asreml(resilRaw~1
          +logSR
          +SH
          +PH
          +SH:logSR
          ,random =~Plot/SH
#         ,weights=wei
          ,data=d)
test.asreml(mm)
par(mfrow=c(1,1))
hist(residuals(mm),n=100)
coefficients(mm)$fixed # Selected communities have higher absolute recovery after flood.

mm=asreml(resilRel~1
          +logSR
          +SH
          +PH
          +SH:logSR
          +PH:logSR
          +PH:SH
          +logSR:SH:PH
          ,random =~Plot/SH 
          ,data=d)
test.asreml(mm)
plot(mm)
mm=asreml(resilRel~1
          +logSR
          +SH
          +PH
          +SH:logSR
          +PH:logSR
          ,random =~Plot/SH
          ,weights=wei
          ,data=d)
test.asreml(mm)
par(mfrow=c(1,1))
hist(residuals(mm),n=100)
coefficients(mm)$fixed

