############################
# V1_Nitrogen
# Isotopomer Data Analysis
# Quentin Schorpp
# 23.05.2015
###########################

# N2O und CO2 on day 3 and day 21 (termin 2 and 8) ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# create data frame
isom.dates <- cbind(n2o.index, n2o.bod[,c(2,8)], co2.bod[,c(2,8)])
colnames(isom.dates)[6:9] <- c("n2o.day3","n2o.day21","co2.day3","co2.day21")

# ANOVA Analysis

Lt <- rep(c(1,1,0,0),each=4,2) 
Lt <- c(Lt, rep(c(1,1,0,0),each=3, 2))
isom.dates$Lt <- Lt

Fc <- rep(c(0,1,1,0),each=4,2) 
Fc <- c(Fc, rep(c(0,1,1,0),each=3, 2))
isom.dates$Fc <- Fc

isom.dates$int <- with(isom.dates, interaction(soil, treat))

isom.dates.exp2 <- isom.dates[isom.dates$experiment=="exp2",]
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Analyses ####

### N2O ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n2o.day3.aov <- aov(n2o.day3 ~ soil*Lt*Fc, isom.dates.exp2); summary(n2o.day3.aov)
n2o.day3.aov <- aov(n2o.day3 ~ soil*Lt + soil*Fc, isom.dates.exp2); summary(n2o.day3.aov)
n2o.day3.aov <- update( n2o.day3.aov,log(.)~ soil*treat, isom.dates.exp2); summary(n2o.day3.aov)

par(mfrow=c(2,2))
plot(n2o.day3.aov)

# Post Hoc Tukey Test

n2o.day3.aov2 <- update(n2o.day3.aov, .~ int, isom.dates.exp2)

dn.tuk <- glht(n2o.day3.aov2, linfct = mcp(int = "Tukey"))
summary(dn.tuk)          # standard display
dn.tuk.cld <- cld(dn.tuk)   # letter-based display
par(mai=c(1,1,1.1,0.2),  # mai specifies margin size in inches
    mfrow=c(1,1),
    mgp=c(2,1,0),
    cex=0.8) 
plot(dn.tuk.cld, cex=0.5, las=2, col="grey", xaxt="n")
axis(1, at=gas15.cum$int,labels=FALSE)
text(gas15.cum$int, labels=gas15.cum$int, par("usr")[3], adj=c(1.2,1.2), xpd=TRUE, srt=45, cex=0.8)

          
n2o.day3 <- with(isom.dates.exp2, list( Mean=round(tapply(n2o.day3, list(treat,soil), mean),1),
                                         SE=round(tapply(n2o.day3, list(treat,soil), se),1)))
data.frame(n2o.day3)[c(1,3,2,4)];rm(n2o.day3)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n2o.day21.aov <- aov(n2o.day21 ~ soil*Lt*Fc, isom.dates.exp2); summary(n2o.day21.aov)
n2o.day21.aov <- aov(n2o.day21 ~ soil*Lt, isom.dates.exp2); summary(n2o.day21.aov)
n2o.day21.aov <- update( n2o.day21.aov,log(.)~ soil*treat, isom.dates.exp2); summary(n2o.day21.aov)

par(mfrow=c(2,2))
plot(n2o.day21.aov)


n2o.day21.aov2 <- update( n2o.day21.aov,.~ int, isom.dates.exp2); summary(n2o.day21.aov)
dn.tuk <- glht(n2o.day21.aov2, linfct = mcp(int = "Tukey"))
summary(dn.tuk)          # standard display
dn.tuk.cld <- cld(dn.tuk)   # letter-based display
par(mai=c(1,1,1.1,0.2),  # mai specifies margin size in inches
    mfrow=c(1,1),
    mgp=c(2,1,0),
    cex=0.8) 
plot(dn.tuk.cld, cex=0.5, las=2, col="grey", xaxt="n")
axis(1, at=gas15.cum$int,labels=FALSE)
text(gas15.cum$int, labels=gas15.cum$int, par("usr")[3], adj=c(1.2,1.2), xpd=TRUE, srt=45, cex=0.8)


n2o.day21 <- with(isom.dates.exp2, list( Mean=round(tapply(n2o.day21, list(treat,soil), mean),1),
                                    SE=round(tapply(n2o.day21, list(treat,soil), se),1)))
data.frame(n2o.day21)[c(1,3,2,4)];rm(n2o.day21)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### CO2 ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
co2.day3.aov <- aov(co2.day3 ~ soil*Lt*Fc, isom.dates.exp2); summary(co2.day3.aov)
co2.day3.aov <- update( co2.day3.aov,log(.)~ soil+Lt+Fc, isom.dates.exp2); summary(co2.day3.aov)
co2.day3.aov <- update( co2.day3.aov,log(.)~ soil*treat, isom.dates.exp2); summary(co2.day3.aov)

par(mfrow=c(2,2))
plot(co2.day3.aov)


co2.day3.aov2 <- update( co2.day3.aov,.~ int, isom.dates.exp2); summary(co2.day3.aov)
dn.tuk <- glht(co2.day3.aov2, linfct = mcp(int = "Tukey"))
summary(dn.tuk)          # standard display
dn.tuk.cld <- cld(dn.tuk)   # letter-based display
par(mai=c(1,1,1.1,0.2),  # mai specifies margin size in inches
    mfrow=c(1,1),
    mgp=c(2,1,0),
    cex=0.8) 
plot(dn.tuk.cld, cex=0.5, las=2, col="grey", xaxt="n")
axis(1, at=gas15.cum$int,labels=FALSE)
text(gas15.cum$int, labels=gas15.cum$int, par("usr")[3], adj=c(1.2,1.2), xpd=TRUE, srt=45, cex=0.8)




co2.day3 <- with(isom.dates.exp2, list( Mean=round(tapply(co2.day3, list(treat,soil), mean),1),
                                        SE=round(tapply(co2.day3, list(treat,soil), se),1)))
data.frame(co2.day3)[c(1,3,2,4)];rm(co2.day3)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
co2.day21.aov <- aov(co2.day21 ~ soil*Lt*Fc, isom.dates.exp2); summary(co2.day21.aov)
co2.day21.aov <- aov(co2.day21 ~ soil+Lt+Fc, isom.dates.exp2); summary(co2.day21.aov)
co2.day21.aov <- update( co2.day21.aov,log(.)~ soil*treat, isom.dates.exp2); summary(co2.day21.aov)

par(mfrow=c(2,2))
plot(co2.day21.aov)


co2.day21.aov2 <- update( co2.day21.aov,.~ int, isom.dates.exp2); summary(co2.day21.aov)
dn.tuk <- glht(co2.day21.aov2, linfct = mcp(int = "Tukey"))
summary(dn.tuk)          # standard display
dn.tuk.cld <- cld(dn.tuk)   # letter-based display
par(mai=c(1,1,1.1,0.2),  # mai specifies margin size in inches
    mfrow=c(1,1),
    mgp=c(2,1,0),
    cex=0.8) 
plot(dn.tuk.cld, cex=0.5, las=2, col="grey", xaxt="n")
axis(1, at=gas15.cum$int,labels=FALSE)
text(gas15.cum$int, labels=gas15.cum$int, par("usr")[3], adj=c(1.2,1.2), xpd=TRUE, srt=45, cex=0.8)



co2.day21 <- with(isom.dates.exp2, list( Mean=round(tapply(co2.day21, list(treat,soil), mean),1),
                                         SE=round(tapply(co2.day21, list(treat,soil), se),1)))
data.frame(co2.day21)[c(1,3,2,4)];rm(co2.day21)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





