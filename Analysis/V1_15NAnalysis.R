############################
# V1_Nitrogen
# 15N Analyses
# Quentin Schorpp
# 25.05.2015
###########################


### Analysis

source("Data/GatherSource/V1_MakeLikeFile.R")

# The data
# Nitox.pool and product ratio are AVERAGES, all other cumulated sums
# omit Loam without L. terrestris
gas15.omit <- gas15.cum[!(gas15.cum$soil=="Loam" & gas15.cum$Lt=="n"),]
binsize1 = sqrt(length(gas15.cum$product.ratio))


# Total [N2+N2O] ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Histogram
hist(gas15.cum$total.Ngas, breaks=binsize1)

# Boxplots
ggplot(gas15.cum, aes(x= Lt, y=total.Ngas),labeller=label_parsed) +
  stat_boxplot(geom="errorbar", coef=1.5, lwd=0.2) + 
  geom_boxplot(fill="light blue", lwd=0.2, outlier.size=0.4, outlier.shape=21) +
  ylab("total.Ngas") +
  xlab("italic(L.~terrestris)") +
  ggtitle("Total N2O") +
  facet_grid(~soil, labeller=label_parsed) + 
  scale_x_discrete(expression(italic(L.~terrestris)), labels=c("absent","present")) +
  mytheme

# Anova without Transformation
gas15.tng.aov <- aov(total.Ngas~soil*Lt*Fc, gas15.cum); summary(gas15.tng.aov)

# ANOVA with Transformation
gas15.tng.aov <- aov(log1p(total.Ngas)~soil*Lt*Fc, gas15.cum); summary(gas15.tng.aov)

# Model Simplification
gas15.tng.aov <- update(gas15.tng.aov,.~soil*Lt, gas15.cum); summary(gas15.tng.aov)

par(mfrow=c(2,2))
par(mar=c(3,3,3,3))
plot(gas15.tng.aov)

fligner.test(log1p(total.Ngas)~interaction(soil,Lt), gas15.cum) # Homoscedasticity

# Post Hoc Tukey Test
gas15.tng.aov2 <- update(gas15.tng.aov, .~ int.SLt)

tng.tuk <- glht(gas15.tng.aov2, linfct = mcp(int.SLt = "Tukey"))
summary(tng.tuk)          # standard display
tng.tuk.cld <- cld(tng.tuk)   # letter-based display
par(mai=c(1,1,1.1,0.2),  # mai specifies margin size in inches
    mfrow=c(1,1),
    mgp=c(2,1,0),
    cex=0.8) 
plot(tng.tuk.cld, cex=0.5, las=2, col="grey", xaxt="n")
axis(1, at=int.SLt,labels=FALSE)
text(int.SLt, labels=int.SLt, par("usr")[3], adj=c(1.2,1.2), xpd=TRUE, srt=45, cex=0.8)

with(gas15.cum, list( tapply(total.Ngas, list(soil,Lt), mean),
                      tapply(total.Ngas, list(soil,Lt), se)))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Total N2O ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Histogram
hist(gas15.cum$nitox.gc, breaks=binsize1)

# Boxplots
ggplot(gas15.cum, aes(x= Lt, y=nitox.gc),labeller=label_parsed) +
  stat_boxplot(geom="errorbar", coef=1.5, lwd=0.2) + 
  geom_boxplot(fill="light blue", lwd=0.2, outlier.size=0.4, outlier.shape=21) +
  ylab(txt.nitox.gc) +
  xlab("italic(L.~terrestris)") +
  ggtitle("Soil derived N2O") +
  facet_grid(~soil, labeller=label_parsed) + 
  scale_x_discrete(expression(italic(L.~terrestris)), labels=c("absent","present")) +
  mytheme

# Anova without Transformation
gas15.nxs.aov <- aov(nitox.gc~soil*Lt*Fc, gas15.cum); summary(gas15.nxs.aov)

# ANOVA with Transformation
gas15.nxs.aov <- aov(log1p(nitox.gc)~soil*Lt*Fc, gas15.cum); summary(gas15.nxs.aov)

# Model Simplification
gas15.nxs.aov <- update(gas15.nxs.aov,.~soil*Lt, gas15.cum); summary(gas15.nxs.aov)

par(mfrow=c(2,2))
par(mar=c(3,3,3,3))
plot(gas15.nxs.aov)

fligner.test(log1p(nitox.gc)~interaction(soil,Lt), gas15.cum) # Homoscedasticity

# Post Hoc Tukey Test
gas15.nxs.aov2 <- update(gas15.nxs.aov, .~ int.SLt)

nxs.tuk <- glht(gas15.nxs.aov2, linfct = mcp(int.SLt = "Tukey"))
summary(nxs.tuk)          # standard display
nxs.tuk.cld <- cld(nxs.tuk)   # letter-based display
par(mai=c(1,1,1.1,0.2),  # mai specifies margin size in inches
    mfrow=c(1,1),
    mgp=c(2,1,0),
    cex=0.8) 
plot(nxs.tuk.cld, cex=0.5, las=2, col="grey", xaxt="n")
axis(1, at=gas15.cum$int.SLt,labels=FALSE)
text(gas15.cum$int.SLt, labels=gas15.cum$int.SLt, par("usr")[3], adj=c(1.2,1.2), xpd=TRUE, srt=45, cex=0.8)

with(gas15.cum, list( tapply(nitox.gc, list(soil,Lt), mean),
                      tapply(nitox.gc, list(soil,Lt), se)))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Litter derived Total [N2+N2O] ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Histogram
hist(gas15.cum$total, breaks=binsize1)

# Boxplots
ggplot(gas15.cum, aes(x= Lt, y=total),labeller=label_parsed) +
  stat_boxplot(geom="errorbar", coef=1.5, lwd=0.2) + 
  geom_boxplot(fill="light blue", lwd=0.2, outlier.size=0.4, outlier.shape=21) +
  ylab(txt.total) +
  xlab("italic(L.~terrestris)") +
  ggtitle("Total N2 + N2O") +
  facet_grid(~soil, labeller=label_parsed) + 
  scale_x_discrete(expression(italic(L.~terrestris)), labels=c("absent","present")) +
  mytheme

# Anova without Transformation
gas15.tot.aov <- aov(total~soil*Lt*Fc, gas15.cum); summary(gas15.tot.aov)

# ANOVA with Transformation
gas15.tot.aov <- aov(log1p(total)~soil*Lt*Fc, gas15.cum); summary(gas15.tot.aov)

# Model Simplification
gas15.tot.aov <- update(gas15.tot.aov,.~soil*Lt, gas15.cum); summary(gas15.tot.aov)

par(mfrow=c(2,2))
par(mar=c(3,3,3,3))
plot(gas15.tot.aov)

fligner.test(log1p(total)~interaction(soil,Lt), gas15.cum) # Heteroscedasticity!

# Post Hoc Tukey Test
gas15.tot.aov2 <- update(gas15.tot.aov, .~ int.SLt)

tot.tuk <- glht(gas15.tot.aov2, linfct = mcp(int.SLt = "Tukey"))
summary(tot.tuk)          # standard display
tot.tuk.cld <- cld(tot.tuk)   # letter-based display
par(mai=c(1,1,1.1,0.2),  # mai specifies margin size in inches
    mfrow=c(1,1),
    mgp=c(2,1,0),
    cex=0.8) 
plot(tot.tuk.cld, cex=0.5, las=2, col="grey", xaxt="n")
axis(1, at=int.SLt,labels=FALSE)
text(int.SLt, labels=int.SLt, par("usr")[3], adj=c(1.2,1.2), xpd=TRUE, srt=45, cex=0.8)

with(gas15.cum, list( tapply(total, list(soil,Lt), mean),
                      tapply(total, list(soil,Lt), se)))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Litter derived N2 ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Histogram
hist(gas15.cum$dinitrogen, breaks=binsize1)

# Boxplots
ggplot(gas15.cum, aes(x= Lt, y=dinitrogen),labeller=label_parsed) +
  stat_boxplot(geom="errorbar", coef=1.5, lwd=0.2) + 
  geom_boxplot(fill="light blue", lwd=0.2, outlier.size=0.4, outlier.shape=21) +
  ylab(txt.dinitrogen) +
  xlab("italic(L.~terrestris)") +
  ggtitle("Litter derived N2") +
  facet_grid(~soil, labeller=label_parsed) + 
  scale_x_discrete(expression(italic(L.~terrestris)), labels=c("absent","present")) +
  mytheme
# Two outlier!

# Anova without Transformation
gas15.dn.aov <- aov(dinitrogen~soil*Lt*Fc, gas15.cum); summary(gas15.dn.aov)

# ANOVA with Transformation
gas15.dn.aov <- aov(log1p(dinitrogen)~soil*Lt*Fc, gas15.cum); summary(gas15.dn.aov)

# Model Simplification
gas15.dn.aov <- update(gas15.dn.aov,.~soil*Lt, gas15.cum); summary(gas15.dn.aov)

par(mfrow=c(2,2))
par(mar=c(3,3,3,3))
plot(gas15.dn.aov)

fligner.test(log1p(dinitrogen)~interaction(soil,Lt), gas15.cum) # Homoscedasticity

# Post Hoc Tukey Test
gas15.dn.aov2 <- update(gas15.dn.aov, .~ int.SLt)

dn.tuk <- glht(gas15.dn.aov2, linfct = mcp(int.SLt = "Tukey"))
summary(dn.tuk)          # standard display
dn.tuk.cld <- cld(dn.tuk)   # letter-based display
par(mai=c(1,1,1.1,0.2),  # mai specifies margin size in inches
    mfrow=c(1,1),
    mgp=c(2,1,0),
    cex=0.8) 
plot(dn.tuk.cld, cex=0.5, las=2, col="grey", xaxt="n")
axis(1, at=gas15.cum$int.SLt,labels=FALSE)
text(gas15.cum$int.SLt, labels=gas15.cum$int.SLt, par("usr")[3], adj=c(1.2,1.2), xpd=TRUE, srt=45, cex=0.8)

with(gas15.cum, list( tapply(dinitrogen, list(soil,Lt), mean),
                      tapply(dinitrogen, list(soil,Lt), se)))

#For Table 3
dinitrogen <- with(gas15.cum, list( Mean=round(tapply(dinitrogen, list(treat,soil), mean),1),
                               SE=round(tapply(dinitrogen, list(treat,soil), se),1)))
data.frame(dinitrogen)[c(1,3,2,4)];rm(dinitrogen)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Litter derived N2O ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Histogram
hist(gas15.cum$nitox, breaks=binsize1)

# Boxplots 
p1 <- ggplot(gas15.cum, aes(x= Lt, y=nitox),labeller=label_parsed) +
        stat_boxplot(geom="errorbar", coef=1.5, lwd=0.2) + 
        geom_boxplot(fill="light blue", lwd=0.2, outlier.size=0.4, outlier.shape=21) +
        ylab(txt.nitox) +
        xlab("italic(L.~terrestris)") +
        ggtitle(expression(paste("Litter derived ", N[2], "O"))) +
        facet_grid(~soil, labeller=label_parsed) + 
        scale_x_discrete(expression(italic(L.~terrestris)), labels=c("absent","present")) +
        mytheme
p1

#ANOVA without transformation
gas15.nx.aov <- aov(nitox~soil*Lt*Fc, gas15.cum); summary(gas15.nx.aov)

# ANOVA with log(x+1) Transformation
gas15.nx.aov <- aov(log1p(nitox)~soil*Lt*Fc, gas15.cum); summary(gas15.nx.aov)

# Model Simplification
gas15.nx.aov <- update(gas15.nx.aov,.~soil*Lt, gas15.cum); summary(gas15.nx.aov)

par(mfrow=c(2,2))
par(mar=c(3,3,3,3))
plot(gas15.nx.aov)

fligner.test(log1p(nitox)~interaction(soil,Lt), gas15.cum) # Heteroscedasticity!

# Post Hoc Tukey Test
gas15.nx.aov2 <- update(gas15.nx.aov, .~ int.SLt)

nx.tuk <- glht(gas15.nx.aov2, linfct = mcp(int.SLt = "Tukey"))
summary(nx.tuk)          # standard display
nx.tuk.cld <- cld(nx.tuk)   # letter-based display
par(mai=c(1,1,1.1,0.2),  # mai specifies margin size in inches
    mfrow=c(1,1),
    mgp=c(2,1,0),
    cex=0.8) 
plot(nx.tuk.cld, cex=0.5, las=2, col="grey", xaxt="n")
axis(1, at=int.SLt,labels=FALSE)
text(int.SLt, labels=int.SLt, par("usr")[3], adj=c(1.2,1.2), xpd=TRUE, srt=45, cex=0.8)

with(gas15.cum, list( tapply(nitox, list(Lt,soil), mean),
                      tapply(nitox, list(Lt,soil), se)))

#For Table 3
nitox <- with(gas15.cum, list( Mean=round(tapply(nitox, list(treat,soil), mean),1),
                               SE=round(tapply(nitox, list(treat,soil), se),1)))
data.frame(nitox)[c(1,3,2,4)];rm(nitox)

boxplot(expm1(fitted(gas15.nx.aov))~gas15.cum$int.SLt)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Product Ratio ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Histogram
hist(gas15.cum$product.ratio, breaks=binsize1) # right skewed
hist(asin(gas15.cum$product.ratio/100)*2/pi, breaks=binsize1) # normalized

# Boxplots 
ggplot(gas15.cum, aes(x= Lt, y=product.ratio),labeller=label_parsed) +
  stat_boxplot(geom="errorbar", coef=1.5, lwd=0.2) + 
  geom_boxplot(fill="light blue", lwd=0.2, outlier.size=0.4, outlier.shape=21) +
  ylab(txt.ratio) +
  xlab("italic(L.~terrestris)") +
  ggtitle(expression(paste("Product Ratio"))) +
  facet_grid(~soil, labeller=label_parsed) + 
  scale_x_discrete(expression(italic(L.~terrestris)), labels=c("absent","present")) +
  mytheme


# ANOVA without transformation
gas15.pr.aov <- aov(product.ratio~soil*Lt*Fc, gas15.cum); summary(gas15.pr.aov)

# ANOVA with Arcsine transformation
# gas15.cum$pr.trans <- asin(gas15.cum$product.ratio/100)*2/pi # Arcsine transformation by Ben Bolker
gas15.cum$pr.trans <- asin(sqrt(gas15.cum$product.ratio))*180/pi
gas15.omit$pr.trans <- asin(sqrt(gas15.omit$product.ratio))*180/pi
gas15.pr.aov <- aov(pr.trans~soil*Lt*Fc, gas15.cum); summary(gas15.pr.aov)

# ANOVA with model simplification
gas15.pr.aov <- update(gas15.pr.aov, .~soil*Lt);summary(gas15.pr.aov)

# Validation 
par(mfrow=c(2,2))
par(mar=c(3,3,3,3))
plot(gas15.pr.aov)

fligner.test(pr.trans ~ interaction(soil,Lt,Fc), gas15.cum) # not significant = Homoscedasticity
with(gas15.cum, shapiro.test(pr.trans)) # if this Test is significant, the data is significantly different from Normal sitribution and the ANOVA assumption is violated


# Post Hoc Tukey Test
gas15.pr.aov2 <- aov(pr.trans~int.SLt, data=gas15.cum)

pr.tuk <- glht(gas15.pr.aov2, linfct = mcp(int.SLt = "Tukey"))
summary(pr.tuk)          # standard display

pr.tuk.cld <- cld(pr.tuk)   # letter-based display
par(mai=c(1,1,1.1,0.2),  # mai specifies margin size in inches
    mfrow=c(1,1),
    mgp=c(2,1,0),
    cex=0.8) 
plot(pr.tuk.cld, cex=0.5, las=2, col="grey", xaxt="n")
axis(1, at=gas15.cum$int.SLt,labels=FALSE)
text(gas15.cum$int.SLt, labels=gas15.cum$int.SLt, par("usr")[3], adj=c(1.2,1.2), xpd=TRUE, srt=45, cex=0.8)


with(gas15.cum, list(round(tapply(product.ratio,interaction(Lt,soil), mean),2),
                     round(tapply(product.ratio,interaction(Lt,soil), se),2)))


# ANOVA wihtout L. terrestris x Loam
gas15.pr.aov3 <- aov(pr.trans~soil*Lt,gas15.omit) ;summary(gas15.pr.aov3)
par(mfrow=c(2,2))
plot(gas15.pr.aov3)
# Leaving out Loamy treatments without L. terrestris because their emissions have almost ambient air concentration, 
# shows significant effect of sandy soil.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# soil derived N2O ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Histogram
hist(gas15.cum$nitox.soil, breaks=binsize1)

# Boxplots
p2 <- ggplot(gas15.cum, aes(x= Lt, y=nitox.soil),labeller=label_parsed) +
        stat_boxplot(geom="errorbar", coef=1.5, lwd=0.2) + 
        geom_boxplot(fill="light blue", lwd=0.2, outlier.size=0.4, outlier.shape=21) +
        ylab(txt.nitox.soil) +
        xlab("italic(L.~terrestris)") +
        ggtitle(expression(paste("Soil derived ", N[2],"O"))) +
        facet_grid(~soil, labeller=label_parsed) + 
        scale_x_discrete(expression(italic(L.~terrestris)), labels=c("absent","present")) +
        mytheme
p2

# Anova without Transformation
gas15.nxt.aov <- aov(nitox.soil~soil*Lt*Fc, gas15.cum); summary(gas15.nxt.aov)

# ANOVA with Transformation
gas15.nxt.aov <- aov(log1p(nitox.soil)~soil*Lt*Fc, gas15.cum); summary(gas15.nxt.aov)

# Model Simplification
gas15.nxt.aov <- update(gas15.nxt.aov,.~soil*Lt, gas15.cum); summary(gas15.nxt.aov)

par(mfrow=c(2,2))
par(mar=c(3,3,3,3))
plot(gas15.nxt.aov)

fligner.test(log1p(nitox.soil)~interaction(soil,Lt), gas15.cum) # Homoscedasticity

# Post Hoc Tukey Test
gas15.nxt.aov2 <- update(gas15.nxt.aov, .~ int.SLt)

nxt.tuk <- glht(gas15.nxt.aov2, linfct = mcp(int.SLt = "Tukey"))
summary(nxt.tuk)          # standard display
nxt.tuk.cld <- cld(nxt.tuk)   # letter-based display
par(mai=c(1,1,1.1,0.2),  # mai specifies margin size in inches
    mfrow=c(1,1),
    mgp=c(2,1,0),
    cex=0.8) 
plot(nxt.tuk.cld, cex=0.5, las=2, col="grey", xaxt="n")
axis(1, at=gas15.cum$int.SLt,labels=FALSE)
text(gas15.cum$int.SLt, labels=gas15.cum$int.SLt, par("usr")[3], adj=c(1.2,1.2), xpd=TRUE, srt=45, cex=0.8)

with(gas15.cum, list( tapply(nitox.soil, list(treat,soil), mean),
                      tapply(nitox.soil, list(treat,soil), se)))

#For Table 3
nitox.soil <- with(gas15.cum, list( Mean=round(tapply(nitox.soil, list(treat,soil), mean),1),
                               SE=round(tapply(nitox.soil, list(treat,soil), se),1)))
data.frame(nitox.soil)[c(1,3,2,4)];rm(nitox.soil)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Fraction of litter derived N2O ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Histogram
hist(gas15.cum$nitox.pool, breaks=binsize1) # right skewed
hist(asin(gas15.cum$nitox.pool/100)*2/pi, breaks=binsize1) # normalized

# Boxplots 
p3 <- ggplot(gas15.cum, aes(x= Lt, y=nitox.pool),labeller=label_parsed) +
        stat_boxplot(geom="errorbar", coef=1.5, lwd=0.2) + 
        geom_boxplot(fill="light blue", lwd=0.2, outlier.size=0.4, outlier.shape=21) +
        ylab(txt.pool) +
        xlab("italic(L.~terrestris)") +
        ggtitle(expression(paste("Fraction of ",N[2],O[Litter]))) +
        facet_grid(~soil, labeller=label_parsed) + 
        scale_x_discrete(expression(italic(L.~terrestris)), labels=c("absent","present")) +
        mytheme
p3


# ANOVA without Transformation
gas15.np.aov <- aov(nitox.pool~soil*Lt*Fc, gas15.cum); summary(gas15.np.aov)

# ANOVA with Arcsine Transformation
gas15.cum$np.trans <- asin(sqrt(gas15.cum$nitox.pool))*180/pi
gas15.np.aov <- aov(np.trans~soil*Lt*Fc, gas15.cum); summary(gas15.np.aov)

# Model Simplification
gas15.np.aov <- update(gas15.np.aov,.~soil*Lt, gas15.cum); summary(gas15.np.aov)

par(mfrow=c(2,2))
par(mar=c(3,3,3,3))
plot(gas15.np.aov)

fligner.test(np.trans~interaction(soil,Lt), gas15.cum) # Heteroscedasticity!

# Post Hoc Tukey Test
gas15.np.aov2 <- aov(np.trans~int.SLt, data=gas15.cum)

np.tuk <- glht(gas15.np.aov2, linfct = mcp(int.SLt = "Tukey"))
summary(np.tuk)          # standard display

np.tuk.cld <- cld(np.tuk)   # letter-based display
par(mai=c(1,1,1.1,0.2),  # mai specifies margin size in inches
    mfrow=c(1,1),
    mgp=c(2,1,0),
    cex=0.8) 
plot(np.tuk.cld, cex=0.5, las=2, col="grey", xaxt="n")
axis(1, at=int.SLt,labels=FALSE)
text(int.SLt, labels=int.SLt, par("usr")[3], adj=c(1.2,1.2), xpd=TRUE, srt=45, cex=0.8)

with(gas15.cum, list( tapply(nitox.pool, list(soil,Lt), mean),
                      tapply(nitox.pool, list(soil,Lt), se)))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Ratio: litter derived N2O/soil derived N2O ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gas15.cum$lsratio <- gas15.cum$nitox/gas15.cum$nitox.gc

# Histogram
hist(gas15.cum$lsratio, breaks=binsize1) # right skewed
hist(asin(gas15.cum$lsratio/100)*2/pi, breaks=binsize1) # not normalized

# Boxplots 
ggplot(gas15.cum, aes(x= Lt, y=lsratio),labeller=label_parsed) +
  stat_boxplot(geom="errorbar", coef=1.5, lwd=0.2) + 
  geom_boxplot(fill="light blue", lwd=0.2, outlier.size=0.4, outlier.shape=21) +
  ylab(txt.pool) +
  xlab("italic(L.~terrestris)") +
  ggtitle(expression(paste("Fraction of ",N[2],O[Litter]))) +
  facet_grid(~soil, labeller=label_parsed) + 
  scale_x_discrete(expression(italic(L.~terrestris)), labels=c("absent","present")) +
  mytheme



# ANOVA without Transformation
gas15.lsr.aov <- aov(lsratio~soil*Lt*Fc, gas15.cum); summary(gas15.lsr.aov)

# ANOVA with Arcsine Transformation
gas15.cum$lsr.trans <- asin(sqrt(gas15.cum$lsratio))*180/pi
gas15.lsr.aov <- aov(lsratio~soil*Lt*Fc, gas15.cum); summary(gas15.lsr.aov)

# Model Simplification
gas15.lsr.aov <- update(gas15.lsr.aov,.~soil*Lt, gas15.cum); summary(gas15.lsr.aov)

par(mfrow=c(2,2))
par(mar=c(3,3,3,3))
plot(gas15.lsr.aov)

fligner.test(lsratio~interaction(soil,Lt), gas15.cum) # Heteroscedasticity!

# Post Hoc Tukey Test
gas15.lsr.aov2 <- aov(lsratio~int.SLt, data=gas15.cum)

lsr.tuk <- glht(gas15.lsr.aov2, linfct = mcp(int.SLt = "Tukey"))
summary(lsr.tuk)          # standard display

lsr.tuk.cld <- cld(lsr.tuk)   # letter-based display
par(mai=c(1,1,1.1,0.2),  # mai specifies margin size in inches
    mfrow=c(1,1),
    mgp=c(2,1,0),
    cex=0.8) 
plot(lsr.tuk.cld, cex=0.5, las=2, col="grey", xaxt="n")
axis(1, at=int.SLt,labels=FALSE)
text(int.SLt, labels=int.SLt, par("usr")[3], adj=c(1.2,1.2), xpd=TRUE, srt=45, cex=0.8)

with(gas15.cum, list( tapply(lsratio, list(soil,Lt), mean),
                      tapply(lsratio, list(soil,Lt), se)))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#####
t.test(total ~ Lt , gas15.cum[gas15.cum$soil=="Sand",]) # n.s.
t.test(nitox ~ Lt , gas15.cum[gas15.cum$soil=="Sand",]) # s.
t.test(nitox.gc ~ Lt , gas15.cum[gas15.cum$soil=="Sand",]) # s.
t.test(dinitrogen ~ Lt , gas15.cum[gas15.cum$soil=="Sand",]) # n.s.


#####
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 1)))
print(p1, vp = vplayout(1, 1))
print(p2, vp = vplayout(2, 1))
print(p3, vp = vplayout(3, 1))



