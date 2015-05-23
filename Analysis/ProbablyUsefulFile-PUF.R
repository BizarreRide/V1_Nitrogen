############################ Analysis of earthworm fresh weight ################################
#*************************************************************************************#
### Earthworm weightchange 
### paired one sided t.test (Rhiziya et.al2007)

setwd("D:/Quentin Schorpp/Arbeitsprozess/Mikrokosmosversuch/Organism growth/Worms_Biomass/ANOVA")
ewfw <- read.delim("D:/Quentin Schorpp/Schreibtisch/2 - Laborversuche/Mikrokosmen/Daten/Tierbesatz/Earthworm Biomass/EW_FrWe.txt")

library(ggplot2)
with(ewfw, list(length(before), range(before),
                range(after)))
ob=origin_bef = 9.25
range_bef= 9.25 - 10.25 = 1
oa=origin_aft = 6.75 
range_aft = 6.75 - 10 = 3.25

bs1 = binsize1 = diff(range(ewfw$before))/26
bs2 = binsize2 = diff(range(ewfw$after))/26

# Histograms,...
ewfw_melt = melt(ewfw[-c(7,8)], id.vars=c(1:4))
subggplot(ewfw_melt, aes(x=value)) + 
  geom_histogram( binwidth = 0.25, fill="white", colour="black") + 
  geom_line(stat="density") +
  xlim(6.5,11) + 
  facet_grid(variable ~ .)
# ... show that we can assume two real means, drawn from normal population

with(ewfw, t.test(before, after, paired=T))
with(ewfw, t.test(diff))
#?? One sided? t.test(before,after,paired=TRUE, mu = mean(before), alternative= "less")

# the varaince of a difference is the average of the summed varainces of each sample pair minus twice the covariance of these partners
# covarainceof two samples is the multiplication of standard deviations of two samples

############################ Analysis of collembolan population ################################
#*************************************************************************************#
setwd("D:/Quentin Schorpp/Arbeitsprozess/Mikrokosmosversuch/Organism growth/Colli_Population")
col.pop <- read.delim("D:/Quentin Schorpp/Schreibtisch/2 - Laborversuche/Mikrokosmen/Daten/Tierbesatz/Collembolan_Population/col_population.txt")

ggplot(col.pop, aes(x=surface)) + 
  geom_histogram(fill="white", colour="black") + 
  geom_density() + 
  facet_grid(soil+treat ~ .)

an2.0 <- aov(surface ~ exp*soil*treat, data=col.pop);  an2.0;  summary(an2.0)
par(mfrow=c(2,2))
plot(an2)

an2.1 <- aov(surface ~ soil*treat, subset(col.pop, exp %in% c("Exp1")));  an2.1;  summary(an2.1)
plot(an2.1)
an2.2 <- aov(surface ~ soil*treat, subset(col.pop, exp %in% c("Exp2")));  an2.2;  summary(an2.2)
plot(an2.2)

############################ Analysis of 15N Enrichment ################################
#*************************************************************************************#

#******************************** L. terrestris ******************************************

setwd("D:/Quentin Schorpp/Arbeitsprozess/Mikrokosmosversuch/15 Enrichment in animal tissues/15N Worms")
iso.worm <- read.delim("D:/Quentin Schorpp/Arbeitsprozess/Mikrokosmosversuch/15 Enrichment in animal tissues/15N Worms/iso_worms.txt")
attach(iso.worm)
MK = as.factor(MK)
op=par(mfrow=c(1,1), mar=c(2,2,2,2), oma=c(0,0,0,0))
hist(atpercent)
hist(atpercent, breaks=12) # Too narrow, breaks appear in Histogram
# what differnce in atpercent is meaningful?
max(atpercent) - min(atpercent) # Spannweite 
#There's a differnce between the Minima and Maxima of 0.86 atom percent

median(atpercent); mean(atpercent) # median and mean aren't exactly the same: mean=0.9 median=0.95 atom percent
icept = with(iso.worm[soil=="L",], min(atpercent))-(with(iso.worm[soil=="S",], min(atpercent))- with(iso.worm[soil=="L",], max(atpercent)))/2
with(iso.worm[soil=="L",], which.min(atpercent))
with(iso.worm[soil=="S",], which.max(atpercent))
maxRows <- by(iso.worm, iso.worm$soil, function(X) X[which.min(X$atpercent),])
do.call("rbind", maxRows)
MK[5]
MK[11]
which.max(atpercent,soil=="S")


iso.worm.plot <- qplot(soil, atpercent, colour=treat) + geom_point(size=2.5) +
  annotate("text", x=2,y=1.01+.02, label="max=1.01", vjust=0, col="red", size=4) +
  annotate("text", x=1,y=0.74+.02, label="min=0.74", vjust=0,col="red", size=4) +
  theme(panel.background= element_rect(fill="white"),
        panel.grid.major = element_line(colour="grey"),
        panel.grid.minor = element_line(colour="grey", linetype="dashed", size=0.2),
        panel.border= element_rect(colour="black", fill=NA, size=1),
        axis.text = element_text(colour="black", size=rel(1.3)),
        axis.text.x = element_text(face="italic"),
        axis.title=element_text(size=rel(1.4)),
        legend.text=element_text(face="italic"))
iso.worm.plot + xlab("Soiltype") + ylab(expression(paste(.^15,"N atom%",sep=""))) + ggtitle(expression(paste("Label in ",italic("L. terrestris"))))+
  scale_x_discrete(breaks=c("L","S"), labels=c("Loam", "Sand"))+
  labs(colour="Treatment")+
  scale_colour_discrete(labels=c("L.terrestris", "L. terrestris + F. candida"))


mod.fc1 <- aov(delta~treat*soil)
summary(mod.fc1)
mod.fc2 <- aov(delta~soil)
summary(mod.fc2)           
anova(mod.fc1,mod.fc2) # zwei df weniger und keine Verschlechterung des Modells für mod.fc2
AIC(mod.fc1,mod.fc2)

regfc1 <- lm(delta~treat*soil)
lm.influence(regfc1)
# point 9 is very influential on the intercept!
regfc2 <- lm(delta~soil)
lm.influence(regfc2)
influence.measures(regfc1)$is.inf
influence.measures(regfc2)$is.inf

t.test(atpercent~soil, iso.worm)

#******************************** F. candida ******************************************

iso.col <- read.delim("D:/Quentin Schorpp/Arbeitsprozess/Mikrokosmosversuch/15 Enrichment in animal tissues/15N Collembolans/iso_colli.txt")
attach(iso.col)
hist(atpercent)
hist(atpercent, breaks=12) # Too narrow, breaks appear in Histogram
# what differnce in atpercent is meaningful?
max(atpercent) - min(atpercent)
#There's a differnce between the Minima and Maxima of 0.896 atom percent

median(atpercent); mean(atpercent) # median and mean are the same: 6.4 atom percent
icept = with(iso.col[soil=="S",], min(atpercent))-(with(iso.col[soil=="S",], min(atpercent))- with(iso.col[soil=="L",], max(atpercent)))/2


iso.col.plot <- qplot(soil, atpercent, colour=treat) + geom_point(size=2.5) + geom_hline(yintercept=icept, linetype="dashed", colour="red") +
  annotate("text", label="Midline Max(Loam) - Min(Sand) y=6.4", x=1, y=6.44, size=4,colour="red") +
  theme(panel.background= element_rect(fill="white"),
        panel.grid.major = element_line(colour="grey"),
        panel.grid.minor = element_line(colour="grey", linetype="dashed", size=0.2),
        panel.border= element_rect(colour="black", fill=NA, size=1),
        axis.text = element_text(colour="black", size=rel(1.3)),
        axis.text.x = element_text(face="italic"),
        axis.title=element_text(size=rel(1.4)),
        legend.text=element_text(face="italic"))
iso.col.plot + xlab("Soiltype") + ylab(expression(paste(.^15,"N atom%",sep=""))) + ggtitle(expression(paste("Label in ",italic("F. candida"))))+
  scale_x_discrete(breaks=c("L","S"), labels=c("Loam", "Sand"))+
  labs(colour="Treatment")+
  scale_colour_discrete(labels=c("F.candida", "F. candida + L. terrestris"))


mod.fc1 <- aov(delta~treat*soil)
summary(mod.fc1)
mod.fc2 <- aov(delta~soil)
summary(mod.fc2)           
anova(mod.fc1,mod.fc2)
AIC(mod.fc1,mod.fc2)

regfc1 <- lm(delta~treat*soil)
lm.influence(regfc1)
# point 9 is very influential on the intercept!
regfc2 <- lm(delta~soil)
lm.influence(regfc2)
influence.measures(regfc1)$is.inf
influence.measures(regfc2)$is.inf

detach(iso.col)

t.test(atpercent~soil, iso.col)



##########Boxplots for both Collembolan and earthworm 15N Enrichment - Figure 1 ################
library(ggplot2)
library(plyr)

candw <- read.delim("D:/Quentin Schorpp/Arbeitsprozess/Mikrokosmosversuch/15 Enrichment in animal tissues/candw.txt")
candw$organism <- revalue(candw$organism, c(L.terrestris="italic('L.')~italic('terrestris')", F.candida = "italic('F.')~italic('candida')"))
candw$organism <- as.character(candw$organism)


# Base graph multipanel Boxplot -------------------------------------------

par(mfcol=c(2,2))
par(mar=c(0,0,0,0), oma = c(4,4,0.5,0.5))
par(cex=0.6)
par(tcl = -0.25)
par(mgp = c(2,0.6,0))
with(candw[(candw$organism=="L.terrestris"),], 
     boxplot(atpercent~soiltype, range=0, xaxt="n"))
box(col="grey60")
with(candw[(candw$organism=="F.candida"),], boxplot(atpercent~soiltype, range=0))
box(col="grey60")

plot(1,1, type="n", axes=F)
mtext("L. terrestris", side = 2, line = -1.5,adj = 0.98, cex = 1,col="black")
box(col="grey60")

plot(1,1, type="n", axes=F)
mtext("F. candida", side = 4, las=2, line = -1.5, adj = 0.98, cex = 1,col="black")
box(col="grey60")



library(lattice)
bwplot(atpercent~soiltype|organism, groups=organism, data=candw, ylim=c(0:7))

par(mfrow=c(1,1))
fig1.1 <- boxplot(atpercent~soiltype*organism, data=candw)
fig1.1$stats
df2 <- cbind(as.data.frame(fig1.1$stats),c("min","lq","m","uq","max"))

library(extrafont)
loadfonts()
fonts()

ggplot(candw, aes(x=soiltype, y=atpercent), family="Arial") + 
  stat_boxplot(geom="errorbar", coef=1.5, lwd=0.2) + 
  geom_boxplot(fill="grey", lwd=0.2, outlier.size=0.4, outlier.shape=21) +  
  facet_grid(organism ~ . , scales="free", space="free",labeller=label_parsed) + 
  ylab("atom% excess") + 
  xlab("soil texture") + 
  # ggtitle("15N enrichment in animal tissues") +
  theme_bw() + 
  theme(strip.background = element_rect(color = "light grey", fill="black", size=0.1),
        strip.text.y = element_text(size=8,  colour="white", face="italic"),
        axis.text.x = element_text(size=7),
        axis.title.x = element_text(size=8,face="bold", family="Times New Roman"),
        axis.text.y = element_text(size=7),
        axis.title.y = element_text(size=8,face="bold", family="Times New Roman"),
        axis.line = element_line(size=0.25),
        axis.ticks = element_line(size=0.25),
        panel.margin = unit(0, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black", size=0.2, fill=NA))

setwd("D:/Quentin Schorpp/Schreibtisch/2 - Laborversuche/Mikrokosmen/Manuskript/Entwurf")
ggsave("Figure1.pdf", width=9, height=7, units="cm", useDingbats=FALSE) 
ggsave("Figure1.svg", width=9, height=7, units="cm") 

############################ Analyse der kumulativen Gas emissionen ################################
#*************************************************************************************#
##Cumulative N2O emissions

setwd("D:/Quentin Schorpp/Arbeitsprozess/Mikrokosmosversuch/Gas Fluxes/N2O/Cumulative N2O emission")
nitox.cum <- read.delim("D:/Quentin Schorpp/Arbeitsprozess/Mikrokosmosversuch/Gas Fluxes/N2O/Cumulative N2O emission/nitox_cum.txt")

std <- par(mfrow=c(1,1))

## Testen
library(lattice)

#******************* Main-experiment: labelled leaf litter*******************

## assumptions
with(subset(nitox.cum,label==T), fligner.test(log(Kum_N2O) ~ treat*soil))
with(subset(nitox.cum,label==T), bartlett.test(log(Kum_N2O) ~ treat*soil))

# model with presence of soil organisms as factor
mod1<- with(subset(nitox.cum,nitox.cum$label==T), aov(log(Kum_N2O) ~ Lt*Fc*soil))
mod2<- with(subset(nitox.cum,nitox.cum$label==T), aov(log(Kum_N2O) ~ treat*soil))


AIC(mod1,mod2)

mod1, mod2
summary(mod1); summary(mod2)
summary.lm(mod1); summary.lm(mod2)


TukeyHSD(mod1)
par(mfrow=c(1,1))
plot(TukeyHSD(mod1), las=1)

TukeyHSD(mod2)
par(mfrow=c(1,1))
plot(TukeyHSD(mod2), las=1)

with(subset(nitox.cum,nitox.cum$label==T),boxplot(Kum_N2O ~ Lt*soil))

#************ Control experiment, cumulative emissions, calculated for period after day 7

# check assumptions
with(subset(nitox.cum,label==F), fligner.test(log(Kum_N2O) ~ treat*soil))
with(subset(nitox.cum,label==F), bartlett.test(log(Kum_N2O) ~ treat*soil))
#assumptions not violated!!!

coplot(Kum_N2O ~ treat|soil, data=nitox.cum[31:54,])

mod3 <- with(subset(nitox.cum,nitox.cum$label==F), aov(log(Kum_N2O) ~ treat*soil))

op <- par(mfrow = c(2,2), mar=c(5,4,1,2))
plot(mod3, add.smooth=F, which=1)
E <- resid(mod3)
hist(E, xlab="residuals", main="")
with(subset(nitox.cum,nitox.cum$label==F), plot(treat, E, xlab="Treatment", ylab="Residuals"))
with(subset(nitox.cum,nitox.cum$label==F), plot(soil, E, xlab="soiltype", ylab="Residuals"))

# ??? graphical assumptions check ???


# Does the food type has a significant effect on cumulative N2O emmissions?
# check with data for all days

modmod1 <- with(subset(nitox.cum,nitox.cum$label==F), aov(log(Kum_N2O) ~ Lt*Fc+soil+food))
summary(modmod1)
drop1(modmod1, test="F")
# warum fehlt die Interaktion Lt:Fc in der summary?

# check with data for days 7+
modmod2 <- with(subset(nitox.cum,nitox.cum$label==F), aov(log(Kum_N2O_7a) ~ Lt*Fc+soil+food))
summary(modmod2)
drop1(modmod2, test="F")

# sehr wiedersprüchlich, warum wird das modell besser,
# wenn der einzige signifikante faktor (soil) gedropped wird??

modmod3 <- with(subset(nitox.cum,nitox.cum$label==F & nitox.cum$soil=="Sand"), aov(log(Kum_N2O_7a) ~ Lt*Fc+food))
summary(modmod3)
drop1(modmod3, test="F")

with(subset(nitox.cum,nitox.cum$label==F & nitox.cum$soil=="Sand"), boxplot(Kum_N2O_7a ~ food))
with(subset(nitox.cum,nitox.cum$label==F & nitox.cum$soil=="Sand"), boxplot(Kum_N2O_7a ~ treat))


# aus Trotz
# Spiegelung der Analyse von mainexperiment:

mod4 <- with(subset(nitox.cum,nitox.cum$label==F), aov(log(Kum_N2O) ~ Lt*Fc*soil))
mod5 <- with(subset(nitox.cum,nitox.cum$label==F), aov(log(Kum_N2O) ~ treat*soil))

AIC(mod4,mod5) # identisch

mod4
summary(mod4)
summary.lm(mod4)

TukeyHSD(mod4)
par(mfrow=c(1,1), mar=c(5,9,3,3))
plot(TukeyHSD(mod4), las=1)

par(mfrow=c(1,1), mar=c(2,2,2,2), oma=c(0,0,0,0))
with(subset(nitox.cum,nitox.cum$label==F),boxplot(log(Kum_N2O) ~ Lt*soil))
with(subset(nitox.cum,nitox.cum$label==F),boxplot(log(Kum_N2O_7a) ~ Fc))
with(subset(nitox.cum,nitox.cum$label==F),boxplot(Kum_N2O ~ Fc*Lt))
summary(subset(nitox.cum,nitox.cum$label==F))

# Ende Spiegelung     

# Diesselbe Analyse mit fixed Variance structure

library(nlme)

mod4.lm <- with(subset(nitox.cum,nitox.cum$label==F), gls(Kum_N2O_7a ~ treat*soil))
vf1 <- varIdent(form ~ 1|soil)
mod4.lm2 <- with(subset(nitox.cum,nitox.cum$label==F), gls(Kum_N2O_7a ~ treat*soil, weights = vf1))
anova(mod4.lm, mod4.lm2)

# Kein Unterschied!! hat funktioniert??



##Cumulative CO2 emissions

setwd("D:/Quentin Schorpp/Arbeitsprozess/Mikrokosmosversuch/Gas Fluxes/CO2")
cardi.cum <- read.delim("D:/Quentin Schorpp/Arbeitsprozess/Mikrokosmosversuch/Gas Fluxes/CO2/cardi_cum.txt")
std <- par(mfrow=c(1,1))

## Testen
library(lattice)

# Main-experiment: labelled leaf litter

## assumptions
with(subset(cardi.cum,label==T), fligner.test(log(Kum_C) ~ treat*soil))
with(subset(cardi.cum,label==T), bartlett.test(log(Kum_C) ~ treat*soil))

# model with presence of soil organisms as factor
mod1<- with(subset(cardi.cum,cardi.cum$label==T), aov(log(Kum_C) ~ Lt*Fc*soil))
mod2<- with(subset(cardi.cum,cardi.cum$label==T), aov(log(Kum_C) ~ treat*soil))


AIC(mod1,mod2)

mod1; mod2
summary(mod1); summary(mod2)
summary.lm(mod1); summary.lm(mod2)


TukeyHSD(mod1)
par(mfrow=c(1,1))
plot(TukeyHSD(mod1), las=1)

TukeyHSD(mod2)
par(mfrow=c(1,1))
plot(TukeyHSD(mod2), las=1)

with(subset(cardi.cum,cardi.cum$label==T),boxplot(Kum_C ~ Lt*soil))

# Control experiment, cumulative emissions, calculated for period after day 7

# check assumptions
with(subset(cardi.cum,label==F), fligner.test(log(Kum_C) ~ treat*soil))
with(subset(cardi.cum,label==F), bartlett.test(log(Kum_C) ~ treat*soil))
#assumptions violated!!!

coplot(Kum_C_7b ~ treat|soil, data=cardi.cum[31:54,])

mod3 <- with(subset(cardi.cum,cardi.cum$label==F), aov(log(Kum_C_7b) ~ treat*soil))

op <- par(mfrow = c(2,2), mar=c(5,4,1,2))
plot(mod3, add.smooth=F, which=1)
E <- resid(mod3)
hist(E, xlab="residuals", main="")
with(subset(cardi.cum,cardi.cum$label==F), plot(treat, E, xlab="Treatment", ylab="Residuals"))
with(subset(cardi.cum,cardi.cum$label==F), plot(soil, E, xlab="soiltype", ylab="Residuals"))

# ??? graphical assumptions check ???

# go on despite violateed assumptions

# Does the food type has a significant effect on cumulative N2O emmissions?
# check with data for all days

modmod1 <- with(subset(cardi.cum,cardi.cum$label==F), aov(log(Kum_C) ~ Lt*Fc+soil+food))
summary(modmod1)
drop1(modmod1, test="F")
# warum fehlt die Interaktion Lt:Fc in der summary?

# check with data for days 7+
modmod2 <- with(subset(cardi.cum,cardi.cum$label==F), aov(log(Kum_C) ~ Lt*Fc+soil+food))
summary(modmod2)
drop1(modmod2, test="F")

# sehr wiedersprüchlich, warum wird das modell besser,
# wenn der einzige signifikante faktor (soil) gedropped wird??

modmod3 <- with(subset(cardi.cum,cardi.cum$label==F & cardi.cum$soil=="Sand"), aov(log(Kum_C) ~ Lt*Fc+food))
summary(modmod3)
drop1(modmod3, test="F")

with(subset(cardi.cum,cardi.cum$label==F & cardi.cum$soil=="Sand"), boxplot(Kum_C ~ food))
with(subset(cardi.cum,cardi.cum$label==F & cardi.cum$soil=="Sand"), boxplot(Kum_C ~ treat))


# aus Trotz
# Spiegelung der Analyse von mainexperiment:

mod4 <- with(subset(cardi.cum,cardi.cum$label==F), aov(log(Kum_C) ~ Lt*Fc*soil))
mod5 <- with(subset(cardi.cum,cardi.cum$label==F), aov(log(Kum_C) ~ treat*soil))

AIC(mod4,mod5) # identisch

mod4
summary(mod4)
summary.lm(mod4)

TukeyHSD(mod4)
par(mfrow=c(1,1), mar=c(5,9,3,3))
plot(TukeyHSD(mod4), las=1)

par(mfrow=c(1,1))
with(subset(cardi.cum,cardi.cum$label==F),boxplot(log(Kum_C) ~ Lt*soil))
with(subset(cardi.cum,cardi.cum$label==F),boxplot(log(Kum_C) ~ Fc))

# Ende Spiegelung     

# Diesselbe Analyse mit fixed Variance structure

library(nlme)

mod4.lm <- with(subset(cardi.cum,cardi.cum$label==F), gls(Kum_C ~ treat*soil))
vf1 <- varIdent(form ~ 1|soil)
mod4.lm2 <- with(subset(cardi.cum,cardi.cum$label==F), gls(Kum_C ~ treat*soil, weights = vf1))
anova(mod4.lm, mod4.lm2)

# Kein Unterschied!! hat funktioniert??

# Figure2, N2O, N2O-Cum, CO2-Cum, Multipanel Plot -------------------------

library(reshape2)
library(ggplot2)
library(grid)
library(plyr)

setwd("D:/Quentin Schorpp/Arbeitsprozess/Mikrokosmosversuch/Gas Fluxes/N2O")
fig2 <- read.delim("D:/Quentin Schorpp/Arbeitsprozess/Mikrokosmosversuch/Gas Fluxes/N2O/Fig_2_data.txt")
fig2  <- melt(fig2, id.vars= c(1:5)) 
fig2  <- dcast(fig2, treat+soil+variable+day+Experiment~measure) 
fig2$variable <- revalue(fig2$variable, c(nitox="paste(N[2],O)~Production~Rate", nitox_cum="Cumulative~paste(N[2],O)", carbon="Cumulative~CO[2]"))
fig2$variable <- factor(fig2$variable, levels=c("paste(N[2],O)~Production~Rate","Cumulative~paste(N[2],O)", "Cumulative~CO[2]"))

txt.cum_nitox = expression(paste(CO[2],"[µg ",CO[2],"-C ",kg^-1,"]","    ",N[2],"O [µg ",N[2],"O-N ",kg^-1,"]","   ", N[2],"O [µg ",N[2],"O-N ",kg^-1, h^-1,"]"))
#library(extrafont)  
#loadfonts()
fonts()
pd = position_dodge(0.2)
ggplot(data=fig2, aes(x=day, y=mean, group=treat, shape=treat)) + 
  geom_point(size=1.5, position=pd, fill="white") + 
  geom_line (size=0.35, position=pd) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=4, size=0.2, position=pd) + 
  xlab("Days [d]") +
  ylab(txt.cum_nitox) +
  facet_grid(variable~soil, scales ="free", labeller=label_parsed) + 
  scale_shape_manual(values=c(21,16,15,17), 
                     labels=c("Control","F. candida","Interaction", "L. terrestris")) + 
  #                   limits=c("L. terrestris","F. candida","Interaction","Control")) +
  labs(shape="Treatment:") +
  theme_bw()+
  theme(
    strip.background = element_rect(fill="black"),
    strip.text.x = element_text(size=8, face="bold", colour="white"),
    strip.text.y = element_text(size=8, face="bold", colour="white"),
    panel.margin = unit(0, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour="black", size=0.2, fill=NA),
    axis.text.x = element_text(size=7),
    axis.title.x = element_text(size=8,face="bold", family="Times New Roman"),
    axis.text.y = element_text(size=7),
    axis.title.y = element_text(size=8,face="bold", family="Times New Roman"),
    axis.ticks = element_line(size=0.25),
    legend.key=element_blank(),
    legend.background=element_blank(),
    legend.text=element_text(size=8,face="italic", family="Times New Roman"),
    legend.title=element_text(size=8),
    legend.position=c(0.55,0.55) )
#legend.position="bottom")





ggsave("Figure 2c.pdf", width=19, height=15, units="cm", useDingbats=FALSE)
ggsave("Figure 2c.svg", width=19, height=15, units="cm")


#color palette for treatment Graphics
trt_pal = c("firebrick3", "dodgerblue2", "darkorchid", "dimgray")

#universal theme
mytheme <- theme_bw() + 
  theme(axis.title = element_text(size = rel(1.25)),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.background = element_rect(fill="black"),
        strip.text.x = element_text(size=11, face="bold", colour="white"),
        strip.text.y = element_text(size=11, face="bold", colour="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),)

# y-Axis labels
txt.total = expression(paste("Total ", N[2],"O+",N[2]," [µg ",N[2],"(O)-N ",kg^-1, h^-1,"]"))

txt.cum_nitox = expression(paste("Dinitrogen ", N[2]," [µg ",N[2],"-N ",kg^-1,"]"))
txt.ratio = expression(paste("Product Ratio ", N[2],"O/",N[2],"+",N[2],"O"))
txt.pool = expression(paste("Plant derived ", N[2],"O","[%]"))

txt.nitox = N[2]O~Production~Rate
txt.nitox = expression(paste(N[2],"O [µg ",N[2],"O-N ",kg^-1, h^-1,"]"))


mytheme + 
  scale_shape_manual(values=c(15,16, 18,17)) +
  #scale_colour_manual(values=trt_pal) +
  scale_colour_grey(start=0, end=0.6) +
  theme(legend.position="none") + 
  theme(panel.margin = unit(0, "lines"))
pfig2


############################ Analyse der 15N Gas-Daten ################################
#*************************************************************************************#

######### Daten einlesen ###########

setwd("D:/Quentin Schorpp/Arbeitsprozess/Mikrokosmosversuch/15N_gas/Analyse_02.04.14")
gas15 <- read.delim("D:/Quentin Schorpp/Schreibtisch/2 - Laborversuche/Mikrokosmen/Daten/Isotope 15N-13C/15N Gas/V1_15Ngas.txt")
str(gas15)

gas15$date <- as.POSIXlt(gas15$date)
gas15$treat <- factor(gas15$treat, levels = c("Lt", "Fc", "LF", "C"))
str(gas15)

gas_time = cbind(day=gas15$day, soil=gas15$soil, treat= gas15$treat, gas15[,20:24])

# subsets
ident = gas15[,1:3]
treatments = gas15[,8:13]
time = gas15[,4:7]
gas = gas15[,19:24]

########## Data summary ##########
gas15.cum2$troil = (with(gas15.cum2, interaction(treat,soil)))
gas15.cum2$fnloss = (with(gas15.cum2, dinitrogen/nitox))
gas15.cum$troil = (with(gas15.cum, interaction(treat,soil)))
gas15.cum$fnloss = (with(gas15.cum, dinitrogen/nitox))
gas15.cum$fnloss[10] = 254
gas15.cum$fnloss[13] = 84
gas15.cum$product.ratio = (with(gas15.cum, nitox/(dinitrogen+nitox)))

with(gas15.cum,plot(fnloss~troil))
with(gas15.cum2,plot(fnloss~troil*week))

with(gas15.cum, tapply(fnloss,troil,mean))
with(gas15.cum, tapply(product.ratio,troil,mean))
with(gas15.cum, boxplot(product.ratio~Lt*Fc*soil))

with(gas15.cum, mean(fnloss))



########## Data Exploration + Graphics ##########
library(ggplot2)
library(reshape2)

#color palette for treatment Graphics
trt_pal = c("firebrick3", "dodgerblue2", "darkorchid", "dimgray")

# y-Axis labels
txt.total = expression(paste("Total ", N[2],"O+",N[2]," [µg ",N[2],"(O)-N ",kg^-1, h^-1,"]"))

txt.dinitrogen = expression(paste("Dinitrogen ", N[2]," [µg ",N[2],"-N ",kg^-1, h^-1,"]"))
txt.ratio = expression(paste("Product Ratio ", N[2],"O/",N[2],"+",N[2],"O"))
txt.pool = expression(paste("Plant derived ", N[2],"O","[%]"))
txt.nitox = expression(paste("Nitrous Oxide ", N[2],"O [µg ",N[2],"O-N ",kg^-1, h^-1,"]"))
# universal theme
mytheme <- theme_bw() + 
  theme(axis.title = element_text(size = rel(1.25)),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        strip.background = element_rect(fill="black"),
        strip.text.x = element_text(size=14, face="bold", colour="white"),
        strip.text.y = element_text(size=14, face="bold", colour="white"))

# create data frmae for plotting
gas_time = cbind(day=gas15$day, soil=gas15$soil, treat= gas15$treat, gas15[,20:24])

# Single Panel plot
p1 = ggplot(data=gas15, aes(x=day, y=total, colour=treat, shape=treat)) +
  facet_grid(~soil) +
  #guides(colour=F, shape=F) +
  geom_point(size=2.5, position=position_dodge(width=1)) +
  scale_shape_manual(values=c(15,16, 18,17)) +
  #scale_colour_manual(values=trt_pal) +
  scale_colour_grey(start=0, end=0.6) +
  scale_x_continuous(name="Days [d]", breaks=gas_time$day) +
  scale_y_continuous(breaks=seq(0,14,2)) +
  ylab(txt.total)
p1 +  mytheme


##### Multi-Panel Dotplot #####
#_______________________________________________________________________________________

gas_time2 <- melt(gas15, id.vars=c("day","treat","soil"), measure.vars=c("total","product.ratio","nitox"))
gas_time2$variable <- factor(gas_time2$variable, levels = c("total", "product.ratio", "nitox"))
# For better names in facet.grid:
levels(gas_time2$variable)[levels(gas_time2$variable)=="total"] <- "N2O+N2"
levels(gas_time2$variable)[levels(gas_time2$variable)=="product.ratio"] <- "Product~Ratio"
levels(gas_time2$variable)[levels(gas_time2$variable)=="nitox"] <- "N2O"

p2 = ggplot(data=gas_time2, aes(x=day, y=value, colour=treat, shape=treat)) +
  facet_grid(variable~soil, scales="free", labeller=label_parsed) +
  #guides(colour=F, shape=F) +
  geom_point(size=2.5, position=position_dodge(width=1)) +
  scale_shape_manual(values=c(15,16, 18,17)) +
  #scale_colour_manual(values=trt_pal) +
  scale_colour_grey(start=0, end=0.6) +
  scale_x_continuous(name="Days [d]", breaks=gas_time$day) +
  ylab("")
p2 + mytheme + theme(legend.position="none")


##### Multi-Panel Lineplot with mean and standard error #####

#**** create dataset with mean +standard error ****
#________________________________________________________________________________________________________________________________

head1 = c("treat","day","soil", "total","product.ratio","nitox","dinitrogen","nitox.pool")

# df for mean over time(day)
gas_mean2     = matrix(0,32,8)
gas_mean2[,1] =rep(c("Lt", "Fc", "LF", "C"),c(4,4,4,4))
gas_mean2[,2] =rep(c(19, 26, 33, 40),4)
gas_mean2[,3] =rep(c("Loam", "Sand"),c(16,16))
gas_mean2[,4] = with(gas15, tapply(total, list(day,treat,soil), mean))
gas_mean2[,5] = with(gas15, tapply(product.ratio, list(day,treat,soil), mean))
gas_mean2[,6] = with(gas15, tapply(nitox, list(day,treat,soil), mean))
gas_mean2[,7] = with(gas15, tapply(dinitrogen, list(day,treat,soil), mean))
gas_mean2[,8] = with(gas15, tapply(nitox.pool, list(day,treat,soil), mean))
colnames(gas_mean2) = head1
gas_mean2 <- as.data.frame(gas_mean2)

# df for standard error over time(day)
gas_se2     = matrix(0,32,8)
gas_se2[,1] =rep(c("Lt", "Fc", "LF", "C"),c(4,4,4,4))
gas_se2[,2] =rep(c(19, 26, 33, 40),4)
gas_se2[,3] =rep(c("Loam", "Sand"),c(16,16))
gas_se2[,4] = with(gas15, tapply(total, list(day,treat,soil), function(x) sqrt(var(x)/length(x))))
gas_se2[,5] = with(gas15, tapply(product.ratio, list(day,treat,soil), function(x) sqrt(var(x)/length(x))))
gas_se2[,6] = with(gas15, tapply(nitox, list(day,treat,soil), function(x) sqrt(var(x)/length(x))))
gas_se2[,7] = with(gas15, tapply(dinitrogen, list(day,treat,soil), function(x) sqrt(var(x)/length(x))))
gas_se2[,8] = with(gas15, tapply(nitox.pool, list(day,treat,soil), function(x) sqrt(var(x)/length(x))))
colnames(gas_se2) = head1
gas_se2 <- as.data.frame(gas_se2)

#**** prepare dataframe for multi panel plot ****
#________________________________________________________________________________________________________________________________
gas_mean2_melted            <- melt(gas_mean2, id.vars= c(1:3)) 
gas_se2_melted              <- melt(gas_se2, id.vars = c(1:3)) 

## combine mean and standard error in one dataframe
gas_meanse_melted           <- cbind(gas_mean2_melted[1:96,], gas_se2_melted[1:96,5])
colnames(gas_meanse_melted)[6] <- "se"

# Reorder Factor levels
gas_meanse_melted$variable  <- factor(gas_meanse_melted$variable, levels = c("total", "product.ratio", "nitox", "dinitrogen", "nitox.pool"))
gas_meanse_melted$treat     <- factor(gas_meanse_melted$treat, levels = c("Lt", "Fc", "LF", "C"))
# For facet.grid:
levels(gas_meanse_melted$variable)[levels(gas_meanse_melted$variable)=="total"] <- "paste(N[2],O)+N[2]"
levels(gas_meanse_melted$variable)[levels(gas_meanse_melted$variable)=="product.ratio"] <- "Product~Ratio"
levels(gas_meanse_melted$variable)[levels(gas_meanse_melted$variable)=="nitox"] <- "paste(N[2],O)"

#Coerce objects in required structure, 
# keep an eye on as.numeric(as.character()), for coercion from factor to numeric!!!
gas_meanse_melted$value     <- as.numeric(gas_meanse_melted$value)
gas_meanse_melted$se        <- as.numeric(as.character(gas_meanse_melted$se))
str(gas_meanse_melted)

#**** Multi Panel Plot ****
#________________________________________________________________________________________________________________________________

txt.total = expression(paste(µg~Total-N~h^-1,kg^-1))
txt.ratio = expression(paste(N[2],O)/paste(N[2],O+N[2]))
txt.nitox = expression(paste(µg~paste(N[2],O)-N~h^-1,kg^-1))
txt.full = expression(paste(paste(µg~paste(N[2],O)-N~kg^-1,h^-1),"       ",frac(paste(N[2],O),paste(N[2],O+N[2])),"        ",paste(µg~Total-N~kg^-1,h^-1)))

pd=position_dodge(0.2)
ggplot(data=gas_meanse_melted, aes(x=day, y=value, group=treat, shape=treat)) +
  facet_grid(variable~soil, scales="free", labeller=label_parsed) +
  #guides(colour=F, shape=F) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), size=0.25, width=.3, position=pd) + 
  geom_line(size=0.3, position=pd) +
  geom_point(size=1.5, position=pd, colour="black", fill="white") +
  scale_shape_manual(values=c(17,16,15,21), labels=c("L. terrestris", "F. candida","Interaction","Control")) +
  #scale_colour_manual(values=trt_pal) +
  #scale_colour_grey(start=0, end=0.6) +
  labs(shape="") +
  xlab("Days [d]") +
  ylab(txt.full) +
  theme_bw() + theme(strip.background = element_rect(fill="black"),
                     strip.text.x = element_text(size=8, face="bold", colour="white"),
                     strip.text.y = element_text(size=8, face="bold", colour="white"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.margin=unit(0,"lines"),
                     panel.border = element_rect(colour="black", size=0.2, fill=NA),  
                     axis.text.x = element_text(size=7),
                     axis.title.x = element_text(size=8,face="bold", family="Times New Roman"),
                     axis.text.y = element_text(size=7),
                     axis.title.y = element_text(size=8,face="bold", family="Times New Roman"),
                     axis.ticks = element_line(size=0.25),
                     legend.position="bottom",
                     legend.key=element_blank(),
                     legend.background=element_blank(),
                     legend.text=element_text(size=8,face="italic", family="Times New Roman"),
                     legend.title=element_text(size=8))


ggsave("Figure3a.pdf",  width=9, height=10.5, unit="cm", useDingbats=FALSE)                             

#### nitox.pool ####
#________________________________________________________________________________________________________________________________

# create dataframe #
# Mean values over MK
gas_meanse_pool <- matrix(0,26,7)
gas_meanse_pool <- as.data.frame(gas_meanse_pool)
gas_meanse_pool[,1] = gas15[1:26,1]
gas_meanse_pool[,2] = factor(rep(c("Lt", "LF", "Fc", "C"), c(4,3,3,3)),levels = c("Lt", "Fc", "LF", "C"))
gas_meanse_pool[,3] = factor(rep(c("Loam", "Sand"), c(13,13)))
gas_meanse_pool[,4] = with(gas15, round(tapply(nitox, MK, mean),3))
gas_meanse_pool[,5] = with(gas15, round(tapply(nitox.pool, MK, mean),3))
gas_meanse_pool[,6] = with(gas15, round(tapply(nitox, MK, function(x) sqrt(var(x)/length(x))),3))
gas_meanse_pool[,7] = with(gas15, round(tapply(nitox.pool, MK, function(x) sqrt(var(x)/length(x))),3))
colnames(gas_meanse_pool) = head2 = c("MK","treat","soil","nitox.mean","nitox.pool.mean","nitox.se","nitox.pool.se")

# gas_meanse_pool$treat <- factor(gas_meanse_pool$treat, levels = c("Lt", "Fc", "LF", "C"))
gas_meanse_pool$Lter <- factor(rep(rep(c("+Lt", "-Lt"), c(7,6)),2))
gas_meanse_pool$col.fac <- interaction(gas_meanse_pool$soil, gas_meanse_pool$Lter, sep="")
str(gas_meanse_pool)

# create Single Panel Plot with ellipse #
#________________________________________________________________________________________________________________________________

# Color palettes
LtSoil_pal=c("gray8", "gray30", "gray48", "gray70")
LtSoil_pal2=c("darkseagreen", "burlywood", "darkolivegreen3", "darkgoldenrod1")
title = expression(paste("Plant derived ",N[2],"O depending on soiltype and presence of L. terrestris"))

# Ellipse Building
# Script from http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo

# For annotation
ellipse.mean=with(gas_meanse_pool, aggregate(gas_meanse_pool[,c(4,5)],list(group=col.fac),mean))

# Special function 
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

# Data frame df_ell contains values to show ellipses. 
# It is calculated with function veganCovEllipse which is hidden in vegan package. 
# This function is applied to each level of NMDS (group) and 
# it uses also function cov.wt to calculate covariance matrix.

df_ell <- data.frame()
for(g in levels(gas_meanse_pool$col.fac)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(gas_meanse_pool[gas_meanse_pool$col.fac==g,],
                                                   veganCovEllipse(cov.wt(cbind(nitox.mean,nitox.pool.mean),wt=rep(1/length(nitox.mean),length(nitox.mean)))$cov,center=c(mean(nitox.mean),mean(nitox.pool.mean)))))
                                ,group=g))
}

# plot #
#________________________________________________________________________________________________________________________________

txt.pool = expression(paste("Plant derived ", N[2],"O","[%]"))
txt.nitox = expression(paste("Nitrous Oxide ", N[2],"O [µg ",N[2],"O-N ",kg^-1, h^-1,"]"))

p4 = ggplot(data=gas_meanse_pool, aes(x=nitox.mean, y=nitox.pool.mean)) +
  #guides(colour=F, shape=F) +
  geom_point(aes(colour=col.fac, shape=treat), size=1.5) +
  geom_path(data=df_ell, aes(x=nitox.mean, y=nitox.pool.mean, colour=group), size=0.5, linetype=2) +
  annotate("text",x=ellipse.mean$nitox.mean+0.1,y=ellipse.mean$nitox.pool.mean,label=ellipse.mean$group, size=2.5) +
  scale_shape_manual(values=c(17,16,15,21)) +
  scale_colour_manual(values=LtSoil_pal) +
  scale_y_continuous(limits=c(0,1)) +
  ylab(txt.pool) +
  xlab(txt.nitox)
p4 +  mytheme +theme(
  strip.background = element_rect(fill="black"),
  strip.text.x = element_text(size=8, face="bold", colour="white"),
  strip.text.y = element_text(size=8, face="bold", colour="white"),
  panel.margin = unit(0, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(colour="black", size=0.2, fill=NA),
  axis.text.x = element_text(size=7),
  axis.title.x = element_text(size=8,face="bold", family="Times New Roman"),
  axis.text.y = element_text(size=7),
  axis.title.y = element_text(size=8,face="bold", family="Times New Roman"),
  axis.ticks = element_line(size=0.25),
  legend.key=element_blank(),
  legend.background=element_blank(),
  legend.text=element_text(size=8,face="italic", family="Times New Roman"),
  legend.title=element_text(size=8),
  legend.position="none")
#legend.position=c(0.08,0.75) )

ggsave("Figure 4a.pdf", width=9, height=7, units="cm", useDingbats=FALSE)
ggsave("Figure 4a.svg", width=9, height=7, units="cm")
#________________________________________________________________________________________________________________________________

# Rudimentary:
# Mittelwerte über Zeit, treat und soil
gas_meanse_pool2       <- cbind(gas_mean2[,1:3],gas_mean2[,c(6,8)], gas_se2[,c(6,8)])
colnames(gas_meanse_pool2)[c(4:7)] <- c("nitox.mean","nitox.pool.mean", "nitox.se", "nitox.pool.se")
gas_meanse_pool2$nitox.mean  <- round(as.numeric(as.character(gas_meanse_pool$nitox.mean)),3)
gas_meanse_pool2$nitox.pool.mean  <- round(as.numeric(as.character(gas_meanse_pool2$nitox.pool.mean)),3)
gas_meanse_pool2$nitox.se    <- round(as.numeric(as.character(gas_meanse_pool2$nitox.se)),3)
gas_meanse_pool2$nitox.pool.se    <- round(as.numeric(as.character(gas_meanse_pool2$nitox.pool.se)),3)
gas_meanse_pool2$day   <- as.integer(as.character(gas_meanse_pool2$day))
gas_meanse_pool2$treat <- factor(gas_meanse_pool2$treat, levels = c("Lt", "Fc", "LF", "C"))
gas_meanse_pool2$col.fac <- interaction(gas_meanse_pool2$treat, gas_meanse_pool2$soil)
str(gas_meanse_pool2)











########## Analysis ##########
#________________________________________________________________________________________________________________________________

### Cumulative N2O, N2, Total

setwd("D:/Quentin Schorpp/Arbeitsprozess/Mikrokosmosversuch/15N_gas/Analyse_02.04.14")
gas15 <- read.delim("D:/Quentin Schorpp/Schreibtisch/2 - Laborversuche/Mikrokosmen/Daten/Isotope 15N-13C/15N Gas/V1_15Ngas.txt")

library(reshape2)
library(ggplot2)

# subsets
treatments = gas15[,8:13]
time = gas15[,4:7]
gas = gas15[,19:24]


# data-frame mit R 
# Cumulative gas between two dates
cumulative = function (x) 168*(x[i:(length(x)-j)]+0.5*(x[(i+j):length(x)]-x[i:(length(x)-j)]))
gas15.cum = matrix(rep(0,78),,3)
i=1
j=26
a = cumulative(gas15$total)
b = cumulative(gas15$dinitrogen)
c = cumulative(gas15$nitox)

# summed up for up to last date
cumulative2 = function(x) x[i:j]+x[(i+26):(j+26)]+x[(i+26+26):(j+26+26)]
gas15.cum[,1] = a = cumulative2(a)
gas15.cum[,2] = b = cumulative2(b)
gas15.cum[,3] = c = cumulative2(c)
colnames(gas15.cum) = colnames(gas15[,22:24])


# data frame for analysis 
gas15.cum = cbind(gas15[1:26,1:3],gas15[1:26,8:13], gas15.cum)
str(gas15.cum)

write.table(gas15.cum, "gas15_cumR.txt", dec=".", sep=";", row.names=F)

#________________________________________________________________________________________________________________________________

# Second dataframe for standard error
cumulative = function (x) 168*(x[i:(length(x)-j)]+0.5*(x[(i+j):length(x)]-x[i:(length(x)-j)]))
gas15.cum2 = matrix(0,78,3)
i=1
j=26
gas15.cum2[,1] = a = cumulative(gas15$total)
gas15.cum2[,2] = b = cumulative(gas15$dinitrogen)
gas15.cum2[,3] = c = cumulative(gas15$nitox)
colnames(gas15.cum2) = colnames(gas15[,22:24])
gas15.cum3 = gas15.cum2
gas15.cum3 = as.data.frame(gas15.cum3, stringsAsFactor=F)
str(gas15.cum3)

for (j in 1:3)  {
  for (i in 27:78) {
    gas15.cum2[i,j] =  gas15.cum2[i,j] + gas15.cum2[i-26,j]
  }
}

gas15.cum2 = cbind(gas15[27:104,1:13], gas15.cum2)
gas15.cum2$treat = factor(gas15.cum2$treat, level=c("Lt","Fc","LF","C"))
str(gas15.cum2)
# gas_cum_test= gas15.cum2[gas15.cum2$day == 40,]


# df for mean over time(day)
gas15.cum2.mean   = matrix(0,24,6)
gas15.cum2.mean[,1] = rep(c("Lt", "Fc", "LF", "C"),c(3,3,3,3))
gas15.cum2.mean[,2] = rep(c(26, 33, 40),4)
gas15.cum2.mean[,3] = rep(c("Loam", "Sand"),c(12,12))
gas15.cum2.mean[,4] = with(gas15.cum2, tapply(total, list(day,treat,soil), mean))
gas15.cum2.mean[,5] = with(gas15.cum2, tapply(nitox, list(day,treat,soil), mean))
gas15.cum2.mean[,6] = with(gas15.cum2, tapply(dinitrogen, list(day,treat,soil), mean))
colnames(gas15.cum2.mean) = c("treat","day","soil","total.mean","nitox.mean","dinitrogen.mean")
gas15.cum2.mean <- as.data.frame(gas15.cum2.mean, stringsAsFactors=F)
str(gas15.cum2.mean)


# df for standard error over time(day)
gas15.cum2.se     = matrix(0,24,6)
gas15.cum2.se[,1] = rep(c("Lt", "Fc", "LF", "C"),c(3,3,3,3))
gas15.cum2.se[,2] = rep(c(26, 33, 40),4)
gas15.cum2.se[,3] = rep(c("Loam", "Sand"),c(12,12))
gas15.cum2.se[,4] = with(gas15.cum2, tapply(total, list(day,treat,soil), function(x) sqrt(var(x)/length(x))))
gas15.cum2.se[,5] = with(gas15.cum2, tapply(nitox, list(day,treat,soil), function(x) sqrt(var(x)/length(x))))
gas15.cum2.se[,6] = with(gas15.cum2, tapply(dinitrogen, list(day,treat,soil), function(x) sqrt(var(x)/length(x))))
colnames(gas15.cum2.se) = c("treat","day","soil","total.se","nitox.se","dinitrogen.se")
gas15.cum2.se <- as.data.frame(gas15.cum2.se, stringsAsFactors=F)
str(gas15.cum2.se)

# Formel für Fehler durch Fehlerfortpflanzung: function (x) sqrt(sum((Kum_StdF[i,1:j])^2

gas15.cum2.se <- melt(gas15.cum2.se, id.vars= c(1:3))
gas15.cum2.se$treat = factor(gas15.cum2.se$treat, level=c("Lt","Fc","LF","C"))
gas15.cum2.se <- dcast(gas15.cum2.se, variable + soil + treat ~ day)
colnames(gas15.cum2.se) [1] = "gas"
gas15.cum2.se.ts <- as.data.frame(lapply(gas15.cum2.se[,4:6], as.numeric))

for (i in 1:24)  {
  for (j in 1:3) {
    a=sqrt(sum((gas15.cum2.se.ts[i,1:j])^2))
    gas15.cum2.se.ts[i,j]=a
  }
}

gas15.cum2.se.ts <- melt(cbind(gas15.cum2.se[,1:3], gas15.cum2.se.ts), id.vars=c(1:3))
gas15.cum2.se.ts <- dcast(gas15.cum2.se.ts, soil+treat+variable~gas)

gas15.cum2.meanse <- cbind(gas15.cum2.mean, gas15.cum2.se.ts[,4:6])
write.table(gas15.cum2.meanse, "gas15_cum2_meanse.txt", dec=".", sep=";", row.names=F)


# Plots
helpframe1 = melt(gas15.cum2.meanse, id.vars=c(1:3), measure.vars=c(4:6), variable.name="mean.fac", value.name="mean")
helpframe2 = melt(gas15.cum2.meanse, id.vars=c(1:3), measure.vars=c(7:9), variable.name="se.fac", value.name="se")
gas15.cum2.meanse.plot <- cbind(helpframe1,helpframe2[,c(4,5)])
gas15.cum2.meanse.plot$treat = factor(gas15.cum2.meanse.plot$treat, levels=c("Lt","Fc","LF","C"))
gas15.cum2.meanse.plot$mean = as.numeric(gas15.cum2.meanse.plot$mean)
str(gas15.cum2.meanse.plot)

pd=position_dodge(0.3)
p5 = ggplot(data=gas15.cum2.meanse.plot, aes(x=day, y=mean, group=treat, colour=treat, shape=treat)) +
  facet_grid(mean.fac~soil, scales="free") +
  #guides(colour=F, shape=F) +
  geom_point(size=2.5, position=pd) +
  geom_line(position=pd) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.3, position=pd) + 
  scale_shape_manual(values=c(15,16, 18,17)) +
  #scale_colour_manual(values=trt_pal) +
  #scale_colour_grey(start=0, end=0.6) +
  #scale_x_continuous(name="Days [d]", breaks=gas15.cum_melted$day) +
  ylab("")
p5 + mytheme + theme(legend.position="none")



# data-frame mit Excel
# gas15.cum2 <- read.delim("D:/Quentin Schorpp/Schreibtisch/2 - Laborversuche/Mikrokosmen/Daten/Isotope 15N-13C/15N Gas/V1_15Ngas_cumulative.txt")

# Tabelle mit Mittelwerten und Standardabweichung für Cumulative Emmission
# Zur Veranschaulichung

gas15.cum$troil = interaction(gas15.cum$treat,gas15.cum$soil)


gas15.cum2.mean = matrix(rep(0,24),8)
for (j in 10:12)
{a=tapply(gas15.cum[,j],gas15.cum[,13],mean)
 gas15.cum2.mean[,(j-9)]=a                        # mit ..cum[,1:3] geht es nicht!!!!!
}
gas15.cum2.mean <- as.data.frame(gas15.cum2.mean)
colnames(gas15.cum2.mean)<- c("total.mean","nitox.pool.mean", "nitox.mean")
gas15.cum2.mean$treat  <- factor(c("Lt", "Fc", "LF", "C"), levels = c("Lt", "Fc", "LF", "C"))
gas15.cum2.mean$soil  <- factor(rep(c("Loam", "Sand"), c(4,4)))
gas15.cum2.mean$troil = interaction(gas15.cum2.mean$treat,gas15.cum2.mean$soil)


gas15.cum2.se=matrix(rep(0,24),8) 
colnames(gas15.cum2.se)<- c("total.mean","nitox.pool.mean", "nitox.mean")
for (j in 10:12)
{a=tapply(gas15.cum[,j],gas15.cum[,13],function(x) sqrt(var(x)/length(x)))
 gas15.cum2.se[,(j-9)]=a
}

gas_se_cum2=matrix(rep(0,24),8)   #Matrix für stdfehler fortpflanzung kum n20 über die zeit
colnames(gas_se_cum)<- c("total.mean","nitox.pool.mean", "nitox.mean")
for (i in 1:16)
{for (j in 1:14)
{a=sqrt(sum((gas_se_cumF[i,1:j])^2)) 
 Kum_Zeit_Fehler[i,j]=a
}
}

#### Statistical Analysis #####

# response = total
ano1 <- with(gas15.cum, aov(total~Lt*Fc*soil))
summary(ano1)
par(mfrow=c(2,2))
plot(ano1)

# response = nitox
ano2 <- with(gas15.cum, aov((nitox~Lt*Fc*soil))
             summary(ano2)
             
             # response = dinitrogen
             ano3 <- with(gas15.cum, aov(dinitrogen~Lt*Fc*soil))
             summary(ano3)
             
             
             #### Pseudoreplication averaged away ####
             # Analysis with mean production rates in µg N / kg h
             
             gas_mean=matrix(rep(0,130),26)   # create a matrix with four columns for the four response variables
             # fill in columns with means of response variables
             gas_mean[,1]= tapply(gas15$total, gas15$MK, mean)           
             gas_mean[,2]= tapply(gas15$dinitrogen, gas15$MK, mean)
             gas_mean[,3]= tapply(gas15$nitox, gas15$MK, mean)
             gas_mean[,4]= tapply(gas15$product.ratio, gas15$MK, mean)
             gas_mean[,5]= tapply(gas15$nitox.pool, gas15$MK, mean)
             
             ## change column/variable names,.... 
             mean_names = c("total", "dinitrogen", "nitox", "product.ratio", "nitox.pool")
             colnames(gas_mean) = mean_names   # give names to responses
             # ....or this way,....
             colnames(gas_mean)[4:5] = colnames(gas15[,20:21])   # give names to responses
             colnames(gas_mean)[1:3] = colnames(gas15[,22:24])   # give names to responses
             
             # bind factors to df
             gas_mean=cbind(ident[1:26,],treatments[1:26,],gas_mean)        
             gas_mean$troil <- with(gas_mean, interaction(soil, treat)) # create new (full) factor
             write.csv(gas_mean, "Mean_values.csv")
             
             
             
             with(gas_mean, hist(product.ratio))
             with(gas_mean, boxplot(product.ratio~Lt*soil))
             
             
             # Analysis response = ratio
             ano4 <- with(gas_mean, aov(product.ratio~Lt*Fc*soil))
             summary(ano4)
             TukeyHSD(ano4)
             
             library(nlme)             
             ano4.1 = glm(product.ratio ~ Lt*Fc*soil, family = "quasibinomial", data=gas_mean)
             anova(ano4.1)             
             summary(ano4.1)             
             
             ano4.2 = aov(product.ratio~Lt*Fc, gas15.cum[gas15.cum$soil=="Sand",])
             summary(ano4.2)
             
             # Analysis response = pool
             ano5 <- with(gas_mean, aov(nitox.pool~Lt*Fc*soil))
             summary(ano5);
             
             boxplot(product.ratio~Lt*soil, data=gas_mean)
             boxplot(product.ratio~Lt, data=gas_mean)
             boxplot(product.ratio~soil, data=gas_mean)
             
             ### Non parametric tests ###
             
             nonpm1 <-  kruskal.test(product.ratio~Lt*Fc*soil, data=gas_mean)
             summary(nonpm1)             
             
             
             with(gas15, hist(sqrt(product.ratio)))
             with(gas15, cor.test(product.ratio,week))
             with(gas15, cor.test(product.ratio,week, method="kendall"))
             with(gas15, cor.test(product.ratio,week, method="spearman"))
             
             
             ### Rudimentary Analysis, Why did i use the gas.cum2 df here? in this all data for course over time is included
             
             # Analysis response = total
             ano1 <- with(gas15.cum2, aov(total~Lt*Fc*soil))
             ano1.1 <- with(gas15.cum2, aov(total~lit_pow*lit_chop*soil))
             summary(ano1); summary(ano1.1)
             anova(ano1, ano1.1)
             AIC(ano1, ano1.1)     # ano1.1 ist schlechter
             
             
             
             # Output
             out <- summary(ano1)
             
             write.csv2(out, file="out.csv")
             
             out2<-capture.output(summary(ano1))
             , file="out.csv")
cat(out,file="out.txt", sep="\n", append=T)

(summary(ano1))
out<-capture.output(vcov(ano1))
cat(out,file="out2.txt",sep="\n",append=TRUE)

par(mfrow=c(1,1))
boxplot(total~Lt)
interaction.plot(Lt, soil, total)

### Regarding pseudoreplication over time ###
# Temporal pseudoreplication, Crawley, The R-Bok p. 695

# 2 -----------------------------------------------------------------------


hist(gas15$total, breaks=20) # Überlegunf: Poisson modell mit lme4
hist(sqrt(gas15$total), breaks=20)
hist(log(gas15$total+1), breaks=20) # +1 da Nullen im Datensatz

library(nlme)
library(lattice)
library(AICcmodavg)
str(gas15)

# create factors out of MK and week
gas15$fMK <- as.factor(gas15$MK)
gas15$fweek <- as.factor(gas15$week)

# 1:1 aus Crawley transponiert, funktioniert nur mit control function
gas2 <- groupedData(total~week|MK, outer=~soil*Lt*Fc,gas)
str(gas2)
ctrl <- lmeControl(opt='optim');
model <-with(gas2, lme(total~Lt*Fc*soil, random=~week|MK, control=ctrl))

# Modellverbeserung (mit Doreen)

model <-lme(sqrt(total)~Lt*Fc*soil*fweek, random=~1|MK, data=gas15)


# Alternativmodelle, die aber alle schlechter sind
model1 <-lme(sqrt(total)~Lt*Fc*soil, random=~1|week/MK, data=gas) # plot(model) = Kurvenverlauf
model2 <-lme(log(total+1)~Lt*Fc*soil*fweek, random=~1|MK, data=gas)
model3 <-lme(sqrt(total)~Lt*Fc*soil*fweek, random=~1|fweek/MK, data=gas)
model3 <-lme(log(total+1)~Lt*Fc*soil, random=~1|MK, subset=week==3, data=gas)



plot(model) # sieht gut aus, weitere plots: Residuen gegen Messwerte,...
summary(model)
ranef(model) # random effects des models
anova(model, type="marginal")  # Zu beginn der Modellvereinfachung stellen wir fest, dass die vierfach Interaktion signifikant ist!!
# Eine Weitere Vereinfachung ist nicht nötig!!

# Create graphics
# barplots mit rücktransformierten prediction values und den dazu gehörigen rücktransformierten Fehlern

td=expand.grid(Lt=levels(gas15$Lt), Fc=levels(gas15$Fc), soil=levels(gas15$soil), fweek=levels(gas15$fweek))
td$p=predict(model, td, level=0)
td$se=predictSE.lme(model,newdata=td, se.fit=T, level=0)$se.fit
td

boxplot(p~Lt*Fc*soil*fweek, data=td)
xyplot(p^2~fweek|Lt*Fc*soil, data=td) # in etwa dieser plot als barplot mit standard error!
# function barplots2 aus gplots
# p^2 weil total in analyse sqrt tranformiert!

xyplot(p~fweek|Lt*Fc*soil, data=td)




#### Pseudoreplication averaged away ####
# Analysis with mean production rates in µg N / kg h

gas_mean=matrix(rep(0,130),26)   # create a matrix with four columns for the four response variables
# fill in columns with means of response variables
gas_mean[,1]= tapply(gas15$total, gas15$MK, mean)           
gas_mean[,2]= tapply(gas15$dinitrogen, gas15$MK, mean)
gas_mean[,3]= tapply(gas15$nitox, gas15$MK, mean)
gas_mean[,4]= tapply(gas15$product.ratio, gas15$MK, mean)
gas_mean[,5]= tapply(gas15$nitox.pool, gas15$MK, mean)

## change column/variable names,.... 
mean_names = c("total", "dinitrogen", "nitox", "product.ratio", "nitox.pool")
colnames(gas_mean) = mean_names   # give names to responses
# ....or this way,....
colnames(gas_mean)[4:5] = colnames(gas15[,20:21])   # give names to responses
colnames(gas_mean)[1:3] = colnames(gas15[,22:24])   # give names to responses

# bind factors to df
gas_mean=cbind(ident[1:26,],treatments[1:26,],gas_mean)        
gas_mean$troil <- with(gas_mean, interaction(soil, treat)) # create new (full) factor
write.csv(gas_mean, "Mean_values.csv")

# gas[MK!=1,]  would Use dataframe without MK1

# Analysis response = total

#global model: Unbalanciertes design, da Futterqualität nicht in ausreichender Wiederholung vorhanden

aov1.1 <- aov(total~Lt*Fc*soil, data=gas_mean)
aov1.2 <- aov(total~lit_pow*lit_chop*soil, data=gas_mean)
summary(aov1.1); summary(aov1.2)
anova(aov1.1, aov1.2)
AIC(aov1.1, aov1.2) # aov1.2 ist deutlich schlechter

# Output
total.out<-capture.output(summary(aov1.1), file="total_mean_out3.csv")

total.out<-capture.output(summary(aov1.1))
cat(total.out,file="total_mean_out.txt",sep="\n",append=TRUE)
total.out<-capture.output(vcov(ano1))
cat(total.out,file="total_mean_out.txt",sep="\n",append=TRUE)

# Graphics
library(car)
par(mfrow=c(1,1))
Boxplot(total~Lt, data=gas_mean, las=1)
Boxplot(total~soil, data=gas_mean, las=1)

library(ggplot2)

gas_mean.plot <- qplot(Lt, log(total), col=soil, data=gas_mean) + geom_point(size=2.5) +
  theme(
    panel.background= element_rect(fill="white"),
    panel.grid.major = element_line(colour="grey"),
    panel.grid.minor = element_line(colour="grey", linetype="dashed", size=0.2),
    panel.border= element_rect(colour="black", fill=NA, size=1),
    axis.text = element_text(colour="black", size=rel(1.3)),
    axis.text.x = element_text(face="italic"),
    axis.title=element_text(size=rel(1.4)),
    legend.text=element_text(face="italic"))
gas_mean.plot

gas_mean.plot + xlab("Soiltype") + ylab(expression(paste(.^15,"N atom%",sep=""))) + ggtitle(expression(paste("Label in ",italic("L. terrestris"))))+
  scale_x_discrete(breaks=c("L","S"), labels=c("Loam", "Sand"))+
  labs(colour="Treatment")+
  scale_colour_discrete(labels=c("L.terrestris", "L. terrestris + F. candida"))

with(gas_mean, interaction.plot(Lt, soil, total)
     
     
     
     
     
     
     
     
     # Analysis response = ratio
     ano2 <- with(gas_mean, aov(product.ratio~Lt*Fc*soil))
     ano2.1 <- with(gas_mean, aov(ratio~soil*food))
     summary(ano2); summary(ano2.1)
     anova(ano2, ano2.1)
     AIC(ano2, ano2.1)
     # Output
     out<-capture.output(summary(ano2))
     cat(out,file="out9.txt",sep="\n",append=TRUE)
     out<-capture.output(vcov(ano2))
     cat(out,file="out4.txt",sep="\n",append=TRUE)
     
     # Analysis response = nitox
     ano3 <- with(gas_mean, aov(nitox~Lt*Fc*soil))
     ano3.1 <- with(gas_mean, aov(nitox~food*soil))
     summary(ano3); summary(ano3.1)
     AIC(ano3, ano3.1)
     # Output
     out<-capture.output(summary(ano3))
     cat(out,file="out5.txt",sep="\n",append=TRUE)
     out<-capture.output(vcov(ano3))
     cat(out,file="out6.txt",sep="\n",append=TRUE)
     
     # Analysis response = pool
     ano4 <- with(gas_mean, aov(nitox.pool~Lt*Fc*soil))
     ano4.1 <- with(gas_mean, aov(pool~food*soil))
     summary(ano4); summary(ano4.1)
     anova(ano4,ano4.1)
     AIC(ano4, ano4.1)
     # Output
     out<-capture.output(summary(ano4))
     cat(out,file="out7.txt",sep="\n",append=TRUE)
     out<-capture.output(vcov(ano4))
     cat(out,file="out8.txt",sep="\n",append=TRUE)
     
     par(mfrow=c(1,1))
     boxplot(ratio~soil*food)
     boxplot(ratio~Lt*Fc*soil)