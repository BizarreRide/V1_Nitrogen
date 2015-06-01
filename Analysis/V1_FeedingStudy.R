############################
# V1_Nitrogen
# Feeding study
# Quentin Schorpp
# 23.05.2015
###########################

source("Data/GatherSource/V1_MakeLikeFile.R")

# 1. Eathworm Biomass ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# earthworm biomass was measured in g after earthworm guts were voided. Mesurements took place before and after the experiment.
# During the Experiment four Individual of L. terrestris died in Microcosms with the interaction Treatment in the 15N experimental part / approach.

## Survival Rates: ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 100*(1-dead worms/all worms)
Exp1  = 100*(1-4/(4*2*2*2)) # = 87.5 %
Exp2  = 100*(1-0/(3*2*2*2)) # = 100 %
Total = 100*(1-4/((4*2*2*2)+(3*2*2*2))) # = 92.86 %
Lt.SurvivalRates <- rbind(Exp1,Exp2,Total)
colnames(Lt.SurvivalRates)[1] <- "%"
rm(Exp1,Exp2,Total)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## Earthworm Biomass Data Exploration ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Ranges
with(ew.bm, list(length(before), range(before),range(after)))
# create binsizes for Histograms
binsize1 = diff(range(ew.bm$before))/26
binsize2 = diff(range(ew.bm$after))/26
binsize3 = diff(range(ew.bm$diff))/length(ew.bm$diff)

# Histograms Before/After
ew.bm.melt = melt(ew.bm[,-c(7,8)], id.vars=c(1:4))
ew.bm.p1 <- ggplot(ew.bm.melt, aes(x=value)) + 
  geom_histogram( binwidth = binsize1, fill="lightblue3", colour="black") + 
  geom_line(stat="density") +
  xlim(6.5,11) + 
  xlab("Biomass [g]") +
  ylab("Frequency") +
  ggtitle("Histograms EW-biomass \n Before/After") +
  facet_grid(variable ~ .) + mytheme

# Histogram of differences
ew.bm.p2 <- ggplot(ew.bm, aes(x=diff)) +
  geom_histogram(binwidth = binsize3, fill="lightblue3", colour="black") + 
  geom_line(stat="density") + 
  xlab("differences [g]") + 
  ylab("Frequency") + 
  ggtitle("Histogram of differences") + mytheme

# One figure in row 1 and two figures in row 2
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(ew.bm.p1, vp = vplayout(1, 1))
print(ew.bm.p2, vp = vplayout(1, 2))

rm(ew.bm.melt, binsize1, binsize2, binsize3)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



## Earthworm Biomass Analysis ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shapiro.test(ew.bm$diff) # kann bei diesem Test eine Signifikanz (p < 0.05) festgestellt werden, so liegt keine Normalverteilung vor.

qqnorm(ew.bm$diff)
qqline(ew.bm$diff)


# t.test ####
with(ew.bm, t.test(after, before, paired=T, alternative="less")) 
# Defining alt="less", means check whether the mean of the values contained in the vector a is less 
# of the mean of the values contained in the vector b 

# ANOVA ####
# Were Changes in L. terrestris Biomass affected by Treatments?
summary(aov(diff ~ treat*soil*exp, data=ew.bm))
# NO it was not.

list(mean=apply(ew.bm[,5:6],2,mean),
     SE=apply(ew.bm[,5:6],2,se))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# 2. Collembolan Population growth ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Collembolans grew very well during the experiment and 
# reached densities hardly observed in agricultural systems.

## Growth Factor ####
initial = 238 #Ind/MC
# How many times more collembolans afterwards?
col.pop$growth = col.pop$surface/initial

mean(col.pop$growth); se(col.pop$growth);range(col.pop$growth)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## Collembolan Population Data Exploration ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

col.pop$exp <- revalue(col.pop$exp, c(Exp1="Experiment~1",Exp2="Experiment~2"))
col.pop$treat <- revalue(col.pop$treat, c(C ="F.~candida",RC="Interaction"))

ggplot(col.pop, aes(x= treat, y=surface),labeller=label_parsed) +
  stat_boxplot(geom="errorbar", coef=1.5, lwd=0.2) + 
  geom_boxplot(fill="light blue",lwd=0.2, outlier.size=0.4, outlier.shape=21) +
  ylab(expression(paste("Collembolan density [Ind ",MC[surface]^-1,"]"))) +
  xlab("Treatment") +
  ggtitle(expression(paste("Boxplots for collembolan density per ",MC[surface]))) +
  facet_grid(exp~soil, labeller=label_parsed) + 
  scale_x_discrete("Treatment", labels=expression(italic(F.~candida), Interaction)) +
  mytheme

# In 3 of 4 cases Presence of Lumbricus terrestris led to increased Collembolan population growth. 
# The only exception is the loamy soil treatment in Exp1. Here the opposite was observed.
# Furthermore loamy soil treatments increased collembolan population growth in both experiments stronger compared to sandy soil (`r tapply(col.pop$surface, col.pop$soil,mean)`, Loam, Sand).
# Hypothesis: Low Nitrogen content of food material leads to decreased collembolan populations in loamy soil if L. terrestris was present.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## Collembolan populations Analysis ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Table of means
with(col.pop, tapply(surface, list(soil, treat, exp), mean))

#### Analysis of variance ####
cp.lm1 <- lm(surface ~ treat*soil*exp, data=col.pop);cp.lm1; summary.aov(cp.lm1);summary.lm(cp.lm1)
cp.lm2 <- lm(surface ~ treat*soil*exp-exp:treat:soil, data=col.pop);cp.lm2; summary.aov(cp.lm2);summary.lm(cp.lm2)
cp.aov <- aov(surface ~ treat*soil*exp-exp:treat:soil, data=col.pop)


# The coefficients table of summary.lm has as many rows, as there are parameters in the model.
# The Intercept is the only mean value in the table. The second column contains the unreliability estimates. 
# The first row contains the standard error of a MEAN. 
# The Other four rows contain the standard errors of  the DIFFERENCE between two means.
# The p-values are misleading, suggesting, wrongly, that there are four significant contrasts for this model.
### All significant differences are differences to Exp1-Loam-F.candida

# Effect Sizes of Main Effects
plot.design(surface ~ soil*exp*treat-exp:treat:soil, data=col.pop)

# The populations were larger in loamy soil
with(col.pop, list(mean=tapply(surface, soil, mean),
                   SE=tapply(surface, soil, se)))

# The populations were larger in Experiment 2
with(col.pop, list(mean=tapply(surface, exp, mean),
                   SE=tapply(surface, exp, se)))

# ANOVA plots
par(mfrow=c(2,2))
par(mar=c(3,4,1.5,2), mgp=c(1.5,0.5,0))
plot(cp.lm2)


### Post Hoc TukeyHSD ####

# Variant 1:
cp.int <- with(col.pop, interaction(exp,soil,treat)) # Triple interaction factor
col.aov2 <- aov(surface ~ cp.int, data=col.pop)
HSD.test(col.aov2, "cp.int", group=TRUE, console=TRUE)

#Variant 2:  
col.tuk1 <- glht(col.aov2, linfct = mcp(cp.int = "Tukey"))
summary(col.tuk1)          # standard display
col.tuk.cld <- cld(col.tuk1)   # letter-based display
par(mai=c(1,1,1.1,0.2),
    mfrow=c(1,1),
    mgp=c(2,1,0),
    cex=0.8) # mai specifiec margin size in inches
plot(col.tuk.cld, cex=0.5, las=2, col="grey", xaxt="n")
axis(1, at=cp.int,labels=FALSE)
text(cp.int, labels=cp.int, par("usr")[3], adj=c(1.2,1.2), xpd=TRUE, srt=45, cex=0.8)


# Effect sizes of interactions
col.pop.effects <- allEffects(cp.lm2)

plot(col.pop.effects, "treat:soil")
plot(col.pop.effects, "treat:exp")
plot(col.pop.effects, "soil:exp")

col.pop.effects <- allEffects(cp.lm1)
plot(col.pop.effects, "treat:soil:exp")

model.tables(cp.aov, "means", se=TRUE)
with(col.pop, tapply(surface, list(treat,soil), mean))
with(col.pop, tapply(surface, list(treat,soil), se))

with(col.pop, tapply(surface, list(treat,exp), mean))
with(col.pop, tapply(surface, list(treat,exp), se))

with(col.pop, tapply(surface, list(soil,exp), mean))
with(col.pop, tapply(surface, list(soil,exp), se))


# Interessieren w?rde mich wie viel gr??er die Populationen im Schnitt in Anwesenheit von L. terrestris waren.
# Datensatz zu differenzen zw. Interaktion und Single-Species Treatments
# PROBLEM: Welche Mikrokosmen bilden Paare? 
# Bilden aller paarweise DIfferenzen und/oder nur f?r selbe Faktoren Kombination aus Exp und soil?
# Wie subtrahiert man Standardfehler voneinander (gibt T-Test Aufschl?sse?)


### Tukey HSD plots Two Way Interactions (optional)
cp.int <- with(col.pop, interaction(treat, soil))
aov.col1 <- aov(surface ~ cp.int, data=col.pop)
HSD.test(aov.col1, "cp.int", group=TRUE, console=TRUE)
tuk <- glht(aov.col1, linfct = mcp(cp.int = "Tukey"))
summary(tuk)          # standard display
tuk.cld <- cld(tuk)   # letter-based display
opar <- par(mai=c(1,1,1.5,1))
plot(tuk.cld, las=2)
par(opar)

cp.int <- with(col.pop, interaction(exp, soil))
aov.col1 <- aov(surface ~ cp.int, data=col.pop)
HSD.test(aov.col1, "cp.int", group=TRUE, console=TRUE)
tuk <- glht(aov.col1, linfct = mcp(cp.int = "Tukey"))
summary(tuk)          # standard display
tuk.cld <- cld(tuk)   # letter-based display
opar <- par(mai=c(1,1,1.5,1))
plot(tuk.cld)
par(opar)

cp.int <- with(col.pop, interaction(exp, treat))
aov.col1 <- aov(surface ~ cp.int, data=col.pop)
HSD.test(aov.col1, "cp.int", group=TRUE, console=TRUE)
tuk <- glht(aov.col1, linfct = mcp(cp.int = "Tukey"))
summary(tuk)          # standard display
tuk.cld <- cld(tuk)   # letter-based display
opar <- par(mai=c(1,1,1.5,1))
plot(tuk.cld)
par(opar)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Analysis of 15N enrichment in animal tissue
## 15N enrichment in L. terrestris
bs4 <- diff(range(ew.iso$atpercent))/length(ew.iso$atpercent)
ggplot(ew.iso, aes(x=atpercent)) + 
  geom_histogram(binwidth = bs4 , fill="lightblue3", colour="black") + 
  xlab(expression(paste(.^15,"N atom%"))) + 
  ylab("frequency") +
  ggtitle(expression(paste("Histogram of", .^15,"N Enrichment in tissue of ",italic(L.~terrestris))))


# Overlap between Loam and Sand in %
100*((with(ew.iso[ew.iso$soil=="S",], max(atpercent)) - with(ew.iso[ew.iso$soil=="L",], min(atpercent)))/with(ew.iso,max(atpercent) - min(atpercent)))

maxRows <- by(ew.iso, ew.iso$soil, function(X) X[which.max(X$atpercent),])
minRows <- by(ew.iso, ew.iso$soil, function(X) X[which.min(X$atpercent),])
do.call("rbind", maxRows)
do.call("rbind", minRows)

ew.iso.plot <-  qplot(soil, atpercent, colour=treat, data=ew.iso) + 
  geom_point(size=2.5) +
  annotate("text", x=2,y=1.01+.02, label="max=1.01", vjust=0, col="red", size=4) +
  annotate("text", x=1,y=0.74+.02, label="min=0.74", vjust=0,col="red", size=4) +
  xlab("Soiltype") + 
  ylab(expression(paste(.^15,"N atom%",sep=""))) + 
  ggtitle(expression(paste("Label in ",italic("L. terrestris"))))+
  scale_x_discrete(breaks=c("L","S"), labels=c("Loam", "Sand"))+
  labs(colour="Treatment")+
  scale_colour_discrete(labels=c("L.terrestris", "L. terrestris + F. candida")) + 
  theme(panel.background= element_rect(fill="white"),
        panel.grid.major = element_line(colour="grey"),
        panel.grid.minor = element_line(colour="grey", linetype="dashed", size=0.2),
        panel.border= element_rect(colour="black", fill=NA, size=1),
        axis.text = element_text(colour="black", size=rel(1.3)),
        axis.text.x = element_text(face="italic"),
        axis.title=element_text(size=rel(1.4)),
        legend.text=element_text(face="italic"))
ew.iso.plot

# Analysis of 15N enrichment in L. terrestris

M.ew.iso <- lm(atpercent~treat*soil, data=ew.iso);anova(M.ew.iso);summary(M.ew.iso)
M.ew.iso <- lm(atpercent~soil, data=ew.iso); anova(M.ew.iso);summary(M.ew.iso)
lm.influence(M.ew.iso)
influence.measures(M.ew.iso)$is.inf

par(mfrow=c(2,2),
    mar=c(3,3,3,3))
plot(M.ew.iso)

shapiro.test(ew.iso$atpercent)
fligner.test(atpercent~soil, data=ew.iso)

t.test(atpercent~soil, ew.iso)
t.test(atpercent~soil,ew.iso, alternative="greater")

with(ew.iso, list( tapply(atpercent, soil, mean),
                   tapply(atpercent, soil, se)))
L. terrestris was significantly more enriched in loamy soil, compared to sandy soil 15N in tissue was with 1.04 +/- 0.08 atom% greater than  0.75 +/- 0.08 in sandy soil.

## 15N enrichment in F.candida
#Histograms
par(mfrow=c(2,1))
hist(col.iso$atpercent)
hist(col.iso$atpercent, main="breaks=12", breaks=12) # Too narrow --> breaks in Histogram

list(range(col.iso$atpercent),
     range(col.iso$atpercent)[2] -
       range(col.iso$atpercent)[1])

#There's a differnce between the Minima and Maxima of 0.896 atom percent
# what differnce in atpercent is meaningful?

with(col.iso, list(median(atpercent),
                   mean(atpercent))) # median and mean are the same: 6.4 atom percent

icept = median(col.iso$atpercent)

col.iso.plot <- 
  qplot(soil, atpercent, colour=treat, data=col.iso) + 
  geom_point(size=2.5) + 
  geom_hline(yintercept=icept, linetype="dashed", colour="red") +
  annotate("text", label="Midline Max(Loam) - Min(Sand) y=6.4", x=1.5, y=6.44, size=4,colour="red") +
  ggtitle(expression(paste("Label in ",italic("F. candida"))))+  
  xlab("soiltexture") + 
  ylab(expression(paste(.^15,"N atom%",sep=""))) + 
  scale_x_discrete(breaks=c("L","S"), labels=c("Loam", "Sand"))+
  scale_colour_discrete(labels=c("F.candida", "F. candida \n + L. terrestris")) + 
  labs(colour="Treatment")+
  theme(panel.background= element_rect(fill="white"),
        panel.grid.major = element_line(colour="grey"),
        panel.grid.minor = element_line(colour="grey", linetype="dashed", size=0.2),
        panel.border= element_rect(colour="black", fill=NA, size=1),
        axis.text = element_text(colour="black", size=rel(1.3)),
        axis.text.x = element_text(face="italic"),
        axis.title=element_text(size=rel(1.4)),
        legend.text=element_text(face="italic"))
col.iso.plot
````

```{r Analysis4}
M.col.iso1 <- with(col.iso, aov(atpercent~treat*soil)); summary(M.col.iso1)
M.col.iso2 <- with(col.iso, aov(atpercent~soil)); summary(M.col.iso2)
# What other nutritional source could F. candida have had in Loam?

anova(M.col.iso1, M.col.iso2)
AIC(M.col.iso1, M.col.iso2)

M.col.iso3 <- lm(atpercent~soil, data=col.iso)
lm.influence(M.col.iso3)
influence.measures(M.col.iso3)$is.inf

t.test(atpercent~soil, col.iso, alternative="less")

with(col.iso, list( tapply(atpercent, soil, mean),
                    tapply(atpercent, soil, se)))

mean(col.iso$atpercent)
se(col.iso$atpercent)


organism <- factor(rep(c("italic(L.~terrestris)", "italic(F.~candida)"), each=12), levels=c("italic(L.~terrestris)", "italic(F.~candida)"))
ew.col.iso <- rbind(ew.iso, col.iso)
ew.col.iso$organism <- organism
ew.col.iso$soil <- revalue(ew.col.iso$soil, c(L="Loam", S="Sand"))

#long cut way to find number of facets
len <- length(levels(ew.col.iso$organism))# *  length(levels(mtcars$am))
vars <- data.frame(organism=levels(ew.col.iso$organism))
y = with(ew.col.iso, tapply(atpercent,organism,max)-0.025)
dat <- data.frame(x = rep(0.5, len), y = y, vars, labs=LETTERS[1:len])

ggplot(ew.col.iso, aes(x=soil, y=atpercent), family="Arial") + 
  stat_boxplot(geom="errorbar", coef=1.5, lwd=0.2) + 
  geom_boxplot(fill="grey", lwd=0.2, outlier.size=0.4, outlier.shape=21) +  
  facet_grid(organism ~ . , scales="free", space="free",labeller=label_parsed) + 
  ylab("atom %") + 
  xlab("soil texture") + 
  geom_text(aes(x, y, label=labs, group=NULL),data=dat) +
  # ggtitle(expression(paste(.^15,"N enrichment in animal tissue"))) +
  mytheme



#setwd("D:/Quentin Schorpp/Schreibtisch/2 - Laborversuche/Mikrokosmen/Manuskript/Entwurf")
#ggsave("Figure1.pdf", width=9, height=7, units="cm", useDingbats=FALSE) 
#embed_fonts("Figure1.pdf", outfile="Figure1.pdf")
#ggsave("Figure1.svg", width=9, height=7, units="cm") 
```
