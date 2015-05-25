


# Analysis
```{r}
# data = n2o.cum

n2o.index$Lt <- revalue(n2o.index$treat, c(Lt=1, int=1, Fc=0, C=0))
n2o.index$Fc <- revalue(n2o.index$treat, c(Lt=0, int=1, Fc=1, C=0))

n2o.cum <- as.data.frame(cbind(n2o.index[c(1:56),][-c(8,22),], n2o.cum))
colnames(n2o.cum)[21] <- "nc"
colnames(n2o.cum)[4] <- "exp"

n2o.cum <- droplevels(n2o.cum)

### Experiment 1

# multimodel averaging
x <- c("Lt","Fc","exp","soil")
nc.0 <- lm(nc~Lt*Fc*soil, n2o.cum[!(n2o.cum$exp=="exp2"),])
#t(combn(x,3))
nc.1 <- lm(nc~Lt*Fc, n2o.cum[!(n2o.cum$exp=="exp2"),])
nc.2 <- lm(nc~Lt*soil, n2o.cum[!(n2o.cum$exp=="exp2"),])
nc.3 <- lm(nc~Fc*soil, n2o.cum[!(n2o.cum$exp=="exp2"),])
#t(combn(x,2))
nc.4 <- lm(nc~Lt, n2o.cum[!(n2o.cum$exp=="exp2"),])
nc.5 <- lm(nc~Fc, n2o.cum[!(n2o.cum$exp=="exp2"),])
nc.6 <- lm(nc~soil, n2o.cum[!(n2o.cum$exp=="exp2"),])

nc.7 <- lm(nc~1, n2o.cum[!(n2o.cum$exp=="exp2"),])

AICctab(nc.0,  nc.1,  nc.2,  nc.3,  nc.4,  nc.5,  nc.6,  nc.7)

#as.character(rep("nc",c(1:15)))
#x <- factor(paste("nc",".",c(0:15),",",sep=""), levels=)
#print(x)
#x <- factor(paste("nc",".",c(0:15)," <- ",sep=""), levels=)
#as.data.frame(x)

summary.aov(nc.2)
par(mfrow=c(2,2))
plot(nc.2)
par(mfrow=c(1,1))
boxcox(nc.2)
#locator(1)

nc.2 <- update(nc.2,I(1/nc^(1/5))~.)
summary.aov(nc.2)
par(mfrow=c(2,2))
plot(nc.2)

nc.aov2 <- aov(nc~Lt*soil, n2o.cum[!(n2o.cum$exp=="exp2"),])
nc.aov2 <- update(nc.aov2,I(1/nc^(1/5))~.)
TukeyHSD(nc.aov2)
nc.int <- with(n2o.cum[!(n2o.cum$exp=="exp2"),], interaction(soil,Lt)) # Double interaction factor
nc.aov3 <- update(nc.aov2,.~nc.int)
HSD.test(nc.aov3, "nc.int", group=TRUE, console=TRUE)

txt.nc = expression(paste("Nitrous Oxide ", N[2],"O [?g ",N[2],"O-N ",kg^-1,"]"))
ggplot(n2o.cum[!(n2o.cum$exp=="exp2"),], aes(x=Lt, y=nc), family="Arial") + 
  stat_boxplot(geom="errorbar", coef=1.5, lwd=0.2) + 
  geom_boxplot(fill="grey", lwd=0.2, outlier.size=0.4, outlier.shape=21) +  
  facet_grid( .~ soil , scales="free", space="free",labeller=label_parsed) + 
  ylab(txt.nc) + 
  xlab("L. terrestris") + 
  scale_x_discrete(labels=c("absent", "present")) +
  ggtitle("Experiment 1") +
  mytheme


nc.int <- with(n2o.cum[!(n2o.cum$exp=="exp2"),], interaction(soil,Lt)) # Double interaction factor
nc.aov3 <- update(nc.aov2,.~nc.int)
HSD.test(nc.aov3, "nc.int", group=TRUE, console=TRUE)

with(n2o.cum[!(n2o.cum$exp=="exp1"),], list(round(tapply(nc, interaction(treat,soil), mean),0),
                                            round(tapply(nc, interaction(treat,soil), se),0)))


with(n2o.cum[!(n2o.cum$exp=="exp2"),], list(round(tapply(nc, interaction(Lt,soil), mean),0),
                                            round(tapply(nc, interaction(Lt,soil), se),0)))


### Experiment 1

# multimodel averaging
x <- c("Lt","Fc","exp","soil")
nc.0 <- lm(nc~Lt*Fc*soil, n2o.cum[!(n2o.cum$exp=="exp1"),])
#t(combn(x,3))
nc.1 <- lm(nc~Lt*Fc, n2o.cum[!(n2o.cum$exp=="exp1"),])
nc.2 <- lm(nc~Lt*soil, n2o.cum[!(n2o.cum$exp=="exp1"),])
nc.3 <- lm(nc~Fc*soil, n2o.cum[!(n2o.cum$exp=="exp1"),])
#t(combn(x,2))
nc.4 <- lm(nc~Lt, n2o.cum[!(n2o.cum$exp=="exp1"),])
nc.5 <- lm(nc~Fc, n2o.cum[!(n2o.cum$exp=="exp1"),])
nc.6 <- lm(nc~soil, n2o.cum[!(n2o.cum$exp=="exp1"),])

nc.7 <- lm(nc~1, n2o.cum[!(n2o.cum$exp=="exp1"),])

AICctab(nc.0,  nc.1,  nc.2,  nc.3,  nc.4,  nc.5,  nc.6,  nc.7)

#as.character(rep("nc",c(1:15)))
#x <- factor(paste("nc",".",c(0:15),",",sep=""), levels=)
#print(x)
#x <- factor(paste("nc",".",c(0:15)," <- ",sep=""), levels=)
#as.data.frame(x)

summary.aov(nc.2)
summary.aov(nc.0)
par(mfrow=c(2,2))
plot(nc.2)
par(mfrow=c(1,1))
boxcox(nc.2)
# locator(1)
fligner.test(nc~interaction(Lt,soil), n2o.cum[!(n2o.cum$exp=="exp1"),])


nc.aov2 <- aov(nc~Lt*soil, n2o.cum[!(n2o.cum$exp=="exp1"),])
TukeyHSD(nc.aov2)
nc.int <- with(n2o.cum[!(n2o.cum$exp=="exp1"),], interaction(soil,Lt)) # Double interaction factor
nc.aov3 <- update(nc.aov2,.~nc.int)
HSD.test(nc.aov3, "nc.int", group=TRUE, console=TRUE)

txt.nc = expression(paste("Nitrous Oxide ", N[2],"O [?g ",N[2],"O-N ",kg^-1,"]"))
ggplot(n2o.cum[!(n2o.cum$exp=="exp1"),], aes(x=Lt, y=nc), family="Arial") + 
  stat_boxplot(geom="errorbar", coef=1.5, lwd=0.2) + 
  geom_boxplot(fill="grey", lwd=0.2, outlier.size=0.4, outlier.shape=21) +  
  facet_grid( .~ soil , scales="free", space="free",labeller=label_parsed) + 
  ylab(txt.nc) + 
  xlab("L. terrestris") + 
  scale_x_discrete(labels=c("absent", "present")) +
  ggtitle("Experiment 2") +
  mytheme


nc.int <- with(n2o.cum[!(n2o.cum$exp=="exp1"),], interaction(soil,Lt,Fc)) # Double interaction factor
nc.aov3 <- update(nc.aov2,.~nc.int)
HSD.test(nc.aov3, "nc.int", group=TRUE, console=TRUE)

with(n2o.cum[!(n2o.cum$exp=="exp1"),], list(round(tapply(nc, interaction(Lt,soil), mean),0),
                                            round(tapply(nc, interaction(Lt,soil), se),0)))


rm(n2o.bod, n2o.prod)


# Calculate mean max values for peaks
n2o.max <- melt(n2o, id.vars=1)
which.max(n2o.max$value)
n2o.max[which(n2o.max$value>15),]



n2o.PR <- n2o
n2o.PR$max.iniloam <- apply(n2o[,-c(1,6:16)],1,max)
n2o.PR$max.ewloam <- apply(n2o[,-c(1:5)],1,max)
n2o.PR$max.sand <- apply(n2o[,-1],1,max)
n2o.PR  <- cbind(n2o.index[1:56,],n2o.PR[,17:19])

with(n2o.PR, tapply(max.sand, list(experiment, soil), mean))
with(n2o.PR, tapply(max.sand, list(experiment, soil), se))