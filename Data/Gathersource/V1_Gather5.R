############################
# V1_Nitrogen
# Isotopomer Data Processing
# Quentin Schorpp
# 23.05.2015
###########################

isom.data.600 <- with(isom.data.all, isom.data.all[N2O.conc > 600,])
isom.data.600$treat <- factor(isom.data.600$treat, levels=c("Lt", "Fc", "Int", "Control")) 

with(isom.data.all, plot(d18o,sp))
with(isom.data.600, plot(d18o,sp))
col=soil, fill=as.factor(Sampleset)),

ggplot(isom.data.600, aes(x=d18o, y=sp, shape=treat),  labeler=label_parsed) +
  #geom_point(aes(x=d18o, y=sp),size=1.6, lwd=0.01,col="black", shape=16, isom.data.600[(isom.data.600$soil=="Loam"),]) +
  geom_point(aes(x=d18o, y=sp),size=1.4, lwd=0.01,col="white", shape=16, isom.data.600[(isom.data.600$soil=="Loam"),]) +
  geom_point(size=2.0, lwd=0.2, col="black", show_guide=TRUE) +
  geom_point(size=1.2, lwd=0.2, aes(colour=as.factor(interaction(Sampleset,soil)))) +
  geom_rect( mapping=aes(xmin=10, xmax=20, ymin=-10, ymax=0), color="grey", alpha=0.01) +
  geom_rect( mapping=aes(xmin=40, xmax=50, ymin=33, ymax=36), color="grey", alpha=0.01) +
  geom_rect( mapping=aes(xmin=30, xmax=40, ymin=34, ymax=37), color="grey", alpha=0.01) +
  geom_abline(intercept = c(-12.5,0), slope = 0.5, alpha=0.1)+
  geom_abline(intercept = c(-20,0), slope = 0.5, alpha=0.1)+
  geom_abline(intercept = c(-5,0), slope = 0.5, alpha=0.1)+
  geom_abline(intercept = c(-8,0), slope = 0.2, alpha=0.1)+
  geom_abline(intercept = c(-14,0), slope = 0.2, alpha=0.1)+
  geom_abline(intercept = c(-2,0), slope = 0.2, alpha=0.1)+
  annotate("text",x=35, y=35.5, label="fungal nitrification", size=1.6, colour="black" ,face="italic") +
  annotate("text",x=45, y=34.5, label="nitrification", size=1.6, colour="black" ,face="italic") +
  annotate("text", x=15, y=-5, label="denitrification", size=1.6, col="black") +
  annotate("text", x=45, y=20, label="y=0.5x+c", size=1.6, col="black") +
  annotate("text", x=40, y=-8, label="y=0.2x+c", size=1.6, col="black") +
  ylab("Site preference [\211]") +
  xlab(expression(paste({delta}^18*O[soil]," [\211]"))) +
  scale_y_continuous(limits=c(-10,40),breaks=seq(-10,40,10)) +
  #scale_shape_manual(values=c(21,23,24,22)) +
  scale_shape_manual(labels=c("L. terrestris","F. candida","Interaction","Control"),values=c(17,16,15,18)) +
  scale_colour_manual("soil texture",labels=c("Loam (day 3)", "Loam (day 21)", "Sand (day 21)"),values=c("grey", "white", "black")) +
  mytheme + theme(legend.position="none", 
                  axis.title.y = element_text(size=8, face="plain", family="Times New Roman"))

+ theme(legend.position=c(0.1,0.75), 
        legend.text=element_text(size=6,face="plain", family="Times New Roman"),
        text=element_text(lineheight=0.4),
        legend.margin=unit(0.2, "cm")) +
  guides(shape = guide_legend(label.theme =element_text(family="Times New Roman", face = "italic", size=6, angle=0)))


#ggsave("Figure 5.pdf", width=9, height=7, units="cm", useDingbats=FALSE)

locator(1)


with(isom.data.600, list(round(tapply(sp, list(treat, soil, Sampleset), mean),0),
                         round(tapply(sp, list(treat, soil, Sampleset), se),0)))



isom.lm1 <- lm(sp~interaction(Sampleset,soil),isom.data.600)
summary(isom.lm1)



isom.test  <-  cbind(isom.data.600$sp, isom.data.600$d18o)
isom.mv1 <- manova(isom.test~interaction(Sampleset,soil),isom.data.600)
summary(isom.mv1, test="Pillai")
summary.aov(isom.mv1)

isom.aov1 <- aov(sp~interaction(Sampleset,soil),isom.data.600)
summary(isom.aov1)
TukeyHSD(isom.aov1)

isom.aov2 <- aov(d18o~interaction(Sampleset,soil),isom.data.600)
summary(isom.aov2)
TukeyHSD(isom.aov2)


with(isom.data.600, by(sp, interaction(soil,Sampleset),summary))
with(isom.data.600, by(d18o, interaction(soil,Sampleset),summary))

# Reduction of n2o to n2, Rayleigh Equation
# Rayleigh <- SP = SP0 + eSP*ln(f) // Sp = Sp0 + eSP*(1-f)
# nach F aufgel?st:
# f = e^((SP-SP0)/eSP) // f = 1-(SP-SP0)/eSP

isom.loam.set2 <- isom.data.600[with(isom.data.600,(Sampleset==2 & soil=="Loam")),]
sp <- mean(isom.loam.set2$sp)+1000
sp0 <- -5+1000
fracfac <- 1000/-5

f  <-  (sp/sp0)^fracfac
round(f,2)
