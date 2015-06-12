############################
# V1_Nitrogen
# 15N Figures
# Quentin Schorpp
# 23.05.2015
###########################



## Variable Explanation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ID: continuous ID
# MK: MikroKosmos, !!!reduced to 3 replicates!!!
# MK_add10: Added 10 to the number of MK, .... sometimes useful
# termin: Nr of measuring date
# day: days after start of the experiment
# week: week of the experiment
# date: date of mesurement
# treat: factor treatment
# Lt: absence/presence of L. terrestris
# Fc: absence/presence of F. candida
# lit_pow: absence/presence of powdered litter
# lit_chop: absence/presence of chopped litter
# soil: factor soil type
# nitox.gc: total N2O of the sample measured with GC [unitµµ]
# nitox.pool.gcms: µµ
# nitox.atperc: 15N atom percent measured in N2O
# total.ppm: litter derived N2+N2O [ppm]
# dinitrogen.ppm: litter derived N2 [ppm]
# nitox.ppm: litter derived N2 [ppm]
# product ratio: N2o/(N2+N2O) product ratio of litter derived N-fluxes
# nitox.pool: fraction of N2O derived from the 15N pool
# total: litter derived N2 + N2O [unitµµ]
# dinitrogen: litter derived N2 [unitµµ]
# nitox: litter derived N2O [unitµµ]
# total.NGas: Dinitrogen + Nitrous Oxide [litter] + Nitrous Oxide [soil]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Data Exploration ####
## Multi-Panel Lineplot with mean and standard error #####
## for Total-N, product ratio and N2O ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#### Data Preparation ####

# df for mean over time(day)
gas15.mean <- aggregate(as.matrix(gas15[,14:25]) ~ day + treat + soil, gas15, mean, na.rm=TRUE)

# df for standard error over time(day)
gas15.se <- aggregate(as.matrix(gas15[,14:25]) ~ day + treat + soil, gas15, se)

gas15.mean.melt  <- melt(gas15.mean, id.vars= c(1:3), value.name="mean") 
gas15.se.melt    <- melt(gas15.se, id.vars = c(1:3),value.name="se") 
gas15.fig3       <- cbind(gas15.mean.melt, se=gas15.se.melt[,5])

# Subset data frame
omit <- c(levels(gas15.fig3$variable)[c(1,2,3,4,5,6,8,10,12)]) # Include only total, product ratio and liter N2O, (cut "1" to include soil derived N2O)
gas15.fig3.subset <- subset(gas15.fig3, !variable %in% omit)

# Reorder Factor levels
gas15.fig3.subset$variable  <- factor(gas15.fig3.subset$variable, levels = c("total", "product.ratio", "nitox", "dinitrogen", "nitox.pool", "nitox.gc"))

# Change notation for pretty plots:

# For facet.grid:
levels(gas15.fig3$variable)[levels(gas15.fig3$variable)=="total"] <- "paste(N[2],O)+N[2]"
levels(gas15.fig3$variable)[levels(gas15.fig3$variable)=="product.ratio"] <- "Product~Ratio"
levels(gas15.fig3$variable)[levels(gas15.fig3$variable)=="nitox"] <- "paste(N[2],O)"

# ylab for Production Rates: nitox, product ratio, total
txt.PR = expression(paste(paste(µg~paste(N[2],O)-N~kg^-1,h^-1),"       ",frac(paste(N[2],O),paste(N[2],O+N[2])),"        ",paste(µg~Total-N~kg^-1,h^-1)))

# color palettes
trt_pal = c("firebrick3", "dodgerblue2", "darkorchid", "dimgray") # for Treatments



### Figure 2 ####

pd=position_dodge(0.2)

Figure2 <- ggplot(gas15.fig3.subset, aes(x=day, y=mean, group=treat, shape=treat)) +
  facet_grid(variable~soil, scales="free", labeller=label_parsed) +
  #guides(colour=F, shape=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), size=0.25, width=.3, position=pd) + 
  geom_line(size=0.3, position=pd) +
  geom_point(size=2.0, position=pd, colour="black", fill="white") +
  scale_shape_manual(values=c(17,16,15,21), labels=c("L. terrestris", "F. candida","Interaction","Control")) +
  #scale_colour_manual(values=trt_pal) +
  #scale_colour_grey(start=0, end=0.6) +
  labs(shape="") +
  xlab("Days [d]") +
  ylab(txt.PR) +
  mytheme + theme(legend.position="bottom")                     
Figure2

## ggsave(Figure3, "Figure 2.pdf",  width=9, height=10.5, unit="cm", useDingbats=FALSE)   
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## For cumulative N2, N2O[litter], N2+N2O, N2O[total] ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

txt.cum.y  = expression(paste( ,"",paste(µg~N[2]-N~kg^-1),"       ",paste(µg~paste(N[2],O)-N~kg^-1),"        ",paste(µg~Total-N~kg^-1),"  ",paste(µg~total~paste(N[2],O)-N~kg^-1)))

helpframe1 = melt(gas15.cum2.meanse, id.vars=c(1:3), measure.vars=c(4:7), variable.name="mean.fac", value.name="mean")
helpframe2 = melt(gas15.cum2.meanse, id.vars=c(1:3), measure.vars=c(8:11), variable.name="se.fac", value.name="se")
gas15.cum2.meanse.plot <- cbind(helpframe1,helpframe2[,c(4,5)])
rm(helpframe1, helpframe2)

gas15.cum2.meanse.plot$treat = factor(gas15.cum2.meanse.plot$treat, levels=c("Lt","Fc","LF","C"))
gas15.cum2.meanse.plot$mean = as.numeric(gas15.cum2.meanse.plot$mean)
gas15.cum2.meanse.plot$mean.fac = revalue(gas15.cum2.meanse.plot$mean.fac, c(total="paste(N[2],O)+N[2]", nitox="paste(N[2],O)", dinitrogen="N[2]",nitox.gc="paste(total~N[2],O)"))


pd=position_dodge(0.2)

Figure3.1 <- ggplot(data=gas15.cum2.meanse.plot, aes(x=day, y=mean, group=treat, shape=treat)) +
  facet_grid(mean.fac~soil, scales="free", labeller=label_parsed) +
  #guides(colour=F, shape=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), size=0.25, width=.3, position=pd) + 
  geom_line(size=0.3, position=pd) +
  geom_point(size=1.5, position=pd, colour="black", fill="white") +
  scale_shape_manual(values=c(17,16,15,21), labels=c("L. terrestris", "F. candida","Interaction","Control")) +
  #scale_colour_manual(values=trt_pal) +
  #scale_colour_grey(start=0, end=0.6) +
  labs(shape="") +
  xlab("Days [d]") +
  ylab(txt.cum.y) +
  mytheme + theme(legend.position="bottom") 
Figure3.1

ggplot(gas15.cum, aes(x= Lt, y=total),labeller=label_parsed) +
  stat_boxplot(geom="errorbar", coef=1.5, lwd=0.2) + 
  geom_boxplot(fill="light blue", lwd=0.2, outlier.size=0.4, outlier.shape=21) +
  ylab(txt.total) +
  xlab("L.~terrestris") +
  facet_grid(.~soil, labeller=label_parsed) + 
  mytheme

ggplot(gas15.cum, aes(x= Lt, y=nitox.gc),labeller=label_parsed) +
  stat_boxplot(geom="errorbar", coef=1.5, lwd=0.2) + 
  geom_boxplot(fill="light blue", lwd=0.2, outlier.size=0.4, outlier.shape=21) +
  ylab(expression(paste("total.nitox.cum"))) +
  xlab("L.~terrestris") +
  facet_grid(.~soil, labeller=label_parsed) + 
  mytheme

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



## Ellipse plot for Fraction of litter-derived nitrous oxide ####
# vs. litter derived N2O
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### Data Preparation ####

# Mean values over MK
gas15.fig4 <- aggregate(as.matrix(gas15[,14:24]) ~ MK + treat + soil, gas15, mean, na.rm=TRUE)

gas15.fig4$Lt <- factor(rep(rep(c(" + Lt", " - Lt"," + Lt", " - Lt"), c(4,3,3,3)),2))
gas15.fig4$fSLt <- interaction(gas15.fig4$soil, gas15.fig4$Lt, sep="")
gas15.fig4 <- droplevels(gas15.fig4)   

# Color palettes
LtSoil_pal=c("gray8", "gray30", "gray48", "gray70")
LtSoil_pal2=c("darkseagreen", "burlywood", "darkolivegreen3", "darkgoldenrod1")

# Expressions
title = expression(paste("Plant derived ",N[2],"O depending on soiltype and presence of L. terrestris"))
txt.pool = expression(paste("Plant derived ", N[2],"O"," [%]"))
txt.nitox = expression(paste(N[2],"O"[litter]," [µg ",N[2],"O-N ",kg^-1, h^-1,"]"))
txt.nitox.gc = expression(paste(N[2],"O"[total]," [µg ",N[2],"O-N ",kg^-1, h^-1,"]"))



### Ellipse Building ####
# Script from http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo

# For annotation
ellipse=with(gas15.fig4, aggregate(gas15.fig4[,c(14,11)],list(group=fSLt),mean))

# Data frame df_ell contains values to show ellipses. 
# It is calculated with function veganCovEllipse which is hidden in vegan package. 
# This function is applied to each level of NMDS (group) and 
# it uses also function cov.wt to calculate covariance matrix.

df_ell <- data.frame()
for(g in levels(gas15.fig4$fSLt)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(gas15.fig4[gas15.fig4$fSLt==g,],
                                                   veganCovEllipse(cov.wt(cbind(nitox,nitox.pool), wt=rep(1/length(nitox),length(nitox)))$cov,
                                                                   center=c(mean(nitox),mean(nitox.pool))))),
                                group=g))
}


### Figure 4 ####

Figure4 = ggplot(data=gas15.fig4, aes(x=nitox, y=nitox.pool)) +
  #guides(colour=F, shape=F) +
  geom_point(aes(colour=fSLt, shape=treat), size=2.0) +
  geom_path(data=df_ell, aes(x=nitox, y=nitox.pool, colour=group), size=0.5, linetype=2) +
  annotate("text",x=ellipse$nitox+0.1,y=ellipse$nitox.pool,label=ellipse$group, size=2.5) +
  scale_shape_manual(values=c(17,16,15,21)) +
  scale_colour_manual(values=LtSoil_pal) +
  scale_y_continuous(limits=c(0,1)) +
  ylab(txt.pool) +
  xlab(txt.nitox) +
  mytheme + theme(legend.position="none" )
Figure4

# ggsave(Figure4, "Figure 4.pdf", width=9, height=7, units="cm", useDingbats=FALSE)
# ggsave("Figure 4.svg", width=9, height=7, units="cm")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## vs. total N2O
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# For annotation
ellipse2=with(gas15.fig4, aggregate(gas15.fig4[,c(4,11)],list(group=fSLt),mean))

df_ell2 <- data.frame()
for(g in levels(gas15.fig4$fSLt)){
  df_ell2 <- rbind(df_ell2, cbind(as.data.frame(with(gas15.fig4[gas15.fig4$fSLt==g,],
                                                     veganCovEllipse(cov.wt(cbind(nitox.gc,nitox.pool),wt=rep(1/length(nitox.gc),length(nitox.gc)))$cov,center=c(mean(nitox.gc),mean(nitox.pool)))))
                                  ,group=g))
}


Figure5 = ggplot(data=gas15.fig4, aes(x=nitox.gc, y=nitox.pool)) +
  #guides(colour=F, shape=F) +
  geom_point(aes(colour=fSLt, shape=treat), size=1.5) +
  geom_path(data=df_ell2, aes(x=nitox.gc, y=nitox.pool, colour=group), size=0.5, linetype=2) +
  annotate("text",x=ellipse2$nitox.gc+0.1,y=ellipse2$nitox.pool,label=ellipse2$group, size=2.5) +
  scale_shape_manual(values=c(17,16,15,21)) +
  scale_colour_manual(values=LtSoil_pal) +
  scale_y_continuous(limits=c(0,1)) +
  ylab(txt.pool) +
  xlab(txt.nitox.gc) +
  mytheme +theme(legend.position="none")
Figure5
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

