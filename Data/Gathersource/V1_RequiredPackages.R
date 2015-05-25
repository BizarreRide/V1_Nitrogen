############################
# V1_Nitrogen
# Load required packages
# Quentin Schorpp
# 23.05.2015
###########################

library(ggplot2)
library(reshape2)
library(grid)
library(plyr)
library(MethComp)
library(agricolae)
library(foreign)
library(multcomp)
library(gridExtra)
library(extrafont)
library(effects)
library(MASS)
library(splitstackshape)
library(bbmle)



#### Additional Functions:

# Function for Bland Altman Plot ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
baplot <- function(m1, m2, ...) {
  # m1 and m2 are the measurements
  means <- (m1 + m2) / 2
  diffs <- m1 - m2
  mdiff <- mean(diffs)
  sddiff <- sd(diffs)
  # Compute the figure limits
  ylimh <- mdiff + 3 * sddiff
  yliml <- mdiff - 3 * sddiff
  # Plot data
  plot(diffs ~ means, xlab = "Average values",
       ylab = "Differences", ylim = c(yliml, ylimh), ...)
  abline(h = mdiff) # Center line
  # Standard deviations lines
  abline(h = mdiff + 1.96 * sddiff, lty = 2)
  abline(h = mdiff - 1.96 * sddiff, lty = 2)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Function for layout ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Function for Standard error ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
se <- function(x) sqrt(var(x)/length(x))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Function to create ellipse in scatter plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#### ggplot2 - Theme ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mytheme = 
  theme_bw() + 
  theme(strip.background = element_rect(color = "grey", fill="black", size=0.1),
        strip.text.x = element_text(size=8,  colour="white", face="italic"),
        strip.text.y = element_text(size=8,  colour="white", face="italic"),
        axis.text.x = element_text(size=7),
        axis.title.x = element_text(size=8,face="bold", family="Times New Roman"),
        axis.text.y = element_text(size=7),
        axis.title.y = element_text(size=8, family="Times New Roman"),
        axis.line = element_line(size=0.25),
        axis.ticks = element_line(size=0.25),
        plot.title = element_text(size=11,face="bold", family="Times New Roman"),
        panel.margin = unit(0, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black", size=0.2, fill=NA),
        legend.key=element_blank(),
        legend.background=element_blank(),
        legend.text=element_text(size=8,face="italic", family="Times New Roman"),
        legend.title=element_text(size=8))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

