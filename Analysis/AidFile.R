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


cp.int <- with(col.pop, interaction(exp, soil))
aov.col1 <- aov(surface ~ cp.int, data=col.pop)
HSD.test(amod, "cp.int", group=TRUE, console=TRUE)

tuk <- glht(aov.col1, linfct = mcp(cp.int = "Tukey"))
summary(tuk)          # standard display
tuk.cld <- cld(tuk)   # letter-based display
opar <- par(mai=c(1,1,1.5,1))
plot(tuk.cld)
par(opar)


library(foreign)
tx <- with(col.pop, interaction(exp, treat))
amod <- aov(surface ~ tx, data=col.pop)
library(agricolae)
HSD.test(amod, "tx", group=TRUE, console=TRUE)


library(multcomp)
tuk <- glht(amod, linfct = mcp(tx = "Tukey"))
summary(tuk)          # standard display
tuk.cld <- cld(tuk)   # letter-based display
opar <- par(mai=c(1,1,1.5,1))
plot(tuk.cld)
par(opar)

library(foreign)
tx <- with(col.pop, interaction(soil, treat))
amod <- aov(surface ~ tx, data=col.pop)
library(agricolae)
HSD.test(amod, "tx", group=TRUE, console=TRUE)

library(multcomp)
tuk <- glht(amod, linfct = mcp(tx = "Tukey"))
summary(tuk)          # standard display
tuk.cld <- cld(tuk)   # letter-based display
opar <- par(mai=c(1,1,1.5,1))
plot(tuk.cld)
par(opar)



