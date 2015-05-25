############################
# V1_Nitrogen
# 15N N2O+N2 Data Processing
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

# Gas15 Data Processing
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gas15$date  <- as.POSIXlt(gas15$date)
gas15$treat <- factor(gas15$treat, levels = c("Lt", "Fc", "LF", "C")) #reorder levels
gas15$total.Ngas <- gas15$nitox.gc+gas15$dinitrogen


# Cumulative Gas15 (Total, N2, N2O-litter, N2O-soil) ####
# Cumulative gas between two dates
cumulative = function (x) 168*(x[i:(length(x)-j)]+0.5*(x[(i+j):length(x)]-x[i:(length(x)-j)]))
i=1
j=26

gas15.cum = matrix(0,78,4)
gas15.cum[,1] = cumulative(gas15$total)
gas15.cum[,2] = cumulative(gas15$dinitrogen)
gas15.cum[,3] = cumulative(gas15$nitox)
gas15.cum[,4] = cumulative(gas15$nitox.gc)
colnames(gas15.cum) = colnames(gas15[,c(22:24,14)])

# save data on groth between two dates
gas15.cum1 = gas15.cum

# New data to sum up
gas15.cum2 = gas15.cum

# sum up 
for (i in 27:78)  {
  for (j in 1:4) {
    gas15.cum2[i,j] =  gas15.cum2[i,j] + gas15.cum2[i-26,j]
  }
}
gas15.cum2 = cbind(gas15[27:104,1:13], gas15.cum2)

# Dataset for Analysis
gas15.cum = cbind(gas15[1:26,1:3],gas15[1:26,8:13], gas15.cum[53:78,])
gas15.cum$total.Ngas <- gas15.cum$dinitrogen + gas15.cum$nitox.gc
gas15.cum$int.SLt <- with(gas15.cum, interaction(soil,Lt)) # Double interaction factor for TukeyHSD Tests
gas15.cum$nitox.soil <- gas15.cum$nitox.gc - gas15.cum$nitox # Double interaction factor for TukeyHSD Tests

# Averages instead of cumulated fractions:
gas15.cum$nitox.pool <- as.numeric(with(gas15, tapply(nitox.pool,MK,mean)))
gas15.cum$product.ratio <- as.numeric(with(gas15, tapply(product.ratio,MK,mean)))

# dataframe for mean  over time(day)
gas15.cum2.mean <- aggregate(as.matrix(gas15.cum2[,14:17]) ~ day + treat + soil, gas15.cum2, mean, na.rm=TRUE)

# Fehler durch Fehlerfortpflanzung:
# Formel: function (x) sqrt(sum((Kum_StdF[i,1:j])^2
gas15.cum2.se <- aggregate(as.matrix(gas15.cum2[,14:17]) ~ day + treat + soil, gas15.cum2, se)
gas15.cum2.se <- melt(gas15.cum2.se, id.vars= c(1:3))
gas15.cum2.se <- dcast(gas15.cum2.se, variable + soil + treat ~ day) # day in columns
colnames(gas15.cum2.se) [1] = "gas"

for (i in 1:24)  {
  for (j in 4:6) {
    a=sqrt(sum((gas15.cum2.se[i,4:j])^2))
    gas15.cum2.se[i,j]=a
  }
}

gas15.cum2.se <- melt(cbind(gas15.cum2.se[,1:3], gas15.cum2.se), id.vars=c(1:3))
gas15.cum2.se <- dcast(gas15.cum2.se, soil+treat+variable~gas)

gas15.cum2.meanse <- cbind(gas15.cum2.mean, gas15.cum2.se[,4:7])
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




