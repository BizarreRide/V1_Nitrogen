############################
# V1_Nitrogen
# CO2 Data Processing
# Quentin Schorpp
# 25.05.2015
###########################

# co2 [ppm] data formatting ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# for plots use co2.cum.mean ans co2.cum.se
# for analysis use co2.cum

# subset raw data 
co2.index <- co2.ppm[1:56,1:4]
co2.index$treat <- factor(co2.index$treat, levels=c("Lt", "Fc", "int", "C")) # reorder factor levels
co2.index$code <- with(co2.index, factor(interaction(treat,soil,experiment), levels= c("Lt.Loam.exp1", "int.Loam.exp1", "Fc.Loam.exp1", "C.Loam.exp1","Lt.Sand.exp1","int.Sand.exp1","Fc.Sand.exp1","C.Sand.exp1", "Lt.Loam.exp2","int.Loam.exp2","Fc.Loam.exp2","C.Loam.exp2","Lt.Sand.exp2","int.Sand.exp2","Fc.Sand.exp2","C.Sand.exp2")))
co2.blanks <- co2.ppm[57:59,5:19]
co2.ppm <- co2.ppm[1:56,5:19]

# mean for blanks 
# mean background concentration
co2.blanks = data.frame(t(apply(co2.blanks,2,mean)))

# subtract blanks from MC measurements
for (i in 1:15) {
  for (j in 1:56) {
    co2.ppm[j,i] = co2.ppm[j,i]-co2.blanks[,i]
  }
}

# set negative values to zero
co2.ppm[co2.ppm<0] = 0
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Calculate co2 Production rate #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# multiplication with gasflow
co2.prod=co2.ppm*gasflow[,5:19]*10^-6*60 #[ml/h]  ### ppm beachten, zeitumrechnen
### umrechnen von ml auf mg
co2.prod=co2.prod*44/22.4*1000 #[mg/h]
### auf das Bodengewicht beziehen;  soil dry weight: 1080g
co2.bod=co2.prod/1.080  #[mg/(kg h)
### calculate co2-N, Anteil an molarer Masse 28/44
co2=co2.bod*12/44  #[mg/(kg h)]
co2 <- cbind(co2.index$code, co2)

# MK8 und MK22 loeschen weil da Lt tot
co2.red=co2[c(-8,-22),]
co2.red <- droplevels(co2.red)

## MW nach Treatment bilden Production Rate ####

co2.mean=matrix(rep(0,240),16)
for (j in 2:16) {
  a=tapply(co2.red[,j],co2.red[,1],mean)
  co2.mean[,(j-1)]=a
}

colnames(co2.mean)<-dimnames(co2[,2:16])[[2]]
co2.mean <- data.frame(co2.mean)
co2.mean <- cbind(unique(co2.index[1:56,-c(1,5)][c(-8,-22),]),co2.mean)

## Standardfehler Production rate ####
co2.se=matrix(rep(0,240),16)
for (j in 2:16) {
  b=tapply(co2.red[,j],co2.red[,1],se)
  co2.se[,(j-1)]=b  
}

colnames(co2.se)<-dimnames(co2[,2:16])[[2]]
co2.se <- data.frame(co2.se)
co2.se <- cbind(unique(co2.index[1:56,-c(1,5)][c(-8,-22),]),co2.se)

## Table with means and standard error co2 ####
co2.mean2=round(co2.mean[,4:18],2)
co2.se2=round(co2.se[,4:18],2)

co2.mean.se=matrix(rep(0,240),16)
for (i in 1:16) {
  for (j in 1:15) {
    a=paste(co2.mean2[i,j],"m",co2.se2[i,j])
    co2.mean.se[i,j]=a    
  }
}

colnames(co2.mean.se)<-dimnames(co2[,2:16])[[2]]
co2.mean.se <- data.frame(co2.mean.se)
co2.mean.se <- cbind(unique(co2.index[1:56,-c(1,5)][c(-8,-22),]),co2.mean.se)
rm(co2.mean2, co2.se2)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Calculation of cumulative co2 ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

code <- co2.red[,1]
co2.red <- co2.red[,-1]
days=c(1,3,5,7,11,15,18,21,25,28,32,35,39,42,48)
co2.rz=rbind(days,co2.red)

## co2 between two sampling dates ####
co2.btw=matrix(rep(0,756),54)
for (i in 2:55){
  for (j in 1:14) {
    Zeit=(co2.rz[1,(j+1)]-co2.rz[1,j])*24  #*24da Zeit d in std umgerechnet werden muss
    a=Zeit*co2.rz[i,j]+0.5*Zeit*(co2.rz[i,(j+1)]-co2.rz[i,j])
    co2.btw[(i-1),j]=a
  }
}

## cumulative co2 for all samples ####
co2.cum=matrix(rep(0,756),54)
for (i in 1:54) {
  for (j in 1:14) {
    a=sum(co2.btw[i,1:j])
    co2.cum[i,j]=a    
  }  
}

colnames(co2.cum)<-days[-1]
cbind(co2.index[-c(8,22,57:59),],co2.cum)

## Means for treatments ####
co2.btw.mean=matrix(rep(0,224),16) 
for (j in 1:14) {
  a=tapply(co2.btw[,j],code,mean)
  co2.btw.mean[,j]=a  
}

## Standard error ####
co2.btw.se=matrix(rep(0,224),16) 
for (j in 1:14) {
  a=tapply(co2.btw[,j],code,se)
  co2.btw.se[,j]=a
}

## summed means ####
co2.cum.mean=matrix(rep(0,224),16)   
for (i in 1:16) {
  for (j in 1:14) {
    a=sum(co2.btw.mean[i,1:j])
    co2.cum.mean[i,j]=a    
  }
}

colnames(co2.cum.mean) <- days[-1]
co2.cum.mean <- cbind(unique(co2.index[1:56,-c(1,5)][c(-8,-22),]),co2.cum.mean)

## Kum- co2: Std.Fehler durch Fehlerfortpflanzung ####
co2.cum.se=matrix(rep(0,224),16)   
for (i in 1:16) {
  for (j in 1:14){
    a=sqrt(sum((co2.btw.se[i,1:j])^2))
    co2.cum.se[i,j]=a
  }
}

colnames(co2.cum.se) <- days[-1]
co2.cum.se <- cbind(unique(co2.index[1:56,-c(1,5)][c(-8,-22),]),co2.cum.se)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(code, days, Zeit, a, b, co2.rz)
