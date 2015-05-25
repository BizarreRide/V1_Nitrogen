############################
# V1_Nitrogen
# N2O Data Processing
# Quentin Schorpp
# 23.05.2015
###########################

# N2O [ppb] data formatting ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# for plots use n2o.cum.mean ans n2o.cum.se
# for analysis use n2o.cum

# subset raw data 
n2o.index <- n2o.ppb[1:56,1:4]
n2o.index$treat <- factor(n2o.index$treat, levels=c("Lt", "Fc", "int", "C")) # reorder factor levels
n2o.index$code <- with(n2o.index, factor(interaction(treat,soil,experiment), levels= c("Lt.Loam.exp1", "int.Loam.exp1", "Fc.Loam.exp1", "C.Loam.exp1","Lt.Sand.exp1","int.Sand.exp1","Fc.Sand.exp1","C.Sand.exp1", "Lt.Loam.exp2","int.Loam.exp2","Fc.Loam.exp2","C.Loam.exp2","Lt.Sand.exp2","int.Sand.exp2","Fc.Sand.exp2","C.Sand.exp2")))
n2o.blanks <- n2o.ppb[57:59,5:19]
n2o.ppb <- n2o.ppb[1:56,5:19]

# mean for blanks 
# mean background concentration
n2o.blanks = data.frame(t(apply(n2o.blanks,2,mean)))

# subtract blanks from MC measurements
for (i in 1:15) {
  for (j in 1:56) {
    n2o.ppb[j,i] = n2o.ppb[j,i]-n2o.blanks[,i]
  }
}

# set negative values to zero
n2o.ppb[n2o.ppb<0] = 0
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Calculate N2O Production rate #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# multiplication with gasflow
n2o.prod=n2o.ppb*gasflow[,5:19]*10^-9*60 #[ml/h]  ### ppb beachten, zeit umrechnen
# umrechnen von ml auf µg
n2o.prod=n2o.prod*44/22.4*1000 #[µg/h]
# auf das Bodengewicht beziehen;  soil dry weight: 1080g
n2o.bod=n2o.prod/1.080  #[µg/(kg h)
# calculate N2O-N, Anteil an molarer Masse 28/44
n2o=n2o.bod*28/44  #[µg/(kg h)]
n2o <- cbind(n2o.index$code, n2o)

# MK8 und MK22 loeschen weil da Lt tot
n2o.red=n2o[c(-8,-22),]
n2o.red <- droplevels(n2o.red)

## MW nach Treatment bilden Production Rate ####
Kopf=c("Tag 1","Tag 3","Tag 5","Tag 7","Tag 11","Tag 15","Tag 18","Tag 21","Tag 25","Tag 28","Tag 32","Tag 35","Tag 39","Tag 42","Tag 48")

n2o.mean=matrix(rep(0,240),16)
for (j in 2:16) {
  a=tapply(n2o.red[,j],n2o.red[,1],mean)
  n2o.mean[,(j-1)]=a
}
 
colnames(n2o.mean)<-dimnames(n2o[,2:16])[[2]]
n2o.mean <- data.frame(n2o.mean)
n2o.mean <- cbind(unique(n2o.index[1:56,-c(1,5)][c(-8,-22),]),n2o.mean)

## Standardfehler Production rate ####
n2o.se=matrix(rep(0,240),16)
for (j in 2:16) {
  b=tapply(n2o.red[,j],n2o.red[,1],se)
  n2o.se[,(j-1)]=b  
}

colnames(n2o.se)<-dimnames(n2o[,2:16])[[2]]
n2o.se <- data.frame(n2o.se)
n2o.se <- cbind(unique(n2o.index[1:56,-c(1,5)][c(-8,-22),]),n2o.se)

## Table with means and standard error N2O ####
n2o.mean2=round(n2o.mean[,4:18],2)
n2o.se2=round(n2o.se[,4:18],2)

n2o.mean.se=matrix(rep(0,240),16)
for (i in 1:16) {
  for (j in 1:15) {
    a=paste(n2o.mean2[i,j],"?",n2o.se2[i,j])
    n2o.mean.se[i,j]=a    
  }
}

colnames(n2o.mean.se)<-dimnames(n2o[,2:16])[[2]]
n2o.mean.se <- data.frame(n2o.mean.se)
n2o.mean.se <- cbind(unique(n2o.index[1:56,-c(1,5)][c(-8,-22),]),n2o.mean.se)
rm(n2o.mean2, n2o.se2)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Calculation of cumulative N2O ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

code <- n2o.red[,1]
n2o.red <- n2o.red[,-1]
days=c(1,3,5,7,11,15,18,21,25,28,32,35,39,42,48)
n2o.rz=rbind(days,n2o.red)

## N2O between two sampling dates ####
n2o.btw=matrix(rep(0,756),54)
for (i in 2:55){
  for (j in 1:14) {
    Zeit=(n2o.rz[1,(j+1)]-n2o.rz[1,j])*24  #*24da Zeit d in std umgerechnet werden muss
    a=Zeit*n2o.rz[i,j]+0.5*Zeit*(n2o.rz[i,(j+1)]-n2o.rz[i,j])
    n2o.btw[(i-1),j]=a
  }
}

## cumulative N2O for all samples ####
n2o.cum=matrix(rep(0,756),54)
for (i in 1:54) {
  for (j in 1:14) {
    a=sum(n2o.btw[i,1:j])
    n2o.cum[i,j]=a    
  }  
}

colnames(n2o.cum)<-days[-1]
cbind(n2o.index[-c(8,22,57:59),],n2o.cum)

## Means for treatments ####
n2o.btw.mean=matrix(rep(0,224),16) 
for (j in 1:14) {
  a=tapply(n2o.btw[,j],code,mean)
  n2o.btw.mean[,j]=a  
}

## Standard error ####
n2o.btw.se=matrix(rep(0,224),16) 
for (j in 1:14) {
  a=tapply(n2o.btw[,j],code,se)
  n2o.btw.se[,j]=a
}

## summed means ####
n2o.cum.mean=matrix(rep(0,224),16)   
for (i in 1:16) {
  for (j in 1:14) {
    a=sum(n2o.btw.mean[i,1:j])
    n2o.cum.mean[i,j]=a    
  }
}

colnames(n2o.cum.mean) <- days[-1]
n2o.cum.mean <- cbind(unique(n2o.index[1:56,-c(1,5)][c(-8,-22),]),n2o.cum.mean)

## Kum- n2o: Std.Fehler durch Fehlerfortpflanzung ####
n2o.cum.se=matrix(rep(0,224),16)   
for (i in 1:16) {
  for (j in 1:14){
    a=sqrt(sum((n2o.btw.se[i,1:j])^2))
    n2o.cum.se[i,j]=a
  }
}

colnames(n2o.cum.se) <- days[-1]
n2o.cum.se <- cbind(unique(n2o.index[1:56,-c(1,5)][c(-8,-22),]),n2o.cum.se)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(code, days, Kopf, Zeit, a, b, n2o.rz)
