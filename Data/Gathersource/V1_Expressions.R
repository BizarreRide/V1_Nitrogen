############################
# V1_Nitrogen
# V1_Expressions
# Quentin Schorpp
# 23.05.2015
###########################


# Expressions
# x and y labs for N2O Production rates: 
txt.pool = expression(paste("Plant derived ", N[2],"O"," [%]"))
txt.nitox = expression(paste(N[2],"O"[litter]," [µg ",N[2],"O-N ",kg^-1, h^-1,"]"))
txt.nitox.gc = expression(paste(N[2],"O"[total]," [µg ",N[2],"O-N ",kg^-1, h^-1,"]"))
txt.nitox.soil = expression(paste(N[2],"O"[soil]," [µg ",N[2],"O-N ",kg^-1, h^-1,"]"))
txt.total      = expression(paste("Total ", N[2],"O+",N[2]," [µg ",N[2],"(O)-N ",kg^-1, h^-1,"]"))
txt.dinitrogen = expression(paste("Dinitrogen ", N[2]," [µg ",N[2],"-N ",kg^-1, h^-1,"]"))
txt.ratio      = expression(paste("Product Ratio ", frac(paste(N[2],O),paste(N[2],O+N[2]))))

# ylab for Production Rates: nitox, product ratio, total
txt.PR = expression(paste(paste(µg~paste(N[2],O)-N~kg^-1,h^-1),"       ",frac(paste(N[2],O),paste(N[2],O+N[2])),"        ",paste(µg~Total-N~kg^-1,h^-1)))

# ylab for cumulated variables: dinitrogen, nitox, total, nitox.gc
txt.cum  = expression(paste( ,"",paste(µg~N[2]-N~kg^-1),"       ",paste(µg~paste(N[2],O)-N~kg^-1),"        ",paste(µg~Total-N~kg^-1),"  ",paste(µg~total~paste(N[2],O)-N~kg^-1)))

# title for ellipse plot of litter derived N2O fractions
title.ellipse = expression(paste("Plant derived ",N[2],"O depending on soiltype and presence of L. terrestris"))





