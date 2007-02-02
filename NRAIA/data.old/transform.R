nms <- c("BOD2","Chloride","Coal","Ethyl","Isom","Leaves","Lipo",
         "Lubricant","Nitren","Nitrite","Oilshale","O.xylene",  
         "PCB","Pinene","Pinene2","Rumford","Sacch2","Saccharin", 
         "sPMMA","Sulfi","Tetra")
base <- "/home/bates/src/Rlibs/NRAIA"
data.old <- file.path(base, "data.old")
for (nm in nms) source(file.path(data.old, paste(nm, ".R", sep='')))
for (nm in nms) save(list=nm, file = file.path(data, paste(nm, ".rda", sep='')))
