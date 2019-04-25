## 
## Florida Coastal Everglades LTER
## Biogeochemistry Working Group
## Grab Sample water quality
##
## Code was compiled by Paul Julian
## contact infor: pjulian@ufl.edu

#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

#Libraries
library(AnalystHelper);#devtools::install_github("SwampThingPaul/AnalystHelper")
library(plyr)
library(reshape)
library(zoo);

#Custom Functions

#Paths
setwd("D:/_GitHub/FCE-BGChemWG/Trends")

paths=paste0(getwd(),c("/Exports/","/Plots/"))
#Folder.Maker(paths);#One and done. Creates folders in working directory.
export.path=paths[1]
plot.path=paths[2]

# custom functions
read.lter=function(PASTA,DOI){
  
}


