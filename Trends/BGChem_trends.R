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
read.lter=function(data.package,PASTA,DOI){
  prefix="http://pasta.lternet.edu/package/data/eml/"
  
  infile1=paste0(prefix,data.package,PASTA,DOI)
  dt1=read.csv(infile1,na.strings=c("-9999","-9999.00","-9999.000"))
  return(dt1)
}


## Data Import
dat.sources=data.frame(Dataset=c("LT_ND_Losada_002","LT_ND_Grahl_002","LT_ND_Rubio_002","LT_ND_Rondeau_002"),
                      PASTA=c("1075/8/","1073/12/","1080/8/","1077/3/"),
                      DOI=c("ac7159e66cbad75abb61bd1992f8d2c0","15e6ff92f875a272ba6d98db3f738026","3a84ab2009eda6b73e375aa3ad56da1a","832e5493edf5a5d79d8a60535e26a012"))

fce.losada=with(dat.sources[1,],read.lter("knb-lter-fce/",PASTA,DOI))
fce.grahl=with(dat.sources[2,],read.lter("knb-lter-fce/",PASTA,DOI))
fce.rubio=with(dat.sources[3,],read.lter("knb-lter-fce/",PASTA,DOI))
fce.rondeau=with(dat.sources[4,],read.lter("knb-lter-fce/",PASTA,DOI))

names(fce.grahl)
names(fce.losada)
names(fce.rondeau)
names(fce.rubio)
fce.grahl=rename(fce.grahl,c("NandN"="N.N"));#renaming NOx column to be consistent with other datasets
fce.wq=rbind(fce.grahl,fce.losada,fce.rondeau,fce.rubio)
fce.wq$Date.EST=date.fun(as.character(fce.wq$DATE),form="%Y-%m-%d")

unique(fce.wq$SITENAME)
fce.wq.sites=data.frame(SITENAME=c(paste0("SRS",c("1a","1c","1d",2,3,4,5,6)),paste0("TS/PH",c("1a","1b",2,3,"6b","6a","7a","7b"))),
                        Station.ID=c(paste0("SRS",c("1a","1c","1d",2,3,4,5,6)),paste0("TS/PH",c("1a","1b",2,3,"6b","6a","7a","7b"))),
                        Region=c(rep("SRS",8),rep("TS",8)))

fce.wq=merge(fce.wq,fce.wq.sites,"SITENAME")
unique(fce.wq$SITENAME)

fce.wq$TP=round((fce.wq$TP*P.mw),4);#convert uM concentration to ug/L
fce.wq$TN=round((fce.wq$TN*N.mw)*0.001,2);#convert uM concentration to mg/L
fce.wq$TOC=round((fce.wq$TOC*C.mw)*0.001,2);#convert uM concentration to mg/L
fce.wq$DOC=round((fce.wq$DOC*C.mw)*0.001,2);#convert uM concentration to mg/L
fce.wq$SRP=round((fce.wq$SRP*P.mw),4);#convert uM concentration to ug/L
fce.wq$SRP=with(fce.wq,ifelse(SRP<=0.3,0.3097,SRP))
fce.wq$N.N=with(fce.wq,ifelse(N.N<=0,0.01,N.N))
fce.wq$DIN=with(fce.wq,round(((N.N+NH4)*N.mw)*0.001,4));#convert uM concentration to mg/L

vars3=c("Station.ID","Region","Date.EST","TP","TN","TOC","DOC","SRP","DIN")
fce.wq=fce.wq[,vars3]
head(fce.wq)

grab.wq=fce.wq

##