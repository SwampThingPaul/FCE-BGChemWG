rm(list=ls(all=T));cat("\014");dev.off()

#Libraries
library(AnalystHelper);#devtools::install_github("SwampThingPaul/AnalystHelper")
library(plyr)
library(reshape)
library(zoo);
library(shiny)

library(tmap)
library(rgdal)

#Paths
wd="D:/_GitHub/FCE-BGChemWG/Trends"

GIS.path=paste0(wd,"/GIS")

#Helper variables 
N.mw=14.0067
P.mw=30.973762
C.mw=12.0107
utm17=CRS("+proj=utm +zone=17 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# custom functions
read.lter=function(data.package,PASTA,DOI){
  prefix="http://pasta.lternet.edu/package/data/eml/"
  
  infile1=paste0(prefix,data.package,PASTA,DOI)
  dt1=read.csv(infile1,na.strings=c("-9999","-9999.00","-9999.000"))
  return(dt1)
}

## GIS data
tmap_mode("view")
FCE.LTER.sites=spTransform(readOGR(GIS.path,"ltersites_utm"),utm17)
region.sites=data.frame(SITE=c("SRS-1a","SRS-1b", "SRS-1c", "SRS-1d", "SRS-2", "SRS-3", "SRS-4", "SRS-5", "SRS-6",  "TS/Ph-1a", "TS/Ph-1b", "TS/Ph-2a","TS/Ph-2b", "TS/Ph-3", "TS/Ph-6a", "TS/Ph-6b","TS/Ph-7a", "TS/Ph-7b","TS/Ph-4", "TS/Ph-5", "TS/Ph-8", "TS/Ph-9","TS/Ph-10", "TS/Ph-11"),
                        Region=c(rep("Shark River Slough",9),rep("Taylor Slough",9),rep("Panhandle",3),rep("Florida Bay",3)))
FCE.LTER.sites=merge(FCE.LTER.sites,region.sites,"SITE")

ENP=spTransform(readOGR(GIS.path,"ENP_boundary"),utm17)

## Data Import
dat.sources=data.frame(Dataset=c("LT_ND_Losada_002","LT_ND_Grahl_002","LT_ND_Rubio_002","LT_ND_Rondeau_002"),
                       PASTA=c("1075/8/","1073/12/","1080/8/","1077/3/"),
                       DOI=c("ac7159e66cbad75abb61bd1992f8d2c0","15e6ff92f875a272ba6d98db3f738026","3a84ab2009eda6b73e375aa3ad56da1a","832e5493edf5a5d79d8a60535e26a012"))

fce.losada=with(dat.sources[1,],read.lter("knb-lter-fce/",PASTA,DOI))
fce.losada$Dataset.ID=dat.sources$Dataset[1]
fce.grahl=with(dat.sources[2,],read.lter("knb-lter-fce/",PASTA,DOI))
fce.grahl$Dataset.ID=dat.sources$Dataset[2]
fce.rubio=with(dat.sources[3,],read.lter("knb-lter-fce/",PASTA,DOI))
fce.rubio$Dataset.ID=dat.sources$Dataset[3]
fce.rondeau=with(dat.sources[4,],read.lter("knb-lter-fce/",PASTA,DOI))
fce.rondeau$Dataset.ID=dat.sources$Dataset[4]

fce.grahl=rename(fce.grahl,c("NandN"="N.N"));#renaming NOx column to be consistent with other datasets

fce.wq=rbind(fce.grahl,fce.losada,fce.rondeau,fce.rubio)
fce.wq$Date.EST=date.fun(as.character(fce.wq$DATE),form="%Y-%m-%d")
fce.wq$WY=WY(fce.wq$Date.EST)
fce.wq$dec.year=decimal_date(fce.wq$Date.EST)
fce.wq$hydro.season=FL.Hydroseason(fce.wq$Date.EST)
fce.wq$SITENAME=with(fce.wq,ifelse(substring(fce.wq$SITENAME,1,2)=="TS",paste0("TS/PH",substring(fce.wq$SITENAME,6,7)),as.character(fce.wq$SITENAME)))

fce.wq.sites=data.frame(SITENAME=c(paste0("SRS",c("1a","1c","1d",2,3,4,5,6)),paste0("TS/PH",c("1a","1b",2,3,"6b","6a","7a","7b")),paste0("TS/PH",c(4,5,8)),paste0("TS/PH",c(9,10,11))),
                        Region=c(rep("SRS",8),rep("TS",8),rep("Ph",3),rep("FB","3")))
fce.wq=merge(fce.wq,fce.wq.sites,all.x=T)
fce.wq=subset(fce.wq,Date.EST<date.fun("2018-05-01"))

## Parameter reversals
fce.wq$TPReversal.abs=with(fce.wq,ifelse(is.na(SRP)==T|is.na(TP)==T,0,ifelse(SRP>TP,1,0)))
fce.wq$TPReversal=with(fce.wq,ifelse(is.na(SRP)==T|is.na(TP)==T,0,ifelse(SRP>(TP*1.3),1,0)));# Reversals identified as 1 reversals consistent with TP rule evaluation
fce.wq$OCReversal=with(fce.wq,ifelse(is.na(DOC)==T|is.na(TOC)==T,0,ifelse(DOC>(TOC*1.3),1,0)));
fce.wq$TNReversal=with(fce.wq,ifelse(is.na(N.N+NH4)==T|is.na(TN)==T,0,ifelse((N.N+NH4)>(TN*1.3),1,0)));

##Convert uM concentrations to ug/L (P-species) or mg/L (N and C species)
fce.wq$TP=round((fce.wq$TP*P.mw),4);#convert uM concentration to ug/L
fce.wq$TN=round((fce.wq$TN*N.mw)*0.001,2);#convert uM concentration to mg/L
fce.wq$TOC=round((fce.wq$TOC*C.mw)*0.001,2);#convert uM concentration to mg/L
fce.wq$DOC=round((fce.wq$DOC*C.mw)*0.001,2);#convert uM concentration to mg/L
fce.wq$SRP=with(fce.wq,ifelse(SRP<0.01,0.01,SRP))
fce.wq$SRP=round((fce.wq$SRP*P.mw),4);#convert uM concentration to ug/L
fce.wq$N.N=with(fce.wq,ifelse(N.N<=0,0.01,N.N))
fce.wq$DIN=with(fce.wq,N.N+NH4)
fce.wq$DIN=with(fce.wq,round((DIN*N.mw)*0.001,4));#convert uM concentration to mg/L

names(fce.wq)
id.vars.val=c("SITENAME","Region","Dataset.ID","Date.EST","WY","hydro.season","dec.year")
param.vars.val=c("TP","SRP","TN","DIN","DOC")

fce.wq.stack=melt(fce.wq[,c(id.vars.val,param.vars.val)],id.vars=id.vars.val,variable_name = "parameter")

rev.list=list(fce.wq[,c("SITENAME","Date.EST","TPReversal")],fce.wq[,c("SITENAME","Date.EST","TNReversal")],fce.wq[,c("SITENAME","Date.EST","OCReversal")])
for(i in 1:3){
colnames(rev.list[[i]])=c("SITENAME","Date.EST","Reversal")
}
names(rev.list)=c("TP","TN","OC")

param.list=c("TP","SRP","TN","DIN","DOC")
rev.param.list=c("TP","TP","TN","TN","OC")

fce.wq.stack.qa=data.frame()
for(i in 1:length(param.list)){
  tmp=subset(fce.wq.stack,parameter==param.list[i])
  tmp=merge(tmp,rev.list[[rev.param.list[i]]],c("SITENAME","Date.EST"))
  fce.wq.stack.qa=rbind(fce.wq.stack.qa,tmp)
}
fce.wq.stack.qa=subset(fce.wq.stack.qa,is.na(value)==F)

WQ.WY=ddply(subset(fce.wq.stack.qa,Reversal==0),c("SITENAME","WY","Region","parameter"),summarise,mean.val=mean(value,na.rm=T),GM.val=exp(mean(log(value),na.rm=T)),N.val=N(value))
WQ.WY$dec.yr=with(WQ.WY,WY+0.328);#decimal date for Apirl 30th 

## for the shiny
param="TP"
site.val="SRS2"
WY.val=2007

tmp.da.dat=subset(fce.wq.stack.qa,parameter==param&SITENAME==site.val&WY==WY.val)

xlim.val=c(date.fun(c(paste0(WY.val-1,"-05-01"),paste0(WY.val,"-05-01"))))
