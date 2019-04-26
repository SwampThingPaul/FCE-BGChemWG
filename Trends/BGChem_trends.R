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

#Helper variables 
N.mw=14.0067
P.mw=30.973762
C.mw=12.0107

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
fce.losada$Dataset.ID=dat.sources$Dataset[1]
fce.grahl=with(dat.sources[2,],read.lter("knb-lter-fce/",PASTA,DOI))
fce.grahl$Dataset.ID=dat.sources$Dataset[2]
fce.rubio=with(dat.sources[3,],read.lter("knb-lter-fce/",PASTA,DOI))
fce.rubio$Dataset.ID=dat.sources$Dataset[3]
fce.rondeau=with(dat.sources[4,],read.lter("knb-lter-fce/",PASTA,DOI))
fce.rondeau$Dataset.ID=dat.sources$Dataset[4]

names(fce.grahl)
names(fce.losada)
names(fce.rondeau)
names(fce.rubio)
fce.grahl=rename(fce.grahl,c("NandN"="N.N"));#renaming NOx column to be consistent with other datasets
fce.wq=rbind(fce.grahl,fce.losada,fce.rondeau,fce.rubio)
fce.wq$Date.EST=date.fun(as.character(fce.wq$DATE),form="%Y-%m-%d")
fce.wq=subset(fce.wq,Date.EST<date.fun("2018-05-01"))

plot(TP~SRP,fce.wq,ylim=c(0,4),xlim=c(0,2.5));abline(0,1,col="red")

unique(fce.wq$SITENAME)
fce.wq$SITENAME=with(fce.wq,ifelse(substring(fce.wq$SITENAME,1,2)=="TS",paste0("TS/PH",substring(fce.wq$SITENAME,6,7)),as.character(fce.wq$SITENAME)))

fce.wq.sites=data.frame(SITENAME=c(paste0("SRS",c("1a","1c","1d",2,3,4,5,6)),paste0("TS/PH",c("1a","1b",2,3,"6b","6a","7a","7b"))),
                        Station.ID=c(paste0("SRS",c("1a","1c","1d",2,3,4,5,6)),paste0("TS/PH",c("1a","1b",2,3,"6b","6a","7a","7b"))),
                        ALIAS=c(paste0("SRS",c(1,1,1,2,3,4,5,6)),paste0("TS/PH",c("1a","1b",2,3,"6b","6a","7a","7b"))),
                        Region=c(rep("SRS",8),rep("TS",8)))

fce.wq=merge(fce.wq,fce.wq.sites,"SITENAME")
unique(fce.wq$SITENAME)

#fce.wq$SRPFlag=with(fce.wq,ifelse(is.na(SRP)|SRP==0,1,0));
#fce.wq$TPFlag=with(fce.wq,ifelse(is.na(TP),1,0));
fce.wq$TPReversal.abs=with(fce.wq,ifelse(is.na(SRP)==T|is.na(TP)==T,0,ifelse(SRP>TP,1,0)))
fce.wq$TPReversal=with(fce.wq,ifelse(is.na(SRP)==T|is.na(TP)==T,0,ifelse(SRP>(TP*1.3),1,0)));# Reversals identified as 1 reversals consistent with TP rule evaluation
fce.wq$OCReversal=with(fce.wq,ifelse(is.na(DOC)==T|is.na(TOC)==T,0,ifelse(DOC>(TOC*1.3),1,0)));
fce.wq$TNReversal=with(fce.wq,ifelse(is.na(N.N+NH4)==T|is.na(TN)==T,0,ifelse((N.N+NH4)>(TN*1.3),1,0)));
sum(fce.wq$TPReversal.abs,na.rm=T)
sum(fce.wq$TPReversal,na.rm=T)
sum(fce.wq$OCReversal,na.rm=T)
sum(fce.wq$TNReversal,na.rm=T)
#write.csv(fce.wq,paste0(export.path,"fce_biogeochem.csv"),row.names=F)

# Check
par(family="serif",oma=c(1,1,1,1),mar=c(4,4,1,1))
plot(TP~SRP,fce.wq,ylim=c(0,4),xlim=c(0,2.5),ylab="TP (\u03BCM)",xlab="SRP (\u03BCM)",pch=21,bg=ifelse(TPReversal==1,"red",NA),col=adjustcolor("grey",0.8));abline(0,1,col="red")
#

fce.wq$TP.ugL=round((fce.wq$TP*P.mw),4);#convert uM concentration to ug/L
fce.wq$TN.mgL=round((fce.wq$TN*N.mw)*0.001,2);#convert uM concentration to mg/L
fce.wq$TOC.mgL=round((fce.wq$TOC*C.mw)*0.001,2);#convert uM concentration to mg/L
fce.wq$DOC.mgL=round((fce.wq$DOC*C.mw)*0.001,2);#convert uM concentration to mg/L
fce.wq$SRP=with(fce.wq,ifelse(SRP<0.01,0.01,SRP))
fce.wq$SRP.ugL=round((fce.wq$SRP*P.mw),4);#convert uM concentration to ug/L
fce.wq$N.N=with(fce.wq,ifelse(N.N<=0,0.01,N.N))
fce.wq$DIN.mgL=with(fce.wq,round(((N.N+NH4)*N.mw)*0.001,4));#convert uM concentration to mg/L
summary(fce.wq)

vars3=c("Dataset.ID","Station.ID","ALIAS","Region","Date.EST","TP.ugL","TN.mgL","TOC.mgL","DOC.mgL","SRP.ugL","DIN.mgL","TPReversal","OCReversal","TNReversal")
grab.wq=fce.wq[,vars3]

grab.wq$month=format(grab.wq$Date.EST,"%m")
grab.wq$CY=as.numeric(format(grab.wq$Date.EST,"%Y"))
grab.wq$WY=WY(grab.wq$Date.EST)
grab.wq$dec.year=decimal_date(grab.wq$Date.EST)
grab.wq=grab.wq[order(grab.wq$Station.ID,grab.wq$Date.EST),]
#grab.wq$SRP.ugL=with(grab.wq,ifelse(SRP.ugL<=0,"2",SRP));#sets minimum to a standard MDL

grab.wq$TP.ugL.scn=with(grab.wq,ifelse(TPReversal==1,NA,TP.ugL))
grab.wq$SRP.ugL.scn=with(grab.wq,ifelse(TPReversal==1,NA,SRP.ugL))
grab.wq$TN.mgL.scn=with(grab.wq,ifelse(TNReversal==1,NA,TN.mgL))
grab.wq$DIN.mgL.scn=with(grab.wq,ifelse(TNReversal==1,NA,DIN.mgL))
grab.wq$TOC.mgL.scn=with(grab.wq,ifelse(OCReversal==1,NA,TOC.mgL))
grab.wq$DOC.mgL.scn=with(grab.wq,ifelse(OCReversal==1,NA,DOC.mgL))


## Annual average
site.val=c(paste0("SRS",c(1,2,3,4,5,6)),paste0("TS/PH",c("1a",2,3,"6a","7a")))

WQ.WY=ddply(grab.wq,c("Station.ID","ALIAS","Region","WY"),summarise,mean.TP=mean(TP.ugL.scn,na.rm=T),N.TP=N(TP.ugL.scn),mean.SRP=mean(SRP.ugL.scn,na.rm=T),N.SRP=N(SRP.ugL.scn),mean.TN=mean(TN.mgL.scn,na.rm=T),N.TN=N(TN.mgL.scn),mean.DIN=mean(DIN.mgL.scn,na.rm=T),N.DIN=N(DIN.mgL.scn),mean.DOC=mean(DOC.mgL.scn,na.rm=T),N.DOC=N(DOC.mgL.scn),mean.TOC=mean(TOC.mgL.scn,na.rm=T),N.DOC=N(TOC.mgL.scn))
WQ.WY$dec.yr=with(WQ.WY,WY+0.328);#decimal date for Apirl 30th 
WQ.WY=subset(WQ.WY,ALIAS%in%site.val)

#building an uniform dataset for ALIAS=="SRS1"
WQ.WY$dup=with(WQ.WY,ifelse(Station.ID=="SRS1c"&WY==2005|Station.ID=="SRS1d"&WY==2006,1,0))
WQ.WY=subset(WQ.WY,dup==0)[,-ncol(WQ.WY)]


## Plot

## Time Series plot
unique(WQ.WY$Station.ID)
site.val=c(paste0("SRS",c(1,2,3,4,5,6)),paste0("TS/PH",c("1a",2,3,"6a","7a")))
WQ.WY$ALIAS=factor(WQ.WY$ALIAS,level=site.val)
#tiff(filename=paste0(plot.path,"TPTNDOC_WQPlots.tiff"),width=9,height=6.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
#png(filename=paste0(plot.path,"TPTNDOC_WQPlots.png"),width=9,height=6.5,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(2.5,1.75,1,0.25),mar=c(1.5,2,0.5,1))
layout(matrix(1:48,6,8,byrow=F),widths=c(1,1,0.3,1,1,0.30,1,1))

xlim.val=c(2001.8,2018.5);by.x=5;xmaj=seq(round(xlim.val[1]),xlim.val[2],by.x);xmin=seq(round(xlim.val[1]),xlim.val[2],by.x/by.x)
ylim.val=c(1,110);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
for(i in 1:length(site.val)){
  plot(mean.TP~WY,WQ.WY,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n",log="y")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(grab.wq,ALIAS==site.val[i]),points(jitter(WY,1),TP.ugL.scn,pch=19,col=adjustcolor("grey",0.75)))
  with(subset(WQ.WY,ALIAS==site.val[i]),points(WY,mean.TP,pch=21,bg=ifelse(N.TP>=4,"indianred1",NA),cex=1.25,lwd=0.2))
  if(i==6|i==11){axis_fun(1,line=-0.5,xmaj,xmin,xmaj)}else(axis_fun(1,xmaj,xmin,NA))
  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=3,site.val[i],cex=0.8)
  if(i==4){text(x=1993.5,y=ylim.val[1],"Total Phosphorus (\u03BCg L\u207B\u00B9)",xpd=NA,srt=90,cex=1.75,adj=0)}
}
for(j in 1:7){plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)}

ylim.val=c(0.1,5);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
for(i in 1:length(site.val)){
  plot(mean.TN~WY,WQ.WY,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n",log="y")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(grab.wq,ALIAS==site.val[i]),points(jitter(WY,1),TN.mgL.scn,pch=19,col=adjustcolor("grey",0.75)))
  with(subset(WQ.WY,ALIAS==site.val[i]),points(WY,mean.TN,pch=22,bg=ifelse(N.TN>=4,"dodgerblue1",NA),cex=1.25,lwd=0.2))
  if(i==6|i==11){axis_fun(1,line=-0.5,xmaj,xmin,xmaj)}else(axis_fun(1,xmaj,xmin,NA))
  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=3,site.val[i],cex=0.8)
  if(i==4){text(x=1993.5,y=ylim.val[1],"Total Nitrogen (mg L\u207B\u00B9)",xpd=NA,srt=90,cex=1.75,adj=0)}
}
for(j in 1:7){plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)}

ylim.val=c(1,55);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
for(i in 1:length(site.val)){
  plot(mean.DOC~WY,WQ.WY,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n",log="y")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(grab.wq,ALIAS==site.val[i]),points(jitter(WY,1),DOC.mgL.scn,pch=19,col=adjustcolor("grey",0.75)))
  with(subset(WQ.WY,ALIAS==site.val[i]),points(WY,mean.DOC,pch=23,bg=ifelse(N.DOC>=4,"darkgoldenrod1",NA),cex=1.25,lwd=0.2))
  if(i==6|i==11){axis_fun(1,line=-0.5,xmaj,xmin,xmaj)}else(axis_fun(1,xmaj,xmin,NA))
  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=3,site.val[i],cex=0.8)
  if(i==5){text(x=1993.5,y=ylim.val[2],"Dissolved Organic Carbon (mg L\u207B\u00B9)",xpd=NA,srt=90,cex=1.75,adj=0)}
}
mtext(side=1,line=0.5,outer=T,"Water Year")
plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)
legend.text=c("Annual Mean TP", "Annual Mean TN","Annual Mean DOC", "Grab Samples")
pt.cols=c("indianred1","dodgerblue1","darkgoldenrod1",adjustcolor("grey",0.75))
legend(0.5,0.75,legend=legend.text,pch =c(21,22,23,19),col=c("black","black","black",pt.cols[4]),lwd=0.2,lty=NA,pt.bg=pt.cols,pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()

#tiff(filename=paste0(plot.path,"SRPDIN_WQPlots.tiff"),width=7,height=6.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
#png(filename=paste0(plot.path,"SRPDIN_WQPlots.png"),width=7,height=6.5,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(2.5,1.75,1,0.25),mar=c(1.5,2,0.5,1))
#layout(matrix(1:12,6,2,byrow=F))
layout(matrix(1:30,6,5,byrow=F),widths=c(1,1,0.3,1,1))

xlim.val=c(2001.8,2018.5);by.x=5;xmaj=seq(round(xlim.val[1]),xlim.val[2],by.x);xmin=seq(round(xlim.val[1]),xlim.val[2],by.x/by.x)
ylim.val=c(0.1,65);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
for(i in 1:length(site.val)){
  plot(mean.SRP~WY,WQ.WY,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n",log="y")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(grab.wq,ALIAS==site.val[i]),points(jitter(WY,1),SRP.ugL.scn,pch=19,col=adjustcolor("grey",0.75)))
  with(subset(WQ.WY,ALIAS==site.val[i]),points(WY,mean.SRP,pch=21,bg=ifelse(N.SRP>=4,"sienna1",NA),cex=1.25,lwd=0.2))
  if(i==6|i==11){axis_fun(1,line=-0.5,xmaj,xmin,xmaj)}else(axis_fun(1,xmaj,xmin,NA))
  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=3,site.val[i],cex=0.8)
  if(i==5){text(x=1995.5,y=ylim.val[2],"Soluable Reactive Phosphorus (\u03BCg L\u207B\u00B9)",xpd=NA,srt=90,cex=1.75,adj=0)}
}
for(j in 1:7){plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)}

ylim.val=c(0.01,3);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
for(i in 1:length(site.val)){
  plot(mean.DIN~WY,WQ.WY,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n",log="y")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(grab.wq,ALIAS==site.val[i]),points(jitter(WY,1),DIN.mgL.scn,pch=19,col=adjustcolor("grey",0.75)))
  with(subset(WQ.WY,ALIAS==site.val[i]),points(WY,mean.DIN,pch=22,bg=ifelse(N.DIN>=4,"darkcyan",NA),cex=1.25,lwd=0.2))
  if(i==6|i==11){axis_fun(1,line=-0.5,xmaj,xmin,xmaj)}else(axis_fun(1,xmaj,xmin,NA))
  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=3,site.val[i],cex=0.8)
  if(i==5){text(x=1995,y=ylim.val[2],"Dissolved Inorganic Nitrogen (mg L\u207B\u00B9)",xpd=NA,srt=90,cex=1.75,adj=0)}
}
mtext(side=1,outer=T,"Water Year")
plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)
legend.text=c("Annual Mean SRP", "Annual Mean DIN", "Grab Samples")
pt.cols=c("sienna1","darkcyan",adjustcolor("grey",0.75))
legend(0.5,0.75,legend=legend.text,pch =c(21,22,19),col=c("black","black",pt.cols[3]),lwd=0.2,lty=NA,pt.bg=pt.cols,pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()