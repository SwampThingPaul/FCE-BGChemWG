## 
## Florida Coastal Everglades LTER
## Biogeochemistry Working Group
## bacterioplankton productivity
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
  dt1=read.csv(infile1,na.strings=c("-9999","-9999.00","-9999.000",-9999))
  return(dt1)
}

## Data Import
bact.prod=read.lter("knb-lter-fce/","1056/7/","c5e7b326feffc37b900f2e042b70a434")
bact.prod[bact.prod==-9999]=NA
bact.prod$Date.EST=date.fun(as.character(bact.prod$Date),form="%Y-%m-%d")
bact.prod$WY=WY(bact.prod$Date.EST)
bact.prod=subset(bact.prod,Date.EST<date.fun("2018-05-01"))

bact.prod$Bact.ab.cellsL=bact.prod$BacteriaAbundance*0.001

unique(bact.prod$SITE)
bact.prod$SITE=with(bact.prod,ifelse(substring(SITE,1,2)=="TS",paste0("TS/PH",substring(SITE,6,7)),as.character(SITE)))

fce.wq.sites=data.frame(SITE=c(paste0("SRS",c("1a","1c","1d",2,3,4,5,6)),paste0("TS/PH",c("1a","1b",2,3,"6b","6a","7a","7b"))),
                        Station.ID=c(paste0("SRS",c("1a","1c","1d",2,3,4,5,6)),paste0("TS/PH",c("1a","1b",2,3,"6b","6a","7a","7b"))),
                        ALIAS=c(paste0("SRS",c(1,1,1,2,3,4,5,6)),paste0("TS/PH",c("1a","1b",2,3,"6b","6a","7a","7b"))),
                        Region=c(rep("SRS",8),rep("TS",8)))
bact.prod=merge(bact.prod,fce.wq.sites,"SITE")

site.val=c(paste0("SRS",c(1,2,3,4,5,6)),paste0("TS/PH",c("1b",2,3,"6a","7a")))

WQ.WY=ddply(bact.prod,c("SITE","ALIAS","Region","WY"),summarise,mean.BP.da=mean(BacteriaProductionDailyRate,na.rm=T),N.BP.da=N(BacteriaProductionDailyRate),mean.BAbund=mean(Bact.ab.cellsL,na.rm=T),N.BAbund=N(Bact.ab.cellsL))
WQ.WY$dec.yr=with(WQ.WY,WY+0.328);#decimal date for Apirl 30th 
WQ.WY=subset(WQ.WY,ALIAS%in%site.val)

#building an uniform dataset for ALIAS=="SRS1"
WQ.WY$dup=with(WQ.WY,ifelse(SITE=="SRS1c"&WY==2005|SITE=="SRS1d"&WY==2006,1,0))
WQ.WY=subset(WQ.WY,dup==0)[,-ncol(WQ.WY)]

## Time Series plot
unique(WQ.WY$SITE)
site.val=c(paste0("SRS",c(1,2,3,4,5,6)),paste0("TS/PH",c("1b",2,3,"6a","7a")))
WQ.WY$ALIAS=factor(WQ.WY$ALIAS,level=site.val)
#tiff(filename=paste0(plot.path,"BP_WQPlots.tiff"),width=7,height=6.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",oma=c(2.5,1.75,1,0.25),mar=c(1.5,2,0.5,1))
layout(matrix(1:30,6,5,byrow=F),widths=c(1,1,0.3,1,1))

xlim.val=c(2001.8,2018.5);by.x=5;xmaj=seq(round(xlim.val[1]),xlim.val[2],by.x);xmin=seq(round(xlim.val[1]),xlim.val[2],by.x/by.x)
ylim.val=c(-20,200);by.y=100;ymaj=seq(max(0,ylim.val[1]),ylim.val[2],by.y);ymin=seq(max(0,ylim.val[1]),ylim.val[2],by.y/2)
for(i in 1:length(site.val)){
  plot(mean.BP.da~WY,WQ.WY,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  #abline(h=0)
  with(subset(bact.prod,ALIAS==site.val[i]),points(jitter(WY,1),BacteriaProductionDailyRate,pch=19,col=adjustcolor("grey",0.75)))
  with(subset(WQ.WY,ALIAS==site.val[i]),points(WY,mean.BP.da,pch=21,bg=ifelse(N.BP.da>=4,"sienna1",NA),cex=1.25,lwd=0.2))
  if(i==6|i==11){axis_fun(1,line=-0.5,xmaj,xmin,xmaj)}else(axis_fun(1,xmaj,xmin,NA))
  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=3,site.val[i],cex=0.8)
  if(i==5){text(x=1995.5,y=ylim.val[2],"Baterial Productivity (\u03BCg C L\u207B\u00B9 d\u207B\u00B9)",xpd=NA,srt=90,cex=1.75,adj=0)}
}
for(j in 1:7){plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)}

ylim.val=c(0,40000);by.y=10000;ymaj=seq(max(0,ylim.val[1]),ylim.val[2],by.y);ymin=seq(max(0,ylim.val[1]),ylim.val[2],by.y/2)
for(i in 1:length(site.val)){
  plot(mean.BAbund~WY,WQ.WY,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  abline(h=0)
  with(subset(bact.prod,ALIAS==site.val[i]),points(jitter(WY,1),Bact.ab.cellsL,pch=19,col=adjustcolor("grey",0.75)))
  with(subset(WQ.WY,ALIAS==site.val[i]),points(WY,mean.BAbund,pch=22,bg=ifelse(N.BAbund>=4,"darkcyan",NA),cex=1.25,lwd=0.2))
  if(i==6|i==11){axis_fun(1,line=-0.5,xmaj,xmin,xmaj)}else(axis_fun(1,xmaj,xmin,NA))
  axis_fun(2,ymaj,ymin,ymaj*0.001);box(lwd=1)
  mtext(side=3,site.val[i],cex=0.8)
  if(i==5){text(x=1995.5,y=ylim.val[2],"Bacterial Abundance (x10\u207B\u00B3 cells L\u207B\u00B9)",xpd=NA,srt=90,cex=1.75,adj=0)}
}
mtext(side=1,line=0.5,outer=T,"Water Year")
plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)
legend.text=c("Annual Mean\nBacterioplankton Productivity","Annual Mean\nBacterial Abundance", "Grab Samples")
pt.cols=c("sienna1","darkcyan",adjustcolor("grey",0.75))
legend(0.5,1,legend=legend.text,pch =c(21,22,19),col=c("black","black",pt.cols[3]),lwd=0.2,lty=NA,pt.bg=pt.cols,pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1.5,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()