## 
## Florida Coastal Everglades LTER
## Biogeochemistry Working Group
## sawgrass AG and BG
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

paths=paste0(getwd(),c("/Exports/","/Plots/","/Data/"))
#Folder.Maker(paths);#One and done. Creates folders in working directory.
export.path=paths[1]
plot.path=paths[2]
data.path=paths[3]

#Helper variables 
N.mw=14.0067
P.mw=30.973762
C.mw=12.0107

# custom functions
read.lter=function(data.package,PASTA,DOI,data.set.ID){
  prefix="http://pasta.lternet.edu/package/data/eml/"
  
  infile1=paste0(prefix,data.package,PASTA,DOI)
  dt1=read.csv(infile1,na.strings=c("-9999","-9999.00","-9999.000"))
  dt1$Dataset.ID=data.set.ID
  return(dt1)
}

## Data Import
fce.wq.sites=data.frame(SITENAME=c(paste0("SRS",c("1a","1b","1c","1d",2,3,4,5,6)),paste0("TS/PH",c("1a","1b",2,3,"6b","6a","7a","7b",4,5))),
                        Station.ID=c(paste0("SRS",c("1a","1b","1c","1d",2,3,4,5,6)),paste0("TS/PH",c("1a","1b",2,3,"6b","6a","7a","7b",4,5))),
                        ALIAS=c(paste0("SRS",c(1,1,1,1,2,3,4,5,6)),paste0("TS/PH",c("1","1",2,3,"6b","6a","7a","7b",4,5))),
                        Region=c(rep("SRS",9),rep("TS",8),rep("Ph",2)))


# Datasets
#dat.sources=data.frame(Dataset=c("LT_ND_Grahl_004"),PASTA=c("1070/9/"),DOI=c("a284cfad90ff5c7d8b6623128b0491fd"))
dat.sources=rbind(dat.sources,
                  data.frame(Dataset=c("LT_ND_Rubio_005"),PASTA=c("1083/7/"),DOI=c("0759acbb8d9136a4f2123aa29009139f")))
dat.sources=rbind(dat.sources,
                  data.frame(Dataset=c("LT_ND_Grahl_005"),PASTA=c("1071/9/"),DOI=c("dcbf60148d6e9e6d0c45a916d0fdfeff")))
dat.sources=rbind(dat.sources,
                  data.frame(Dataset=c("LT_ND_Rubio_004"),PASTA=c("1082/9/"),DOI=c("6accea872c3cc53add878bf3bafabbec")))

sg.grahl004=with(dat.sources[1,],read.lter("knb-lter-fce/",PASTA,DOI,Dataset))
sg.rubio005=with(dat.sources[2,],read.lter("knb-lter-fce/",PASTA,DOI,Dataset))
sg.grahl005=with(dat.sources[3,],read.lter("knb-lter-fce/",PASTA,DOI,Dataset))
sg.rubio004=with(dat.sources[4,],read.lter("knb-lter-fce/",PASTA,DOI,Dataset))

###


#Data formatting
sg.CN=rbind(sg.grahl004,sg.rubio005)
sg.CN$Date=date.fun(sg.CN$Date,form="%Y-%m-%d")
names(sg.CN)

plot(mgC.g~Date,subset(sg.CN,SITENAME=="TS/Ph3"&Type=="Aboveground"))

idvars=c("RecordNum", "SITENAME", "Date", "Type","Dataset.ID")
sg.CN.melt=melt(sg.CN,id.vars=idvars)

plot(value~Date,subset(sg.CN.melt,variable=="mgC.g"&SITENAME=="TS/Ph3"&Type=="Aboveground"))

sg.P=rbind(sg.grahl005,sg.rubio004)
sg.P=rename(sg.P,c("TypeTP"="Type"))
sg.P$Date=date.fun(sg.P$Date,form="%Y-%m-%d")
sg.P.melt=melt(sg.P,id.vars=idvars)

sg.CNP.melt=rbind(sg.CN.melt,sg.P.melt)
sg.CNP.melt$Type[sg.CNP.melt$Type=="Aboveground"]="Above"
sg.CNP.melt$Type[sg.CNP.melt$Type=="Belowground"]="Below"
plot(value~Date,subset(sg.CNP.melt,variable=="mgC.g"&SITENAME=="TS/Ph3"&Type=="Above"))

sg.CNP.melt$SITENAME=with(sg.CNP.melt,ifelse(substring(SITENAME,1,2)=="TS",paste0("TS/PH",substring(SITENAME,6,7)),as.character(SITENAME)))
sg.CNP.melt=merge(sg.CNP.melt,fce.wq.sites,"SITENAME",all.x=T)
sum(is.na(sg.CNP.melt$ALIAS))
plot(value~Date,subset(sg.CNP.melt,variable=="mgC.g"&SITENAME=="TS/PH3"&Type=="Above"))

sg.CNP=data.frame(cast(sg.CNP.melt,ALIAS+Date+Type~variable,value="Value",fun.aggregate = function(x)mean(x,na.rm=T)))
sg.CNP$CY=as.numeric(format(sg.CNP$Date,"%Y"))
ddply(sg.CNP,c("ALIAS","Type","CY"),summarise,N.val=N(mgC.g))

unique(sg.CNP$ALIAS)
plot(mgC.g~Date,subset(sg.CNP,ALIAS=="TS/PH3"&Type=="Above"))

cols=c(rgb(237/255,125/255,49/255),rgb(0/255,112/255,192/255),rgb(0/255,176/255,80/255),rgb(255/255,192/255,0/255))

#tiff(filename=paste0(plot.path,"SRS_sawgrass.tiff"),width=4.3,height=6,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",oma=c(1,2,1,0.25),mar=c(3,2,0.5,1))
layout(matrix(1:6,3,2,byrow=F),widths=c(1,0.3))

xlim.val=date.fun(c("2002-09-01","2018-03-01"));xmaj=seq(xlim.val[1],xlim.val[2],"1 year");xmin=seq(xlim.val[1],xlim.val[2],"1 year")

site.vals=c("SRS1","SRS2","SRS3")
ylim.val=c(350,550);by.y=50;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y)
plot(mgC.g~Date,sg.CNP,ylim=ylim.val,xlim=xlim.val,yaxs="i",xaxs="i",yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n")
axis_fun(2,ymaj,ymin,ymaj);axis_fun(1,xmaj,xmin,NA);box(lwd=1)
text(xmaj,ylim.val[1]-by.y*0.35,format(xmaj,"%Y"),xpd=NA,srt=90,adj=1)
for(i in 1:3){
  with(subset(sg.CNP,ALIAS==site.vals[i]&Type=="Above"),lines(Date,mgC.g,lty=1,lwd=2,col=cols[i]))
  with(subset(sg.CNP,ALIAS==site.vals[i]&Type=="Below"),lines(Date,mgC.g,lty=3,lwd=2,col=cols[i]))
}
mtext(side=2,line=2.5,"mg C/g")

ylim.val=c(0,15);by.y=5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y)
plot(mgC.g~Date,sg.CNP,ylim=ylim.val,xlim=xlim.val,yaxs="i",xaxs="i",yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n")
axis_fun(2,ymaj,ymin,ymaj);axis_fun(1,xmaj,xmin,NA);box(lwd=1)
text(xmaj,ylim.val[1]-by.y*0.3,format(xmaj,"%Y"),xpd=NA,srt=90,adj=1)
for(i in 1:3){
  with(subset(sg.CNP,ALIAS==site.vals[i]&Type=="Above"),lines(Date,mgN.g,lty=1,lwd=2,col=cols[i]))
  with(subset(sg.CNP,ALIAS==site.vals[i]&Type=="Below"),lines(Date,mgN.g,lty=3,lwd=2,col=cols[i]))
}
mtext(side=2,line=2.5,"mg N/g")

ylim.val=c(0,1000);by.y=200;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y)
plot(mgC.g~Date,sg.CNP,ylim=ylim.val,xlim=xlim.val,yaxs="i",xaxs="i",yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n")
axis_fun(2,ymaj,ymin,ymaj);axis_fun(1,xmaj,xmin,NA);box(lwd=1)
text(xmaj,ylim.val[1]-by.y*0.35,format(xmaj,"%Y"),xpd=NA,srt=90,adj=1)
for(i in 1:3){
  with(subset(sg.CNP,ALIAS==site.vals[i]&Type=="Above"),lines(Date,ugP.g,lty=1,lwd=2,col=cols[i]))
  with(subset(sg.CNP,ALIAS==site.vals[i]&Type=="Below"),lines(Date,ugP.g,lty=3,lwd=2,col=cols[i]))
}
mtext(side=2,line=2.5,"\u03BCg P/g")

plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)
legend.text=sort(c(paste(site.vals,"AG"),paste(site.vals,"BG")))
legend(0.25,0.5,legend=legend.text,col=c(rep(cols[1],2),rep(cols[2],2),rep(cols[3],2)),lwd=2,lty=c(1,3),ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)
plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)
dev.off()

#tiff(filename=paste0(plot.path,"TS_sawgrass.tiff"),width=4.3,height=6,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",oma=c(1,2,1,0.25),mar=c(3,2,0.5,1))
layout(matrix(1:6,3,2,byrow=F),widths=c(1,0.3))

xlim.val=date.fun(c("2002-03-01","2018-03-01"));xmaj=seq(xlim.val[1],xlim.val[2],"1 year");xmin=seq(xlim.val[1],xlim.val[2],"1 year")

site.vals=c("TS/PH1","TS/PH2","TS/PH3","TS/PH6b")
ylim.val=c(350,550);by.y=50;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y)
plot(mgC.g~Date,sg.CNP,ylim=ylim.val,xlim=xlim.val,yaxs="i",xaxs="i",yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n")
axis_fun(2,ymaj,ymin,ymaj);axis_fun(1,xmaj,xmin,NA);box(lwd=1)
text(xmaj,ylim.val[1]-by.y*0.35,format(xmaj,"%Y"),xpd=NA,srt=90,adj=1)
for(i in 1:4){
  with(subset(sg.CNP,ALIAS==site.vals[i]&Type=="Above"&is.na(mgC.g)==F),lines(Date,mgC.g,lty=1,lwd=2,col=cols[i]))
  with(subset(sg.CNP,ALIAS==site.vals[i]&Type=="Below"&is.na(mgC.g)==F),lines(Date,mgC.g,lty=3,lwd=2,col=cols[i]))
}
mtext(side=2,line=2.5,"mg C/g")

ylim.val=c(0,15);by.y=5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y)
plot(mgC.g~Date,sg.CNP,ylim=ylim.val,xlim=xlim.val,yaxs="i",xaxs="i",yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n")
axis_fun(2,ymaj,ymin,ymaj);axis_fun(1,xmaj,xmin,NA);box(lwd=1)
text(xmaj,ylim.val[1]-by.y*0.3,format(xmaj,"%Y"),xpd=NA,srt=90,adj=1)
for(i in 1:4){
  with(subset(sg.CNP,ALIAS==site.vals[i]&Type=="Above"&is.na(mgN.g)==F),lines(Date,mgN.g,lty=1,lwd=2,col=cols[i]))
  with(subset(sg.CNP,ALIAS==site.vals[i]&Type=="Below"&is.na(mgN.g)==F),lines(Date,mgN.g,lty=3,lwd=2,col=cols[i]))
}
mtext(side=2,line=2.5,"mg N/g")

ylim.val=c(0,1000);by.y=200;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y)
plot(mgC.g~Date,sg.CNP,ylim=ylim.val,xlim=xlim.val,yaxs="i",xaxs="i",yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n")
axis_fun(2,ymaj,ymin,ymaj);axis_fun(1,xmaj,xmin,NA);box(lwd=1)
text(xmaj,ylim.val[1]-by.y*0.35,format(xmaj,"%Y"),xpd=NA,srt=90,adj=1)
for(i in 1:4){
  with(subset(sg.CNP,ALIAS==site.vals[i]&Type=="Above"&is.na(ugP.g)==F),lines(Date,ugP.g,lty=1,lwd=2,col=cols[i]))
  with(subset(sg.CNP,ALIAS==site.vals[i]&Type=="Below"&is.na(ugP.g)==F),lines(Date,ugP.g,lty=3,lwd=2,col=cols[i]))
}
mtext(side=2,line=2.5,"\u03BCg P/g")

plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)
legend.text=sort(c(paste(site.vals,"AG"),paste(site.vals,"BG")))
legend(0.25,0.5,legend=legend.text,col=c(rep(cols[1],2),rep(cols[2],2),rep(cols[3],2),rep(cols[4],2)),lwd=2,lty=c(1,3),ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)
plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)
dev.off()

#Periphyton
peri=read.lter("knb-lter-fce/","1107/9/","ac394fca14259073329ffe19ecf096f6","LT_PP_Gaiser_005")
