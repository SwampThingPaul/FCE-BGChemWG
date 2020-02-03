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
fce.wq=subset(fce.wq,Date.EST<date.fun("2019-05-01"))

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

plot(TP.ugL.scn~Date.EST,grab.wq)
head(grab.wq)


grab.wq.month=ddply(subset(grab.wq,ALIAS=="SRS6"),c("ALIAS","CY","month"),summarise,mean.TP=mean(TP.ugL.scn,na.rm=T),mean.SRP=mean(SRP.ugL.scn,na.rm=T))
grab.wq.month=subset(grab.wq.month,CY>=2001)

date.fill=seq(date.fun("2001-01-01"),date.fun("2018-12-01"),"months")
grab.wq.month=merge(grab.wq.month,data.frame(month=format(date.fill,"%m"),CY=format(date.fill,"%Y"),fill=1),all.y=T)
grab.wq.month$month.CY=with(grab.wq.month,date.fun(paste(CY,month,"01",sep="-")))
grab.wq.month

#tiff(filename=paste0(plot.path,"SRS6_TP_TS.tiff"),width=9,height=4,units="in",res=200,type="windows",bg="white",compression=c("lzw"))
#png(filename=paste0(plot.path,"SRS6_TP_TS.png"),width=9,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(2.5,1.75,1,0.25),mar=c(1.5,2,0.5,1))

xlim.val=date.fun(c("2001-05-01","2019-05-01"));xmaj=seq(xlim.val[1],xlim.val[2],"5 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
ylim.val=c(4,110);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")

plot(mean.TP~month.CY,grab.wq.month,axes=F,type="n",xlim=xlim.val,ylim=ylim.val,log="y")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(grab.wq.month,pt_line(month.CY,mean.TP,2,"dodgerblue1",1,21,"indianred1"))
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"))
axis_fun(2,ymaj,ymin,ymaj)
box(lwd=1)
mtext(side=1,line=1.5,"Date (Month-Year)")
mtext(side=2,line=2,"Total P (\u03BCg L\u207B\u00B9)")
text(xlim.val[1]-ddays(180),ylim.val[2]+30,"Monthly Surface Water TP concentration at SRS6",xpd=NA,adj=0,cex=1.25)
dev.off()

library(ggplot2)
ggplot(grab.wq.month)+
  geom_path(aes(month.CY,mean.SRP),color="indianred",size=1)+
  scale_y_continuous("Total P (\u03BCg L\u207B\u00B9)",breaks=seq(0,40,5))+
  theme_bw()


ggplot(grab.wq.month)+
  geom_point(mapping=aes(month.CY,mean.TP),shape=21,bg="grey")+
  geom_path(aes(month.CY,mean.TP))+
  theme_bw()


library(gganimate)
#library(tidyverse)
#library(lubridate)

#frame_count=nrow(grab.wq.month)#diff(range(grab.wq.month$month.CY))/duration(1,"months")
#cycle_length=12

#sample3=map_df(seq_len(frame_count),~grab.wq.month,.id="id")%>%
#  mutate(id=as.integer(id))%>%
#  mutate(view_date=as.Date(min(as.Date(grab.wq.month$month.CY))+duration(id-1,"month")))%>%
#  filter(month.CY <= view_date)%>%
#  mutate(months_ago=month(view_date)-month(month.CY),
#         phase_dif=(months_ago %% cycle_length)/cycle_length,
#         x_pos=-sin(2*pi*phase_dif),
#         nearness=cos(2*pi*phase_dif))

#cols=wes_palette("Zissou1",nrow(grab.wq.month),"continuous")

#b=ggplot(sample3)+
#  geom_path(aes(x_pos,mean.TP,alpha=nearness,color=months_ago,size=-months_ago*0.5))+
#  #ylim(0,40)+
#  scale_size(range=c(0,5))+
#  transition_manual(id)+theme_void()+
#  guides(size="none",alpha="none",color="none")+
#  ggtitle("Monthly SRS1 TP concentration")

#b
#animate(b,fps=25,duration=15,width=300,height=150)


## None tidy version
frame.count=nrow(grab.wq.month)
cycle.length=12

sample1=data.frame()
for(i in 1:frame.count){
  tmp=grab.wq.month
  tmp$id=as.integer(i)
  sample1=rbind(sample1,tmp)
}
sample1$view.date=min(sample1$month.CY)+months(sample1$id-1)
sample1=subset(sample1,month.CY<=view.date)

sample1$months.ago=with(sample1,as.numeric(round((view.date-month.CY)*3.80517e-7)))
sample1$phase.dif=with(sample1,(months.ago%%cycle.length)/cycle.length)
sample1$x.pos=with(sample1,-sin(2*pi*phase.dif))
sample1$nearness=with(sample1,cos(2*pi*phase.dif))

#https://github.com/thomasp85/gganimate/wiki/Temperature-time-series
#https://twitter.com/JustTheSpring/status/1049468485334523906?s=20

#library(wesanderson)
#cols=wes_palette("Zissou1",frame.count,"continuous")
#cols.df=data.frame(colors.val=as.character(cols),id=1:frame.count)
#sample1=merge(sample1,cols.df,"id")

b=ggplot(sample1)+
  geom_path(aes(x.pos,mean.TP,alpha=nearness,color=months.ago,size=-months.ago),lineend="round",show.legend = F)+
  #geom_segment(aes(xend=0.95,yend=mean.TP),linetype=2,colour="grey")+
  #geom_point(aes(x.pos,mean.TP),size=months.ago/12)+
  #geom_text(aes(x=1.25,y=mean.TP,label=month.CY),hjust=0)+
  geom_label(aes(x=0.9,y=100,label=format(month.CY,"%m-%Y")),label.r=unit(0,"lines"),label.size=1)+
  scale_y_continuous("Total P (\u03BCg L\u207B\u00B9)",breaks=seq(0,100,20),limits=c(0,100))+
  scale_x_continuous(breaks = NULL,limits=c(-1,1))+
  theme_minimal(base_size=15,base_line_size = 1,base_rect_size = 1)+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.5,0.25,0.5,1), "cm"),
        text=element_text(family="serif"))+
  transition_manual(view.date)+
  ggtitle("Monthly Surface Water TP concentration at SRS6")

animate(b,fps=15)

anim_save(paste0(plot.path,"SRS6_TP.gif"))


grab.wq.month=subset(grab.wq.month,is.na(mean.TP)==T)
b2=ggplot(grab.wq.month,aes(month.CY,mean.TP))+
  geom_line(colour="dodgerblue1",linetype=1,size=1)+
  geom_segment(aes(xend=date.fun("2020-12-01"),yend=mean.TP),linetype=2,colour="grey")+
  geom_point(shape=19,colour="indianred",size=4)+
  geom_text(aes(x=date.fun("2020-12-01"),y=mean.TP,label=format(month.CY,"%m-%Y")),hjust=0)+
  transition_reveal(month.CY,keep_last=T)+
  scale_y_continuous("Total P (\u03BCg L\u207B\u00B9)",breaks=seq(0,100,20),limits=c(0,100))+
  scale_x_datetime("Date", limits=date.fun(c("2001-05-01","2021-12-01")))+
  theme_minimal(base_size=15,base_line_size = 1,base_rect_size = 1)+
  coord_cartesian(clip = 'off') + 
  theme(plot.margin = unit(c(0.5,1,0.5,1), "cm"),
        text=element_text(family="serif"))+
  ggtitle("Monthly Surface Water TP concentration at SRS6")

animate(b2,fps=3)
anim_save(paste0(plot.path,"SRS6_TP_timeseries.gif"))
