# 

# MARSHBIRD MODEL AVERAGING

# Paul J. Taillie
# North Carolina State University
# 12/12/17
  #Re-written 09/13/2018

# load packages
library(tidyverse)
library(unmarked)

# #clear environment
remove(list=ls())

# load packages
library(tidyverse)
library(unmarked)
library(AICcmodavg)
library(MuMIn)
library(gridExtra)

#standarize helper function
standardize<-function(xx){  
  xx.mean=mean(xx,na.rm=T)
  xx.sd=sd(xx,na.rm=T)
  xx.standardized<-(xx-xx.mean)/xx.sd
  return(xx.standardized)
}

# Read in data
BLRA<-read.csv("data/BLRA.csv")
CLRA<-read.csv("data/CLRA.csv")
LEBI<-read.csv("data/LEBI.csv")
VIRA<-read.csv("data/VIRA.csv")
SESP<-read.csv("data/SESP.csv")
covs<-read.csv("data/covs.csv")
date<-read.csv("data/date.csv")
time<-read.csv("data/time.csv")
type<-read.csv("data/type.csv")

# Standardize Observation Covariates and store in list
obs.covs=list(
  time=standardize(as.matrix(time[,2:9])),
  date=standardize(as.matrix(date[,2:9])),
  #type<-as.matrix(type.in%>%select(V1:V8)%>%mutate_all(funs(factor)))
  type=as.matrix(type[,2:9])
)
SPEC=SESP
# -------------- function to use model averaged 
# --------------coefficients to predict response (occupancy)  -----------
plot.occ<-function(SPEC,obs.covs,covs){
  # 
  if(is.na(SPEC[1,3])==T){
    index.2017<-which(covs$Year==2017)
    SPEC<-SPEC%>%
      separate(X,c("site","year"))%>%
      filter(year==2017)
    covs<-covs%>%
      filter(Year==2017) 
    obs.covs[[1]]<-obs.covs[[1]][index.2017,]
    obs.covs[[2]]<-obs.covs[[2]][index.2017,]
    obs.covs[[3]]<-obs.covs[[3]][index.2017,]
  }
    
  # Standardize Site Covariates
  siteCovs<- covs%>%
    select(Herb_hits:Shrub.Woody,Dist_to_For,num.fires)%>%  #select covariate cols
    mutate_all(funs(standardize))
  #siteCovs$year=factor(covs$Year) #add year back on as factor
  

    
  # set up unmarked frame
  umf <- unmarkedFrameOccu(y=SPEC[,4:11],
                           obsCovs=obs.covs,
                           siteCovs=siteCovs)
  
  #fit global model and dredge
  global.mod<-occu(~type+time+date+I(time^2)+I(date^2)~Herb_hits+
                     Woody_height+
                     Cladium+
                     Juncus+
                     #Year+
                     Dist_to_For+
                     num.fires,umf)
  dredge.out<-dredge(global.mod,rank=AIC)
  #get model averaged coefficients
  top.mods<-get.models(dredge.out,subset=delta<2)
  avgm<-model.avg(top.mods)
  results.out<-as.data.frame(avgm$coefArray[1,1:2,])
  #set up prediction data
  newdata.blank<-data.frame(
    Herb_hits=rep(0,20),
    Woody_height=rep(0,20),  
    Dist_to_For=rep(0,20),
    Cladium=rep(0,20),
    Juncus=rep(0,20),
    num.fires=rep(0,20),
    year=factor(rep("2017",20),levels=c("2016","2017")))
  
  #woody
  newdata.woody<-newdata.blank
  Woody_height.X<-seq(min(covs$Woody_height),max(covs$Woody_height),length.out=20)
  newdata.woody$Woody_height<-standardize(Woody_height.X)  
  woody.pred<-as.data.frame(predict(avgm,type="state",newdata=newdata.woody,se.fit=TRUE))
  woody.pred$Woody_height=Woody_height.X
  
  woody.plot<-ggplot(woody.pred)+
    geom_ribbon(aes(x=Woody_height,ymin=fit-se.fit,ymax=fit+se.fit),fill="grey69")+
    geom_line(aes(x=Woody_height,y=fit),size=1.3)+
    ylim(0,1)+
    theme(axis.ticks=element_blank(),axis.text.x=element_blank(),legend.position="none")+
    ylab("Occupancy Probability")+
    xlab("Woody vegetation height")
  
  #Distance to forest
  newdata.Dist_to_For<-newdata.blank
  Dist_to_For.X<-seq(min(covs$Dist_to_For),max(covs$Dist_to_For),length.out=20)
  newdata.Dist_to_For$Dist_to_For<-standardize(Dist_to_For.X)  
  Dist_to_For.pred<-as.data.frame(predict(avgm,type="state",newdata=newdata.Dist_to_For,se.fit=TRUE))
  Dist_to_For.pred$Dist_to_For=Dist_to_For.X
  
  Dist_to_For.plot<-ggplot(Dist_to_For.pred)+
    geom_ribbon(aes(x=Dist_to_For,ymin=fit-se.fit,ymax=fit+se.fit),fill="grey69")+
    geom_line(aes(x=Dist_to_For,y=fit),size=1.3)+
    ylim(0,1)+
    theme(axis.ticks=element_blank(),axis.text.x=element_blank(),legend.position="none")+
    ylab("Occupancy Probability")+
    xlab("Distance to Forest (m)")
  
  #Herb hits
  newdata.Herb<-newdata.blank
  Herb_hits.X<-seq(min(covs$Herb_hits),max(covs$Herb_hits),length.out=20)
  newdata.Herb$Herb_hits<-standardize(Herb_hits.X)  
  Herb.pred<-as.data.frame(predict(avgm,type="state",newdata=newdata.Herb,se.fit=TRUE))
  Herb.pred$Herb_hits=Herb_hits.X
  
  Herb.plot<-ggplot(Herb.pred)+
    geom_ribbon(aes(x=Herb_hits,ymin=fit-se.fit,ymax=fit+se.fit),fill="grey69")+
    geom_line(aes(x=Herb_hits,y=fit),size=1.3)+
    ylim(0,1)+
    theme(axis.ticks=element_blank(),axis.text.x=element_blank(),legend.position="none")+
    ylab("Occupancy Probability")+
    xlab("Herbaceous Vegetation (hits)")
  #Cladium
  newdata.cladium<-newdata.blank
  Cladium.X<-seq(min(covs$Cladium),max(covs$Cladium),length.out=20)
  newdata.cladium$Cladium<-standardize(Cladium.X)  
  cladium.pred<-as.data.frame(predict(avgm,type="state",newdata=newdata.cladium,se.fit=TRUE))
  cladium.pred$Cladium=Cladium.X
  
  cladium.plot<-ggplot(cladium.pred)+
    geom_ribbon(aes(x=Cladium,ymin=fit-se.fit,ymax=fit+se.fit),fill="grey69")+
    geom_line(aes(x=Cladium,y=fit),size=1.3)+
    ylim(0,1)+
    theme(axis.ticks=element_blank(),legend.position="none")+
    ylab("Occupancy Probability")+
    xlab("Cladium Cover")
  #Juncus
  newdata.juncus<-newdata.blank
  Juncus.X<-seq(min(covs$Juncus),max(covs$Juncus),length.out=20)
  newdata.juncus$Juncus<-standardize(Juncus.X)  
  juncus.pred<-as.data.frame(predict(avgm,type="state",newdata=newdata.juncus,se.fit=TRUE))
  juncus.pred$Juncus=Juncus.X
  
  juncus.plot<-ggplot(juncus.pred)+
    geom_ribbon(aes(x=Juncus,ymin=fit-se.fit,ymax=fit+se.fit),fill="grey69")+
    geom_line(aes(x=Juncus,y=fit),size=1.3)+
    ylim(0,1)+
    theme(axis.ticks=element_blank(),legend.position="none")+
    ylab("Occupancy Probability")+
    xlab("Juncus Cover")
  
  #Fire
  newdata.num.fires<-newdata.blank[1:5,]
  num.fires.X<-seq(1,5)
  newdata.num.fires$num.fires<-standardize(num.fires.X)  
  num.fires.pred<-as.data.frame(predict(avgm,type="state",newdata=newdata.num.fires,se.fit=TRUE))
  num.fires.pred$num.fires=num.fires.X
  
  num.fires.plot<-ggplot(num.fires.pred)+
    geom_ribbon(aes(x=num.fires,ymin=fit-se.fit,ymax=fit+se.fit),fill="grey69")+
    geom_line(aes(x=num.fires,y=fit),size=1.3)+
    ylim(0,1)+
    theme(axis.ticks=element_blank(),legend.position="none")+
    ylab("Occupancy Probability")+
    xlab("Number of Fires")
  
  #Combined Gradient plot
  newdata.combo<-newdata.blank
  newdata.combo$Juncus      <-rev(standardize(Juncus.X))  
  newdata.combo$Dist_to_For <-rev(standardize(Dist_to_For.X))
  newdata.combo$Woody_height<-standardize(Woody_height.X) 
  newdata.combo$Cladium    <-standardize(Cladium.X)
  combo.pred<-as.data.frame(predict(avgm,type="state",newdata=newdata.combo,se.fit=TRUE))
  combo.pred$gradient=seq(1:20)
  
  
  
  
  return(list(woody.plot,
              cladium.plot,
              juncus.plot,
              Herb.plot,
              Dist_to_For.plot,
              num.fires.plot,
              combo.pred,
              results.out))
}############  --------  end function   -------------------------------

#run function on each species
SESP.results<-plot.occ(SESP,obs.covs,covs)
CLRA.results<-plot.occ(CLRA,obs.covs,covs)
LEBI.results<-plot.occ(LEBI,obs.covs,covs)
VIRA.results<-plot.occ(VIRA,obs.covs,covs)





