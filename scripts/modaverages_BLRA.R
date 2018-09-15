# 

# MARSHBIRD MODEL AVERAGING

# Paul J. Taillie
# North Carolina State University
# 12/12/17
#Re-written 09/13/2018


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

# Read in black rail observations
blra<-read.csv(file="data/blra_ARU.csv")

# Read in observation covariates, standardize, and put in list
date<-read.csv(file="data/date_aru.csv",row.names=1)
time<-read.csv(file="data/time_aru.csv",row.names=1)
type<-read.csv(file="data/type_aru.csv",row.names=1)

obs.covs=list(
  time=standardize(as.matrix(time[,2:49])),
  date=standardize(as.matrix(date[,2:49])),
  type=as.matrix(type[,2:49])
)

#Standardize Site Covariates
covs<-read.csv(file="data/covs.csv")
siteCovs<- covs%>%
  select(Herb_hits:Shrub.Woody,Dist_to_For,num.fires)%>%  #select covariate columns
  mutate_all(funs(standardize))
siteCovs$year=factor(covs$Year)

umf <- unmarkedFrameOccu(y=blra[,3:50],
                                  obsCovs=obs.covs,
                                  siteCovs=siteCovs)
#fit global model and dredge
global.mod<-occu(~type+time+date+I(time^2)+I(date^2)~Herb_hits+
                   Woody_height+
                   Cladium+
                   Juncus+
                   year+
                   Dist_to_For+
                   num.fires,umf)

dist.mod<-occu(~type+time+date+I(time^2)+I(date^2)~Dist_to_For+
                 year,umf)
clad.mod<-occu(~type+time+date+I(time^2)+I(date^2)~Cladium+
                 year,umf)              
               
               
               
               
               