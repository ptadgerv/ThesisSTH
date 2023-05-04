rm(list=ls())
library(rstan)
library(parallel)
library(loo)
library(tidyverse)
#library(rstanmulticore) #already implemented in rstan by now
library(coda)
library(shinystan)
library(boot)
require(Rcpp)
library(foreign)
library(reshape2)
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())
options(shinystan.rstudio = TRUE)


load(file="Create_data/simulated.worms.l.Rdata")
getwd()
res<-as.tibble(res)
res$ref<-1:nrow(res)
res_tidy<-res%>%drop_na()#drop observations with NA in them
res_tidy<-res_tidy%>%gather(inf_obs11:inf_obs33,key="scheme",value="inf_obs")%>%separate(scheme,into=c("trash","scheme"),sep=-2)%>%select(-trash)%>%separate(scheme,into=c("days","samples"),sep=1,convert=T)
res_tidy<-res_tidy%>%gather(prev_obs11:prev_obs33,key="scheme",value="prev_obs")%>%separate(scheme,into=c("trash","scheme"),sep=-2)%>%select(-trash)%>%separate(scheme,into=c("days1","samples1"),sep=1,convert=T)%>%filter(days==days1)%>%filter(samples==samples1)
res_tidy<-res_tidy%>%gather(s_kk11:s_kk33,key="scheme",value="s_kk")%>%separate(scheme,into=c("trash","scheme"),sep=-2)%>%select(-trash)%>%separate(scheme,into=c("days2","samples2"),sep=1,convert=T)%>%filter(days==days2)%>%filter(samples==samples2)
res_tidy<-res_tidy%>%gather(geom_obs_11:geom_obs_33,key="scheme",value="geom_obs")%>%separate(scheme,into=c("trash","scheme"),sep=-2)%>%select(-trash)%>%separate(scheme,into=c("days3","samples3"),sep=1,convert=T)%>%filter(days==days3)%>%filter(samples==samples3)
res_tidy<-res_tidy%>%select(-days2,-days3,-samples2,-samples3)
res_tidy<-res_tidy%>%unite(days1,samples1,col="scheme",sep="")
res_tidy<-res_tidy%>%mutate(inf_obs_pos=inf_obs*1.02/(prev_obs+0.01))
res_tidy<-res_tidy%>%mutate(geom_obs_pos=geom_obs^(1/prev_obs)-1)%>%mutate(geom_obs=geom_obs-1)#order of transformation matters and I'm not really sure why

names(res_tidy)<-c("True inf. of pos.","N","w_agg","pop_var","True prev.","Mean worms","Ref","Days","Samples","Tot. arith. mean inf.","Scheme","Obs. prev.","Tot. KK sens. (%)","Geom. mean of tot.","Arith. mean of pos.","Geom. mean of pos.")
res_tidy<-res_tidy%>%gather(`True inf. of pos.`,`Tot. arith. mean inf.`,`Arith. mean of pos.`,`Geom. mean of pos.`,`Geom. mean of tot.`,key="inf_type",value="Inf. Int. (EPS)")
res_tidy$Scheme<-as.factor(res_tidy$Scheme)
res_tidy$inf_type<-as.factor(res_tidy$inf_type)


save(res_tidy,file="Create_data/simulated.worms.l.tidy.Rdata")
#go to plot tidy