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

load(file="Create_data//simulated.worms.l.tidy.Rdata")
#go to simple/tidy_fit_N.R after

levels(res_tidy$inf_type)<-c("Arithmetic mean (+)","Geometric mean (+)","Geometric mean","Arithmetic mean","True inf. of pos.")
res_tidy$inf_type<-factor(res_tidy$inf_type,levels=c("Arithmetic mean (+)","Geometric mean (+)","Arithmetic mean","Geometric mean","True inf. of pos."))

levels(res_tidy$Scheme)<-c("1,1","1,2","2,1","2,2","3,2","3,3")
plt<-ggplot(res_tidy%>%filter(inf_type=="Geometric mean"),aes(x=`Obs. prev.`*100,y=`True prev.`*100,colour=`Inf. Int. (EPS)`))+geom_point(size=0.5)+facet_grid(`Scheme`~as.factor(N))
plt+scale_y_continuous(limits = c(0,100))+scale_color_gradient2(low = "#00FF00",mid="#FFFF00", high = "#FF0000",midpoint = 2,limits=c(0,6), space = "Lab",na.value = "grey50", guide = "colourbar")+labs(colour="Geometric mean \n (EPS)",y="True prevalence (%)",x="Observed prevalence (%)")
ggsave("Fig2.pdf",width=190,height=NA,units="mm",dpi=300)
ggsave("Fig2.png",width=190,height=NA,units="mm",dpi=300)

plt<-ggplot(res_tidy%>%filter(inf_type=="Geometric mean")%>%filter(N==50),aes(x=`Obs. prev.`*100,y=`True prev.`*100,colour=`Inf. Int. (EPS)`))+geom_point(size=0.5)+facet_grid(`Scheme`~.)
plt<-plt+scale_color_gradient2(low = "#00FF00",mid="#FFFF00", high = "#FF0000",midpoint = 2,limits=c(0,6), space = "Lab",na.value = "grey50", guide = "colourbar")+labs(colour="Geometric mean \n (EPS)",y="True prevalence(%)",x="Observed prevalence (%)")
plt+theme(legend.position = "bottom")
ggsave("obs-true2.pdf",width=100,height=NA,units="mm",dpi=300)
ggsave("obs-true2.png",width=100,height=NA,units="mm",dpi=300)

plt<-ggplot(res_tidy%>%filter(inf_type=="Geometric mean")%>%filter(`True prev.`<0.9)%>%filter(N==50),aes(x=`Obs. prev.`*100,y=`True prev.`*100,colour=`Inf. Int. (EPS)`))+geom_point(size=0.5)+facet_grid(`Scheme`~.)
plt+scale_color_gradient2(low = "#00FF00",mid="#FFFF00", high = "#FF0000",midpoint = 0.03,limits=c(0,0.06), space = "Lab",na.value = "grey50", guide = "colourbar")+labs(colour="Geometric mean \n (EPS)",y="True prevalence(%)",x="Observed prevalence (%)")


plt<-ggplot(res_tidy%>%filter(`Inf. Int. (EPS)`<50)%>%filter(`inf_type`%in% c("Geometric mean","Arithmetic mean"))%>%filter(N==50)%>%filter(`Obs. prev.`>0.05)%>%sample_frac(0.2),aes(x=`Inf. Int. (EPS)`,y=`Tot. KK sens. (%)`*100,colour=`Obs. prev.`*100))+geom_point(size=0.5)+facet_grid(`Scheme`~`inf_type`,scales = "free_x")
plt+scale_color_gradient2(low = "#00FF00",mid="#FFFF00", high = "#FF0000",midpoint = 50,limits=c(0,100), space = "Lab",na.value = "grey50", guide = "colourbar")+
  theme(legend.position = "bottom")+labs(colour="Observed prevalence (%)",y="Mean Kato-Katz sensitivity (%)",x="Mean infection intensity (EPS)")
ggsave("Fig1.pdf",width=100,height=NA,units="mm",dpi=300)
ggsave("Fig1.png",width=100,height=NA,units="mm",dpi=300)

plt<-ggplot(res_tidy%>%filter(`Inf. Int. (EPS)`<50)%>%filter(`inf_type`%in% c("Arithmetic mean (+)","Geometric mean (+)"))%>%filter(N==50)%>%filter(`Obs. prev.`>0.05)%>%sample_frac(0.2),aes(x=`Inf. Int. (EPS)`,y=`Tot. KK sens. (%)`*100,colour=`Obs. prev.`*100))+geom_point(size=0.5)+facet_grid(`Scheme`~`inf_type`,scales = "free_x")
plt+scale_color_gradient2(low = "#00FF00",mid="#FFFF00", high = "#FF0000",midpoint = 50,limits=c(0,100), space = "Lab",na.value = "grey50", guide = "colourbar")+
  theme(legend.position = "bottom")+labs(colour="Observed prevalence (%)",y="Mean Kato-Katz sensitivity (%)",x="Mean infection intensity (EPS)")
ggsave("S4_Fig.pdf",width=100,height=NA,units="mm",dpi=300)
ggsave("S4_Fig.png",width=100,height=NA,units="mm",dpi=300)


plt<-ggplot(res_tidy%>%filter(inf_type=="Arithmetic mean")%>%sample_frac(0.1),aes(x=`Inf. Int. (EPS)`,y=`Tot. KK sens. (%)`,colour=`Obs. prev.`))+geom_point()+facet_grid(Scheme~N,scales = "free")
plt+scale_color_gradient2(low = "#00FF00",mid="#FFFF00", high = "#FF0000",midpoint = 0.5,limits=c(0,1), space = "Lab",na.value = "grey50", guide = "colourbar")

ggsave("sim_overview_arith.pdf",width=190,height=NA,units="mm",dpi=300)
ggsave("sim_overview_arith.png",width=190,height=NA,units="mm",dpi=300)

plt<-ggplot(res_tidy%>%filter(inf_type=="Geometric mean")%>%sample_frac(0.1),aes(x=`Inf. Int. (EPS)`,y=`Tot. KK sens. (%)`,colour=`Obs. prev.`))+geom_point()+facet_grid(Scheme~N,scales = "free")
plt+scale_color_gradient2(low = "#00FF00",mid="#FFFF00", high = "#FF0000",midpoint = 0.5,limits=c(0,1), space = "Lab",na.value = "grey50", guide = "colourbar")

ggsave("sim_overview_geom.pdf",width=190,height=NA,units="mm",dpi=300)
ggsave("sim_overview_geom.png",width=190,height=NA,units="mm",dpi=300)


plt<-ggplot(res_tidy%>%filter(inf_type=="Geometric mean")%>%sample_frac(0.1),aes(x=`Obs. prev.`,y=`Tot. KK sens. (%)`,colour=`Mean worms`))+geom_point()+facet_grid(Scheme~N,scales = "free")
plt+scale_color_gradient2(low = "#00FF00",mid="#FFFF00", high = "#FF0000",midpoint = 50,limits=c(0,100), space = "Lab",na.value = "grey50", guide = "colourbar")

ggsave("sim_overview_prev.pdf",width=190,height=NA,units="mm",dpi=300)
ggsave("sim_overview_prev.png",width=190,height=NA,units="mm",dpi=300)


plt<-ggplot(res_tidy%>%filter(inf_type=="True inf. of pos.")%>%sample_frac(0.1),aes(y=`Obs. prev.`,x=`Inf. Int. (EPS)`,colour=`Tot. KK sens. (%)`))+geom_point()+facet_grid(Scheme~N,scales = "free")
plt+scale_color_gradient2(low = "#00FF00",mid="#FFFF00", high = "#FF0000",midpoint = 0.5,limits=c(0,1), space = "Lab",na.value = "grey50", guide = "colourbar")

#got to tidy_fit
ggplot(res_tidy)+geom_point(aes(x=`Mean worms`,y=`Tot. KK sens. (%)`))+facet_grid( `Scheme` ~. ,scales = "free")
ggplot(res_tidy)+geom_point(aes(x=`Mean worms`,y=`True prev.`))



ggplot()+geom_point(data=res%>%filter(`Geom. mean of tot.`<10),aes(x=`True prev.`,y=`Tot. KK sens. (%)`),size=0)
ggplot()+geom_point(data=res_tidy%>%spread(key=inf_type,value=`Inf. Int. (EPS)`)%>%filter(`Geom. mean of tot.`<10),aes(x=`Obs. prev.`,y=`Tot. KK sens. (%)`),size=0)+facet_grid(Scheme~.)
temp<-res_tidy%>%filter(inf_type=="Geometric mean")%>%select(pop_var)
hist(temp$pop_var)


te<-res_tidy%>%filter(inf_type=="Geometric mean")%>%filter(`Obs. prev.`<0.05)%>%filter(N==5000)
te%>%group_by(Scheme)%>%summarise(count=n())
