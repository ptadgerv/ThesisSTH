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
getwd()
load(file="Create_data//simulated.worms.l.tidy.Rdata")

scheme<-data.frame(s=c(1,2,1,2,2,3),d=c(1,1,2,2,3,3))
Numbs<-c(30,50)
fitanddat<-list()
for(i in 1:6){
  for(j in 1:2){
  s<-scheme$s[i]
  d<-scheme$d[i]
  Numb<-Numbs[j]

#%>%filter(N==Numb)
res<-res_tidy%>%filter(N%in%Numb)%>%filter(Samples==s)%>%filter(Days==d)%>%filter(inf_type %in% c("Geom. mean of tot.","Tot. arith. mean inf.","True inf. of pos."))%>%spread(key=inf_type,value=`Inf. Int. (EPS)`)#%>%sample_n(1000)
ggplot()+geom_point(data=res%>%filter(`Geom. mean of tot.`<10),aes(x=`Geom. mean of tot.`,y=`Tot. KK sens. (%)`),size=0)+facet_grid(N~.)



stan_data<-list(N=res$N,obs_prev=res$`Obs. prev.`,true_prev=res$`True prev.`,N_ind=n_distinct(res$Ref),inf_eps=res$`Tot. arith. mean inf.`,inf_eps_g=res$`Geom. mean of tot.`,s_kk=res$`Tot. KK sens. (%)`,prev_obs=res$`Obs. prev.`)
#tidy full model includes variation with N, and some variation with infection intensity
#newnu includes more complex variation with infection intensity that first increases and then decreases but only for geometric mean so far
fit_obs<-stan(file='Create_data/model_exp_fin_obs.stan',pars=c("a0","a1","b0","b1"),data=stan_data,warmup=500,iter=2000,chains=4,thin=1,control = list(adapt_delta = 0.95,max_treedepth=16))#stepsize=0.04?
fit_g<-stan(file='Create_data/model_exp_fin.stan',pars=c("a0","a1","a2","b0","b1"),data=stan_data,warmup=500,iter=2000,chains=4,thin=1,control = list(adapt_delta = 0.95,max_treedepth=16))#stepsize=0.04?
fit_a<-stan(file='Create_data/model_exp_fin_a.stan',pars=c("a0","a1","a2","b0","b1"),data=stan_data,warmup=500,iter=2000,chains=4,thin=1,control = list(adapt_delta = 0.95,max_treedepth=16))#stepsize=0.04?

fitanddat[[2*(i-1)+j]]=list(fit_g=fit_g,fit_a=fit_a,fit_obs=fit_obs,res=res,stan_data=stan_data,s=s,d=d,Numb=Numb)
}
}
save(fitanddat,file=paste(paste("E:/TPH/Saved_simulations/pop_sen/fitanddat_geom_",format(Sys.time(), "%Y-%m-%d_%H-%M"),sep="_"), "RData", sep = "."))


# analyze_short<-as.shinystan(fit,pars=c("a0","a1","a2","a3","b0","b1","b2","b3","a0_g","a1_g","a2_g","a3_g","b0_g","b1_g","b2_g","b3_g"))
# analyze<-as.shinystan(fit,pars=c("a0","a1","a2","a3","b0","b1","a0_g","a1_g","a2_g","a3_g","b0_g","b1_g","b2_g","alpha","beta","mu","nu","alpha_g","beta_g","mu_g","nu_g"))
analyze_short<-as.shinystan(fit_obs,pars=c("a0","a1","b0","b1"))
# analyze_short<-as.shinystan(fit_12,pars=c("a0","a1","a2","a3","b0","b1","b2","b3","b4","b5","a5"))


fit_shiny<-launch_shinystan(analyze_short)


setwd("C:\\Users\\32498\\Documents\\Master EPI\\SEM6\\Thesis\\Create_data\\")
list.files()

"S2 Stan code of egg-count model.stan"
"S3 Stan code of statistical model fit to geometric mean simulation.stan"
"S4 Stan code of statistical model fit to arithmetic mean simulation.stan"
"S5 Stan code of statistical model fit to observed prevalence simulation.stan"
"S6 Stan code of egg-count model used for validation including informative priors.stan"

"simulated.worms.l.Rdata"
"simulated.worms.l.tidy.Rdata"
