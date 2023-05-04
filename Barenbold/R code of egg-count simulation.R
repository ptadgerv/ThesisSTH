
gm_mean0 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x[x > 0]))
}
gm_mean1 = function(x, na.rm=TRUE){
  exp(sum(log(x+1), na.rm=na.rm) / length(x))
}

library(tidyverse)
# load(file="trans5.RData")
#select 500 steps from prior samples
numb.p<-20
# priors<-as.data.frame(priors)
# priors.short<-priors[round(seq(from=1,to=nrow(priors),len=numb.p)),]
# w_agg_0<-rep(-2.89,numb.p)#priors.short$w_agg_0
w_egg<-rnorm(numb.p,0.2,0.05)
d_var_1<-rnorm(numb.p,1.17,0.08)
count_var_1<-rnorm(numb.p,6.23,1.07)#taken from priors
# count_var_1<-priors.short$count_var_1
min_infect<-0.05
#mu unnormalized (factor 2 is roughly assumed so we go 5 to 300)
numb.w<-30
mean_worms<-exp(seq(from=log(10),to=log(400),length=numb.w))
#sample sizes testes
N.sample<-c(30,50,5000)##50 or 5000
#number of repeats we will do which matches also the number of priors
numb.r<-5
Ntot<-numb.w*numb.r*numb.p*length(N.sample)
res<-data.frame(inf_obs11=rep(NA,Ntot),inf_obs12=rep(NA,Ntot),inf_obs21=rep(NA,Ntot),inf_obs22=rep(NA,Ntot),inf_obs32=rep(NA,Ntot),inf_obs33=rep(NA,Ntot),inf=rep(NA,Ntot)
                ,prev_obs11=rep(NA,Ntot),prev_obs12=rep(NA,Ntot),prev_obs21=rep(NA,Ntot),prev_obs22=rep(NA,Ntot),prev_obs32=rep(NA,Ntot),prev_obs33=rep(NA,Ntot)
                ,s_kk11=rep(NA,Ntot),s_kk12=rep(NA,Ntot),s_kk21=rep(NA,Ntot),s_kk22=rep(NA,Ntot),s_kk32=rep(NA,Ntot),s_kk33=rep(NA,Ntot),N=rep(NA,Ntot)
                ,geom_obs_11=rep(NA,Ntot),geom_obs_12=rep(NA,Ntot),geom_obs_21=rep(NA,Ntot),geom_obs_22=rep(NA,Ntot),geom_obs_32=rep(NA,Ntot),geom_obs_33=rep(NA,Ntot),
                w_agg=rep(NA,Ntot),pop_var=rep(NA,Ntot),prev=rep(NA,Ntot),mean_worms=rep(NA,Ntot))
z<-1
for(i in 1:length(N.sample)){#
  for(j in 1:numb.w){
    for(l in 1:numb.p){
     for(k in 1:numb.r){
      res$mean_worms[z]<-mean_worms[j]
      res$w_agg[z]=0.2+mean_worms[j]/600*0+rnorm(1,mean=0,sd=0.03)
      q_w<-0.3
      q_m<-1-q_w
      res$prev[z]=1-(1+q_w*mean_worms[j]/res$w_agg[z])^(-res$w_agg[z])-(1+q_m*mean_worms[j]/res$w_agg[z])^(-res$w_agg[z])+(1+mean_worms[j]/res$w_agg[z])^(-res$w_agg[z])
      temp<-rnbinom(N.sample[i],mu=res$mean_worms[z],size=res$w_agg[z])
      pairs<-temp
      for(m in 1:N.sample[i]){
        t<-rbinom(1,size=temp[m],prob=q_w)
        pairs[m]<-min(t,temp[m]-t)
      }
      temp<-pairs*0.2
      infected<-rep(0,N.sample[i])
      infected[which(temp>0)]<-1
      res$inf[z]<-mean(temp)
      var_inf<-var(temp)
      res$pop_var[z]=res$inf[z]/var_inf
      res$N[z]<-N.sample[i]
      Y<-array(data=NA,dim=c(N.sample[i],3,3))
      P<-array(data=NA,dim=c(N.sample[i],3,3))
      tempday<-array(data=NA,dim=c(N.sample[i],3))

      for(d in 0:(N.sample[i]-1)){
      tempsd<-d_var_1[l]
      tempmean<--tempsd^2/2
      tempday[d+1,]<-temp[d+1]*exp(rnorm(3,mean=tempmean,sd=tempsd))
      Y[d+1,1,]<-rnbinom(3,mu=tempday[d+1,1],size=count_var_1[l]+(1-infected[d+1])+infected[d+1]*0.001)
      Y[d+1,2,]<-rnbinom(3,mu=tempday[d+1,2],size=count_var_1[l]+(1-infected[d+1])+infected[d+1]*0.001)
      Y[d+1,3,]<-rnbinom(3,mu=tempday[d+1,3],size=count_var_1[l]+(1-infected[d+1])+infected[d+1]*0.001)
      P[d+1,1,]<-as.numeric(as.logical(Y[d+1,1,]))
      P[d+1,2,]<-as.numeric(as.logical(Y[d+1,2,]))
      P[d+1,3,]<-as.numeric(as.logical(Y[d+1,3,]))
      }
      res$inf_obs11[z]<-mean(Y[,1,1])
      res$prev_obs11[z]<-mean(P[,1,1])

      res$inf_obs12[z]<-mean(Y[,1,1:2])
      res$prev_obs12[z]<-mean(apply(P[,1,1:2],1,max))

      res$inf_obs21[z]<-mean(Y[,1:2,1])
      res$prev_obs21[z]<-mean(apply(P[,1:2,1],1,max))

      res$inf_obs22[z]<-mean(Y[,1:2,1:2])
      res$prev_obs22[z]<-mean(apply(P[,1:2,1:2],1,max))

      res$inf_obs32[z]<-mean(Y[,1:3,1:2])
      res$prev_obs32[z]<-mean(apply(P[,1:3,1:2],1,max))

      res$inf_obs33[z]<-mean(Y[,1:3,1:3])
      res$prev_obs33[z]<-mean(apply(P[,1:3,1:3],1,max))

      # res$geom_obs_11[z]<-gm_mean0(Y[,1,1])
      # res$geom_obs_12[z]<-gm_mean0(apply(Y[,1,1:2],1,mean))
      # res$geom_obs_21[z]<-gm_mean0(apply(Y[,1:2,1],1,mean))
      # res$geom_obs_22[z]<-gm_mean0(apply(Y[,1:2,1:2],1,mean))
      # res$geom_obs_32[z]<-gm_mean0(apply(Y[,1:3,1:2],1,mean))
      # res$geom_obs_33[z]<-gm_mean0(apply(Y[,1:3,1:3],1,mean))
      #meaning I take the geometric mean of all results adding 1
      res$geom_obs_11[z]<-gm_mean1(Y[,1,1])
      res$geom_obs_12[z]<-gm_mean1(as.vector(Y[,1,1:2]))
      res$geom_obs_21[z]<-gm_mean1(as.vector(Y[,1:2,1]))
      res$geom_obs_22[z]<-gm_mean1(as.vector(Y[,1:2,1:2]))
      res$geom_obs_32[z]<-gm_mean1(as.vector(Y[,1:3,1:2]))
      res$geom_obs_33[z]<-gm_mean1(as.vector(Y[,1:3,1:3]))

      res$s_kk11[z]<-mean(1-(count_var_1[l]/(tempday[which(as.logical(infected)),1]+count_var_1[l]))^(count_var_1[l]))
      res$s_kk12[z]<-mean(1-(count_var_1[l]/(tempday[which(as.logical(infected)),1]+count_var_1[l]))^(2*count_var_1[l]))
      res$s_kk21[z]<-1-mean((count_var_1[l]/(tempday[which(as.logical(infected)),1]+count_var_1[l]))^(count_var_1[l])*(count_var_1[l]/(tempday[which(as.logical(infected)),2]+count_var_1[l]))^(count_var_1[l]))
      res$s_kk22[z]<-1-mean((count_var_1[l]/(tempday[which(as.logical(infected)),1]+count_var_1[l]))^(2*count_var_1[l])*(count_var_1[l]/(tempday[which(as.logical(infected)),2]+count_var_1[l]))^(2*count_var_1[l]))
      res$s_kk32[z]<-1-mean((count_var_1[l]/(tempday[which(as.logical(infected)),1]+count_var_1[l]))^(2*count_var_1[l])*(count_var_1[l]/(tempday[which(as.logical(infected)),2]+count_var_1[l]))^(2*count_var_1[l])*(count_var_1[l]/(tempday[which(as.logical(infected)),3]+count_var_1[l]))^(2*count_var_1[l]))
      res$s_kk33[z]<-1-mean((count_var_1[l]/(tempday[which(as.logical(infected)),1]+count_var_1[l]))^(3*count_var_1[l])*(count_var_1[l]/(tempday[which(as.logical(infected)),2]+count_var_1[l]))^(3*count_var_1[l])*(count_var_1[l]/(tempday[which(as.logical(infected)),3]+count_var_1[l]))^(3*count_var_1[l]))
      z<-z+1
      }
    }
  }
}


hist(res$pop_var)
plot(res$pop_var,res$prev)
plot(res$inf,res$pop_var)
save(res,file="Create_data/simulated.worms.l.Rdata")
#go to tidy res new
ggplot(res)+geom_point(aes(x=mean_worms,y=w_agg))
hist(res$prev,100)
hist(res$prev_obs11,100)
hist(res$prev_obs12,100)
