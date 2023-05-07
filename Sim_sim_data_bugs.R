#######################################################################################
################################# prevalence structure  ###############################
#################################for school_j & prevalence_i ##########################
#######################################################################################
rm(list = ls())
library(dplyr)
set.seed(2023)
gridd<-expand.grid(schools=c(10,20), prev.STH=c(0.01,0.03),var.STH=0.001,pooln1=10, CV.i=1.5,mu=1000, CV.d=0.75,   pooln2=25, CV.d.10=0.6,CV.d.25=0.5, CV.s=0.25, plotss=FALSE, TPR=0.88*0.88,TNR=0.97*0.97, n.s=100, n.days=c(1,2, 3), n.sample=c(1,2,3))
gridd$nnn=gridd$schools*gridd$n.s
gridd$Scenario<-rownames(gridd)
gridd$name<-paste0("Number of shools:",gridd$schools,". STH Prevalence:",gridd$prev.STH,". Number of days:",gridd$n.days,". Number of samples:",gridd$n.sample)
gridd$daysXsample<-gridd$n.days*gridd$n.sample
gridd<-gridd[order(gridd$daysXsample),]
gridd$scenario<-1:nrow(gridd)
gridd$Info_level<-ifelse(gridd$daysXsample==1,"Low",ifelse(gridd$daysXsample==2,"Low Medium",ifelse(gridd$daysXsample==3,"Medium",ifelse(gridd$daysXsample==4,"Medium High",ifelse(gridd$daysXsample==6,"High", ifelse(gridd$daysXsample==9,"Very High",NA))))))
library(openxlsx)
openxlsx::write.xlsx(gridd, "gridd.xlsx")
gridd<-openxlsx::read.xlsx("C:\\Users\\32498\\Downloads\\gridd.xlsx")
gridd<-gridd[gridd$Decission=="keep",]
grid.sim<-vector(mode='list', length=nrow(gridd))
for (ji in 1:nrow(gridd)) {
grid.sim[[ji]]<-gridd[ji,]
}
library(dplyr)
library(tidyverse)
library(Hmisc)
library(vcd)
library(VGAM)

#params=grid.sim


start_time <- Sys.time()
simSTH <- function(datatest, params = list()) {
datatest<-vector(mode='list', length=length(params))
#ji=1
for (ji in 1:length(params)) {

  nnn=params[[ji]][["nnn"]]
  schools=params[[ji]][["schools"]]
  pool=params[[ji]][["pooln1"]]
  prev.STH =params[[ji]][["prev.STH"]]
  var.STH=params[[ji]][["var.STH"]]


  n.s=nnn/schools  #subjects in each school
  S.label<-1:schools
  datatest[[ji]]<-rep(S.label,each=n.s)
  datatest[[ji]]<-as.data.frame(datatest[[ji]])
  colnames(datatest[[ji]])<-"School.name"
  datatest[[ji]]$ind.n<-NA
  datatest[[ji]]$ind.n<- rep(1:n.s,each=1, times=schools)

  n.pool<-1:(n.s/pool)
  datatest[[ji]]$pool10<-rep(n.pool,each=pool)
  #Using the idea of the nomogram of Fagan to get miun and max of the prevalence

  #  odds=prev.STH/(1-prev.STH)
  #  ratio=3  #Similar to likelihood ratio
  #  odd.p.i<-exp(log(odds)-log(ratio))
  #  odd.p.s<-exp(log(odds)+log(ratio))
  #temp1<-runif(n=schools, min = odd.p.i/(1+odd.p.i), max = odd.p.s/(1+odd.p.s))

  estBetaParams <- function(mu, var) {
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return(params = list(alpha = alpha, beta = beta))
  }
  #mu.neta ∈ [0,1]
  #var.beta ∈(0,0.5**2)
  #2%=2/100
  temp1<-rbeta(n=schools,
               shape1=(estBetaParams(mu=prev.STH,var=var.STH)$alpha),
               shape2 =(estBetaParams(mu=prev.STH,var=var.STH)$beta ) )   #Prevalence for each school
  datatest[[ji]]$prev<-rep(temp1,each=n.s)      #Prevalence for each school repeated for each subject
  temp2<-list()
  datatest[[ji]]$DX<-NA
  for (ii in 1:schools) {
    temp2[[ii]]<-rbinom(n=n.s,size=1,prob=temp1[ii]) # Number of binomial opsbervation= n.s # subj per school, 1 trial, prevalence of the school
    datatest[[ji]]$DX[ datatest[[ji]]$School.name==ii ]<-  temp2[[ii]]
  }
  datatest[[ji]]<-datatest[[ji]]%>%group_by(School.name,ind.n) %>%dplyr::mutate(ID.sub = cur_group_id())

  n.s=params[[ji]]$n.s
  n.days=params[[ji]]$n.days
  n.sample=params[[ji]]$n.sample
  nnn=params[[ji]]$nnn
  TPR=params[[ji]]$TPR
  TNR=params[[ji]]$TNR
  schools=params[[ji]]$schools
  ####Check-up: n.s should be a multiple of ss1 and ss2 ()
  ss1=n.s/params[[ji]]$pooln1 #Number of subsets for pool1
  ss2=n.s/params[[ji]]$pooln2 #Number of subsets for pool2
  if (!n.s%%ss1==0 & !n.s%%ss2==0 ) {
    warning("Number of subjects in dataset should be a multiple of the numbers of the pooled samples (pooln1 and pooln2)")
    stop()
  }
  if (!n.s%%ss1==0 ) {
    warning("Number of subjects in each school should be a multiple of the numbers of the pooled samples (pooln1)")
    stop()
  }
  if (!n.s%%ss2==0 ) {
    warning("Number of subjects in each shool  should be a multiple of the numbers of the pooled samples (pooln2)")
    stop()
  }
  ##  1 Variability in mean egg intensity between individuals, gamma dist, params$CV.i=1.5
  # I am selecting one Ki value, but not sure if I need to pick them all
  # The DGM can explored a grid with different values of K (or params$CV.i)
  #number of subjects
  Ki=1/params[[ji]]$CV.i^2
  datatest[[ji]]$Mu<- params[[ji]]$mu
  datatest[[ji]]$CV.d<- params[[ji]]$CV.d


  #replicating the datasets n.days times.
  datatest[[ji]]<- map(seq_len(n.days),~datatest[[ji]]) %>% bind_rows(.id="day")
  #datatest<-datatest.copy
  datatest[[ji]]<-datatest[[ji]][order(datatest[[ji]]$School.name,datatest[[ji]]$ind.n,datatest[[ji]]$day),]
  datatest[[ji]]$day<-as.numeric(datatest[[ji]]$day)

  #replicating the datasets n.sample times.
  datatest[[ji]]<- map(seq_len(n.sample),~datatest[[ji]]) %>% bind_rows(.id="sample")
  datatest[[ji]]<-datatest[[ji]][order(datatest[[ji]]$School.name,datatest[[ji]]$ind.n,datatest[[ji]]$day,datatest[[ji]]$sample),]
  datatest[[ji]]$sample<-as.numeric(datatest[[ji]]$sample)


  # Mid: for the days
  #nnn*n.days*n.sample

  datatest[[ji]]$NegDX<-ifelse(datatest[[ji]]$DX==1,0,1)
  prueb3<-datatest[[ji]]%>%group_by(School.name,day,sample) %>%  summarise(across(DX, sum))
  prueb4<-datatest[[ji]]%>%group_by(School.name,day,sample) %>%  summarise(across(NegDX, sum))
  colnames(prueb3)<-c("School.name","day","sample","Sick")
  colnames(prueb4)<-c("School.name","day","sample","Healthy")
  datatest[[ji]]<-merge(datatest[[ji]],prueb3,by = c("School.name","day","sample"))
  datatest[[ji]]<-merge(datatest[[ji]],prueb4,by = c("School.name","day","sample"))
  datatest[[ji]]$IDID<-1:nrow(datatest[[ji]])
  min.inf=0  #to produced a shifted gamma

  datatest[[ji]]<-datatest[[ji]][order(datatest[[ji]]$ID.sub,datatest[[ji]]$day,datatest[[ji]]$sample),]
  datatest[[ji]]$Mi<-datatest[[ji]]$Mid1<-datatest[[ji]]$Mids1<-NA
  for (ii in 1:(nnn)) {
    #o#o#       ii=1
    x.i = n.days * n.sample * (ii - 1) + 1 #BEGINING OF THE SUBJECT LINE
    x.f = x.i + n.days * n.sample - 1 #END OF THE SUBJECT LINE

    ###################### Shifted gamma using minimum infection: corresponding to the lowest possible infection with one worm pair. Barenbold (21) uses min_infect<-0.05.
    ###################### Mean=shape*scale, shape=Mean/scale
    datatest[[ji]]$Mi[datatest[[ji]]$ID.sub==ii] <- ifelse(
      mean(datatest[[ji]]$DX[datatest[[ji]]$ID.sub==ii]) == 1,
      rep(rgamma(1, shape = Ki, scale = (params[[ji]]$mu-min.inf)/Ki),
          each = 1,
          times = n.days * n.sample),
      rep(0, each = 1, times = n.days * n.sample))
    # datatest$Mi[x.i:x.f] <- rep(rgamma(1, shape = Ki, scale = (params$mu-min.inf)/Ki),each = 1,times = n.days * n.sample) Alternative to ZI-gamma > gamma
    for (dd in 1:n.days) {
      #o#o#             dd=2
      x.d = x.i + n.sample * (dd - 1) #-1
      #x.d.f = x.i + n.days * dd - 1
      x.d.f = x.d + n.sample - 1
      datatest[[ji]]$Mid1[x.d:x.d.f] = rep(
        rgamma(1,shape = params[[ji]]$CV.d ^ -2,
               scale = datatest[[ji]]$Mi[x.d:x.d.f] * params[[ji]]$CV.d ^ 2),
        each = 1,times = n.sample)

      for (ss in 1:n.sample) {
        #o#o# ss=1
        datatest[[ji]]$Mids1[x.d + ss - 1] = rgamma(1,shape = params[[ji]]$CV.d ^ -2,
                                                    scale = datatest[[ji]]$Mid1[x.d + ss - 1] * params[[ji]]$CV.d ^ 2)
      }
    }
  }

  hey<-datatest[[ji]]%>%group_by(School.name,day,sample, pool10) %>%dplyr::mutate(ID = cur_group_id())
  temp<-hey%>%group_by(ID)%>%summarise(Mids10=mean(Mids1))
  datatest[[ji]]<-merge(hey,temp, by="ID")

  # Mids: for the samples
  #It will be good to increase number of FN using a ZIP. Sen is now 95, it should be around 70%. DONE!
  #Specificity is perfect, then we should include FP

  datatest[[ji]]<-datatest[[ji]]%>%mutate(TP.ideal=round(Sick*TPR),TN.ideal=round(Healthy*TNR))
  datatest[[ji]]<-datatest[[ji]]%>%mutate(FN.ideal=Sick-TP.ideal,FP.ideal=Healthy-TN.ideal)
  datatest[[ji]]<-datatest[[ji]]%>%mutate(sens.ideal=TP.ideal/Sick,spec.ideal=TN.ideal/Healthy) #Only to check accuracy is ok
  datatest[[ji]]$DX.Miss<-datatest[[ji]]$DX

  for (SSS in 1:schools) {
    # Choose the ID of the missclasified subjects for each shcool
    subconj<-subset(datatest[[ji]],School.name==SSS)
    suj.0.ID<-subconj$IDID[subconj$DX==0] #All N (Negatives)
    suj.1.ID<-subconj$IDID[subconj$DX==1] #All P (Positives)
    ind.FP<-sample(suj.0.ID, subconj$FP.ideal[1], replace=FALSE) #Sample w/o replace FP-times
    ind.FN<-sample(suj.1.ID, subconj$FN.ideal[1], replace=FALSE) #Sample w/o replace FN-times
    datatest[[ji]]$DX.Miss[datatest[[ji]]$IDID %in% ind.FP]<-rep(1,each=1,times=length(ind.FP)) #FP:0 is replaced with 1
    datatest[[ji]]$DX.Miss[datatest[[ji]]$IDID %in% ind.FN]<-rep(0,each=1,times=length(ind.FN)) #FN:1 is replaced with 0
  }

  #Introduce FP: Introducing missclasification
  #Increase FN (to decrease sensitivity): using a ZIP
  datatest[[ji]]$Class2x2.a.ideal<-ifelse(((datatest[[ji]]$DX.Miss==datatest[[ji]]$DX)),"T","F")
  datatest[[ji]]$Class2x2.b.ideal<-ifelse(datatest[[ji]]$DX.Miss==1,"P","N")
  datatest[[ji]]$Class2x2.ideal<-paste0(datatest[[ji]]$Class2x2.a,datatest[[ji]]$Class2x2.b)

#  mids.count.error<-5
#  for (ii in 1:(nnn*n.days*n.sample)) {
#    datatest[[ji]]$count1[ii]<-ifelse(datatest[[ji]]$Class2x2.ideal[ii]=="FP",
#                                      rpois(1,  lambda=(datatest[[ji]]$Mids1[ii]+mids.count.error)),
#                                      rzipois(1,  lambda=datatest[[ji]]$Mids1[ii],pstr0 = 0.15))
#    datatest[[ji]]$count10[ii]<-rpois(1, lambda=datatest[[ji]]$Mids10[ii])
#  }


  mids.count.error<-5
  for (ii in 1:(nnn*n.days*n.sample)) {
    datatest[[ji]]$count1[ii]<-rpois(1,  lambda=(datatest[[ji]]$Mids1[ii]+mids.count.error))
      ifelse(datatest[[ji]]$Class2x2.ideal[ii]=="FP",
                                      rpois(1,  lambda=(datatest[[ji]]$Mids1[ii]+mids.count.error)),
                                      rpois(1,  lambda=(datatest[[ji]]$Mids1[ii])))
    datatest[[ji]]$count10[ii]<-rpois(1, lambda=datatest[[ji]]$Mids10[ii])
  }




  #data.table::fwrite(datatest[[ji]],"datatestdatatest.csv")
  datatest[[ji]]$TestR1<-ifelse(datatest[[ji]]$count1==0,0,1)
  datatest[[ji]]$TestR10<-ifelse(datatest[[ji]]$count10==0,0,1)

  datatest[[ji]]$Class2x2.a<-ifelse(((datatest[[ji]]$TestR1==datatest[[ji]]$DX)),"T","F")
  datatest[[ji]]$Class2x2.b<-ifelse(datatest[[ji]]$TestR1==1,"P","N")
  datatest[[ji]]$Class2x2<-paste0(datatest[[ji]]$Class2x2.a,datatest[[ji]]$Class2x2.b)
  #o#  ver<-data.frame(Class=datatest$Class2x2,TestR1=datatest$TestR1,DX=datatest$DX)
  #o#  tabla.d<-table(datatest$Class2x2,datatest$day)
  #o#  tabla.s<-table(datatest$Class2x2,datatest$sample)
  #o#  tabla.ds<-xtabs(~ Class2x2+day+sample, data=datatest)
  #o#  ftable(tabla.ds)%>%pander::pander(style = 'rmarkdown')
  #o#  colnames(tabla.d)<-c("day 1","day 3","day 3")
  #o#  colnames(tabla.s)<-c("Sample 1","Sample 3","Sample 3")
  #o#  colnames(tabla.ds)<-c("day 1","day 3","day 3")
  #o#  print.table(tabla.d)
  #o#  print.table(tabla.s)
  #o#  tabla.d<-rbind(tabla.d)
  #o#  tabla.s<-rbind(tabla.s)
  #o#  tabla.d.t<-as.data.frame(t(tabla.d))
  #o#  tabla.s.t<-as.data.frame(t(tabla.s))
  datatest[[ji]]<-subset(datatest[[ji]],select = -c(Class2x2.a,Class2x2.b,Class2x2, Class2x2.a.ideal,Class2x2.b.ideal))

  #o#   #sensibility
  #o#   print("sensibility day ")
  #o#   print(round((100*  tabla.d.t$TP/(  tabla.d.t$TP+  tabla.d.t$FN))[1:3],digits=2))
  #o#   #Specificity
  #o#   print("Specificity day ")
  #o#   print( round((100*  tabla.d.t$TN/(  tabla.d.t$TN+  tabla.d.t$FP)   )[1:3] ,digits=2)  )
  #o#   #sensibility
  #o#   print("sensibility sample ")
  #o#   print(round((100*tabla.s.t$TP/(  tabla.s.t$TP+  tabla.s.t$FN) )[1:3] ,digits=2))
  #o#   #Specificity
  #o#   print("Specificity sample ")
  #o#   print( round((100*tabla.s.t$TN/(  tabla.s.t$TN+  tabla.s.t$FP)[1:3]) ,digits=2))

  datatest[[ji]]<-select(datatest[[ji]],-prev)

  temp.1<-datatest[[ji]]%>%group_by(ID.sub)%>%summarise(count1.mean=round(mean(count1),digits=0))
  temp.10<-datatest[[ji]]%>%group_by(School.name,pool10)%>%summarise(count10.mean=round(mean(count10),digits=0))

  datatest[[ji]]<-merge(datatest[[ji]],temp.1, by="ID.sub")
  datatest[[ji]]<-merge(datatest[[ji]],temp.10, by=c("School.name","pool10"))


  var.labels = c(ID="Identifier for each School.name,day,sample and pool10",
                 ID.sub="Unique identifier for each subject",
                 sample="Sample number",
                 day="Day number",
                 Mu="Mean EPG by Country or Province",
                 School.name="School number",
                 ind.n="Subject ID per school (not unique)",
                 pool10="Pool ID",
                 DX="True Diagnostic status",
                 DX.Miss="Diagnostic status with missclasification",
                 Mi="Mean EPG in each subject",
                 Mids1="Mean EPG in each subject per day and sample",
                 Mids10="Mean EPG in each pool per day and sample",
                 count1="Count of EPG in thick smear by subject per day and sample",
                 count10="Count of EPG in thick smear by pool per day and sample",
                 TestR1="KK binary result by subject per day and sample",
                 TestR10="KK binary result by pool 10 per day and sample",
                 count1.mean="Average Count of EPG in thick smear by subject",
                 count10.mean="Average Count of EPG in thick smear by pool")

  label(datatest[[ji]]) <- as.list(var.labels[match(names(datatest[[ji]]), names(var.labels))])

  }

  return(datatest)
}


DF<-list()
reps=1:100
for (kk in 1:length(grid.sim)) {
#for (kk in 1:1) {
DF[[grid.sim[[kk]][["name"]] ]]<-lapply(
  X = reps,FUN = simSTH,
  params = list(grid.sim[[kk]]))
}
end_time <- Sys.time()
end_time - start_time

((end_time - start_time)*((24-9)/24))/(9/24)
#forecast time for 24 secanrios
end_time - start_time+((end_time - start_time)*((24-9)/24))/(9/24)

asas<-DF[["Number of shools:10. STH Prevalence:0.01. Number of days:1. Number of samples:1"]][[1]][[1]]

save.image(file = "final_100rep.RData")
#save.image(file = "SesionSimu.RData")
getwd()
load("Scenario1_100rep.RData")
load("C:\\Users\\32498\\Documents\\Master EPI\\SEM6\\Thesis\\SesionSimu.RData")

#rep1_scen1<-as.data.frame(DF[["Number of shools:20. STH Prevalence:0.01. Number of days:1. Number of samples:1"]][[1]])
rep1_scen1<-as.data.frame(DF[["Number of shools:20. STH Prevalence:0.02. Number of days:2. Number of samples:1"]][[1]])


temp<-rep1_scen1%>%group_by(ID.sub)%>%summarise(count1.mean.extra=mean(count1.mean), School.name.extra=mean(School.name), ID.extra=mean(ID.sub), DX.extra=mean(DX))



#Prepare dataset for models with only one measure Rogen-Gladen and MoM
#We can use the maximum count or the average for all days and samples. Then we should do a round
tempo<-DF
for (kk in 1:length(DF)) {
  cat("Dataset number ");cat(kk, sep="\n")
  for (mm in reps) {
  tempo[[kk]][[mm]]<-DF[[grid.sim[[kk]][["name"]]]][[mm]]%>%as.data.frame%>%group_by(ID.sub)%>%summarise(count1.mean.extra=mean(count1.mean), School.name.extra=mean(School.name), ID.extra=mean(ID.sub), DX.extra=mean(DX))
  }
  cat("\n")
}
save.image(file = "DFandTEMPO.RData")

load(file = "DFandTEMPO.RData")
library(prevalence)

results <- DF

class(DF[["Number of shools:20. STH Prevalence:0.01. Number of days:1. Number of samples:1"]][[1]])

dataDF=DF
gridsDF=grid.sim
model= "Bayesian Rogen-Gladen"
length(grid.sim)
fit<-DF
fit_models <- function(dataDF=DF, gridsDF=grid.sim,  model) {
  # Fit model
  if (model == "Bayesian Rogen-Gladen") {
    dffit<-vector(mode='list', length=length(grid.sim))
    #in 1:length(grid.sim)
    #grid.sim[1]
    #grid.sim[18]

    for (kk in c(1,10,18) ) {  #KK IS THE NUMBER SCENARIOS
      cat("Scenario ",kk,"\n")
      dffit[[kk]]<-data.frame(matrix(nrow = length(grid.sim), ncol = 10))
      colnames(dffit[[kk]])<-c("theta.sd", "theta","theta.var","theta.se","LIC", "UIC", "n", "Scenario","Dataset n", "model")
      for (ii in 1:5) {  #REPS IS THE NUMBER OF REPETITION FOR EACH SCENARIO
        fit[[kk]][[ii]] <- prevalence::truePrev(x=sum(dataDF[[ gridsDF[[kk]][["name"]] ]][[ii]][[1]][["TestR1"]] ),
                                                n=length(dataDF[[ grid.sim[[kk]][["name"]] ]][[ii]][[1]][["TestR1"]]))
        cat("Scenario ",kk," Rep #:",ii,"n: ",length(dataDF[[ grid.sim[[kk]][["name"]] ]][[ii]][[1]][["TestR1"]]),"\n")
        as<-fit[[kk]][[ii]]
        dffit[[kk]][ii,]<-c(theta.sd=as.numeric(summary(as)$TP[3,4]),theta=as.numeric(summary(as)$TP[3,1]),
                            theta.var=as.numeric(summary(as)$TP[3,5]), theta.se=as.numeric(sqrt(summary(as)$TP[3,5])),
                            LIC=as.numeric(summary(as)$TP[3,6]),UIC=as.numeric(summary(as)$TP[3,7]), n=as@par[["n"]], name=gridsDF[[kk]][["name"]],dataset.n=ii, model=model)


      }
    }
#    relhaz2 <- do.call(
#      rbind.data.frame,
#      dffit
#    )
#    view(relhaz2)

  }else if (model == "Method of Moment") {
    fit <- rstpm2::stpm2(Surv(t, d) ~ x, dataDF = dataDF, df = 2)
  } else { #Barenbold model
    fit <- eha::phreg(Surv(t, d) ~ x, dataDF = dataDF, dist = "weibull", shape = 1)
  }
  # Return relevant coefficients
  data.frame(
    dataset = unique(dataDF$dataset),
    n = unique(dataDF$n),
    baseline = unique(dataDF$baseline),
    theta = coef(fit)["x"],
    se = sqrt(ifelse(model == "Exp", fit$var["x", "x"], vcov(fit)["x", "x"])),
    model = model,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

relhaz <- do.call(
  rbind.data.frame,
  dffit
)
view(relhaz2)
class(relhaz)
relhaz<-relhaz[!is.na(relhaz$theta),]
relhaz$theta<-as.numeric(relhaz$theta)
relhaz$theta.se<-as.numeric(relhaz$theta.se)
relhaz$theta.sd<-as.numeric(relhaz$theta.sd)
relhaz$theta.var<-as.numeric(relhaz$theta.var)
relhaz$LIC<-as.numeric(relhaz$LIC)
relhaz$UIC<-as.numeric(relhaz$UIC)
relhaz$n<-as.numeric(relhaz$n)
class(relhaz$n)
str(relhaz)
relhaz2<-subset(relhaz,select = c(theta,theta.se,model,n,Scenario))
library(rsimsum)
s <- rsimsum::simsum(data = relhaz2, estvarname = "theta", se = "theta.se", true = 0.02, methodvar = "model", ref = "Bayesian Rogen-Gladen", by = c("n", "Scenario"))



for (kk in 1:1) {
  #cat("Dataset number ");cat(kk, sep="\n")


  #  x	  The apparent number of positive samples
  #  n	  The sample size

  results[[ grid.sim[[kk]][["name"]] ]] <- do.call(
    rbind.data.frame,
    lapply(
      X = DF[[ grid.sim[[kk]][["name"]] ]],
      FUN = prevalence::truePrev
    )
  )
  #cat("\n")
}

class(results)
relhaz <- do.call(
  rbind.data.frame,
  results)
row.names(relhaz) <- NULL


#fff<-DF[[grid.sim[[kk]][["name"]]]][[2]]
# count1.mean and count10.mean is the rounded average for each subject or pool10 through all days and samples.
# This is requiered to produce data for the models that only use only one daya and sample like Rogen Gladen  and Levecke MoM
#  if (params$plotss[1]){
#    par(mfrow=c(4,2))
#    hist(datatest$Mi, main="Mean egg intensity between individuals", xlab = expression(mu))
#    hist(datatest$Mi1, main="Mean egg intensity for individual sample", xlab = expression(mu[i1]))
#    hist(datatest$Mi10, main="Mean egg intensity for 10-pooled sample", xlab = expression(mu[i10]))
#    hist(datatest$Mids1, main="Mean egg intensity for Mids1", xlab = expression(mu[ids1]))
#    hist(datatest$Mids10, main="Mean egg intensity for Mids10", xlab = expression(mu[ids10]))
#    hist(datatest$count1)
#    hist(datatest$count10)
#  }








