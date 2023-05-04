#######################################################################################
################################# prevalence structure  ###############################
#################################for school_j & prevalence_i ##########################
#######################################################################################
rm(list = ls())
library(dplyr)
set.seed(2023)
gridd<-expand.grid(schools=c(20,50), prev.STH=c(0.01,0.02,0.03),var.STH=0.001,pooln1=10, CV.i=1.5,mu=1000, CV.d=0.75,   pooln2=25, CV.d.10=0.6,CV.d.25=0.5, CV.s=0.25, plotss=FALSE, TPR=0.88*0.88,TNR=0.97*0.97, n.s=150, n.days=c(1,2, 3), n.sample=c(1,2,3))
gridd$nnn=gridd$schools*gridd$n.s
gridd$Scenario<-rownames(gridd)
gridd$name<-paste0("Number of shools:",gridd$schools,". STH Prevalence:",gridd$prev.STH,". Number of days:",gridd$n.days,". Number of samples:",gridd$n.sample)



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


  #################################XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  #################################XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  #################################XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  #################################XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  #################################XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  #################################XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
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

  mids.count.error<-5
  for (ii in 1:(nnn*n.days*n.sample)) {
    datatest[[ji]]$count1[ii]<-ifelse(datatest[[ji]]$Class2x2.ideal[ii]=="FP",
                                      rpois(1,  lambda=(datatest[[ji]]$Mids1[ii]+mids.count.error)),
                                      rzipois(1,  lambda=datatest[[ji]]$Mids1[ii],pstr0 = 0.15))
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
  temp.10<-datatest[[ji]]%>%group_by(School.name,pool10)%>%summarise(count10.mean=round(mean(count10),digits=0))  #NOT SURE by whom should I round

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
#for (kk in 1:length(grid.sim)) {
for (kk in 1:1) {
#cat("Dataset number ");cat(kk, sep="\n")
DF[[grid.sim[[kk]][["name"]] ]]<-lapply(
  X = reps,FUN = simSTH,
  params = list(grid.sim[[kk]]))
#cat("\n")
}
save.image(file = "Scenario1_100rep.RData")
#save.image(file = "SesionSimu.RData")
getwd()
load("Scenario1_100rep.RData")

rep1_scen1<-as.data.frame(DF[["Number of shools:20. STH Prevalence:0.01. Number of days:1. Number of samples:1"]][[1]])
temp<-rep1_scen1%>%group_by(School.name,ID.sub)%>%summarise(Dx.mean=formatC(mean(DX), format = "e", digits = 2)) ### DX  is consistent

rep1_scen1%>%group_by(School.name,ID.sub)%>%mutate(count1.mean.extra = min(count1.mean)) %>%ungroup()
rep1_scen1.2<-rep1_scen1
rep1_scen1.2%>%group_by(School.name,ID.sub)%>%mutate(count1.mean.min = min(count1.mean),count1.mean.mean = mean(count1.mean) )

rep1_scen1<-as.data.frame(DF[["Number of shools:20. STH Prevalence:0.01. Number of days:1. Number of samples:1"]][[1]])
temp<-rep1_scen1%>%group_by(School.name,ID.sub)%>%summarise(count1.mean.extra=mean(count1.mean), School.name.extra=mean(School.name), ID.extra=mean(ID.sub), DX.extra=mean(DX))


# count1.mean and count10.mean is the rounded average for each subject or pool10 through all days and samples.
# This is requiered to produce data for the models that only use only one daya and sample like Rogen Gladen  and Levecke MoM

dataDF1<-as.data.frame(DF[["Number of shools:20. STH Prevalence:0.01. Number of days:1. Number of samples:1"]][[1]])
write.csv(dataDF1, "dataDF1.csv")

#sum(!(DF$ID.sub %in% 1:nrow(datatest)))
datatest$count1.mean<-datatest[[1]]$count10.mean<-NA
for (jj in 1:nrow(datatest)) {
  datatest$count1.mean[datatest$ID.sub==jj]<-DF$count1.mean[DF$ID.sub==jj][1]
  datatest$count10.mean[datatest$ID.sub==jj]<-DF$count10.mean[DF$ID.sub==jj][1]
}

sum(!(DF$ID.sub %in% datatest$ID.sub))
data.table::fwrite(DF,"DF3.30APR23.csv")



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

#Prepare dataset for models with only one measure Rogen-Gladen and MoM
#We can use the maximum count or the average for all days and samples. Then we should do a round


library(haven)
haven::write_xpt(datatest,"datatest.sas7bdat")
haven::write_sav(datatest,"datatest.sav")




##### For the different DGM meachanism
reps<-1:100 #We generate 100 (reps) datasets for each data-generating mechanism
data <- list()
data[["n = 50, baseline = NB, mu = 500"]] <- lapply(
  X = reps,
  FUN = simSTH,
  n = 50,
  params = list(  CV.i=1.5,mu=500, CV.d=0.75, pooln1=10, pooln2=25, CV.d.10=0.6,CV.d.25=0.5, CV.s=0.25, plotss=FALSE)
)

data[["n = 250, baseline = NB, mu = 500"]] <- lapply(
  X = reps,
  FUN = simSTH,
  n = 250,
  #baseline = "NegativeBinomial",
  params = list(  CV.i=1.5,mu=500, CV.d=0.75, pooln1=10, pooln2=25, CV.d.10=0.6,CV.d.25=0.5, CV.s=0.25, plotss=FALSE)
)

data[["n = 50, baseline = NB, mu = 1000"]] <- lapply(
  X = reps,
  FUN = simSTH,
  n = 50,
  params = list(  CV.i=1.5,mu=1000, CV.d=0.75, pooln1=10, pooln2=25, CV.d.10=0.6,CV.d.25=0.5, CV.s=0.25, plotss=FALSE)
)

data[["n = 250, baseline = NB, mu = 1000"]] <- lapply(
  X = reps,
  FUN = simSTH,
  n = 250,
  #baseline = "NegativeBinomial",
  params = list(  CV.i=1.5,mu=1000, CV.d=0.75, pooln1=10, pooln2=25, CV.d.10=0.6,CV.d.25=0.5, CV.s=0.25, plotss=FALSE)
)




#The next function is to write datsets to bne use in Winbugs

if (1==2) {  #NOT ALWAYS REQUIERED


  "writeDatafileR" <-
    function(DATA, towhere = "toWinBUGS.txt", fill = 80)
    {
      #
      # Writes from R to file "towhere" text defining a list containing "DATA" in a form compatable with WinBUGS.
      # Required arguments:
      # DATA - either a data frame or else a list consisting of any combination of scalars, vectors, arrays or data frames (but not lists).
      #   If a list, all list elements that are not data.frames must be named. Names of data.frames in DATA are ignored.
      # Optional arguments:
      # towhere - file to receive output. Is "toWinBUGS.txt" by default.
      # fill - If numeric, number of columns for output. (Default is 80.) If FALSE, output will be on one line. If TRUE, number of
      #   columns is given by .Options$width.
      # Value:
      # Text defining a list is output to file "towhere".
      # Details:
      #  The function performs considerable checking of DATA argument. Since WinBUGS requires numeric input, no factors or character vectors
      # are allowed. All data must be named, either as named elements of DATA (if it is a list) or else using the names given in data frames.
      # Data frames may contain matrices.
      # Arrays of any dimension are rearranged to be in row-major order, as required by WinBUGS. Scientific notation is also handled properly.
      # In particular, the number will consist of a mantissa _containing a decimal point_ followed by "E", then either "+" or "-", and finally
      # a _two-digit_ number.
      # Written by Terry Elrod. Disclaimer: This function is used at the user's own risk.
      # Please send comments to Terry.Elrod@UAlberta.ca.
      # Revision history: 2003-11-14: Fixed to handle missing values properly. (Thanks to Kjetil Halvorsen.)
      #					2003-11-14:	Tests for valid Winbugs names. Forces single precision for all numbers.
      formatDataR <-
        #
        # Prepared DATA for input to WinBUGS.
        function(DATA)
        {
          testWinbugsNames <-
            #
            # Checks to see that all names are valid...
            function(na){
              baseTestString <- c(
                "The following variable names are invalid in R: ",
                "The following variable names are used more than once: ",
                "The following variable names have more than 8 characters: ",
                "The following variable names contain two or more periods in a row: ",
                "The following variable names end in a period: ")
              # Testing for invalid R names ...
              nameTest1 <- make.names(na, unique = FALSE)
              nameTest1 <- (nameTest1 != na)
              # Testing for duplicate names....
              nameTest2 <- make.names(na, unique = TRUE)
              nameTest2 <- (nameTest2 != na)
              # Testing for excess length...
              nameTest3 <- substring(na, 1, 8)
              nameTest3 <- (na != nameTest3)
              # Testing for presence of two or more successive periods ...
              nameTest4 <- regexpr("\\.\\.", na)
              nameTest4 <- (nameTest4 > 0)
              # Testing for presence of ending period ...
              nameTest5 <- regexpr("\\.$", na)
              nameTest5 <- (nameTest5 > 0)
              # Assembling tests and reporting results...
              nameTest <- cbind(nameTest1, nameTest2, nameTest3, nameTest4, nameTest5)
              if(any(nameTest)){
                nameTestInd <- apply(nameTest, 2, any)
                whichTest <- seq(along=nameTestInd)[nameTestInd]
                testString <- "There were problems with names of one or more variables:"
                if(nameTestInd[1])
                  testString <- paste(testString, paste(baseTestString[1], paste(na[nameTest[,1]], collapse = ", "), sep=""), sep="\n")
                if(nameTestInd[2])
                  testString <- paste(testString, paste(baseTestString[2], paste(unique(na[nameTest[,2]]), collapse = ", "), sep=""), sep="\n")
                if(nameTestInd[3])
                  testString <- paste(testString, paste(baseTestString[3], paste(na[nameTest[,3]], collapse = ", "), sep="") ,sep="\n")
                if(nameTestInd[4])
                  testString <- paste(testString, paste(baseTestString[4], paste(na[nameTest[,4]], collapse = ", "), sep="") ,sep="\n")
                if(nameTestInd[5])
                  testString <- paste(testString, paste(baseTestString[5], paste(na[nameTest[,5]], collapse = ", "), sep="") ,sep="\n")
                stop(testString)
              }
              invisible(0)
            }
          toSingle <-
            #
            # Takes numeric vector, adds period to mantissa in scientific notation (if necessary),
            #	converts "e" to "E", expresses mantissa with at most 10 characters,
            #	and eliminates trailing zeros from mantissa.
            function(x)
            {
              myRegMatchPos <-
                #
                # Duplicates regMatchPos in the S4 engine...
                function(w, txt)
                {
                  st <- regexpr(txt, w)
                  pplusind <- (st > 0)
                  fin <- st + attr(st, "match.length") - 1
                  pplus <- cbind(st, fin)
                  pplus[!pplusind,  ] <- NA
                  pplus
                }
              xdim <- dim(x)
              x <- as.single(x)
              x <- sapply(x,function(y) format(y, digits=7, trim=TRUE))
              # First to look for positives:
              pplus <- myRegMatchPos(x, "e\\+0")
              pplusind <- apply(pplus, 1, function(y)
                (!any(in.sa(y))))
              if(any(pplusind)) {
                # Making sure that periods are in mantissa...
                init <- substring(x[pplusind], 1, pplus[
                  pplusind, 1] - 1)
                #...preceeding exponent
                pper <- myRegMatchPos(init, "\\.")
                pperind <- apply(pper, 1, function(y)
                  (all(in.sa(y))))
                if(any(pperind))
                  init[pperind] <- paste(init[pperind],
                                         ".0", sep = "")
                # Changing the format of the exponent...
                x[pplusind] <- paste(init, "E+", substring(
                  x[pplusind], pplus[pplusind, 2] + 1),
                  sep = "")
              }
              # Then to look for negatives:
              pminus <- myRegMatchPos(x, "e\\-0")
              pminusind <- apply(pminus, 1, function(y)
                (!any(in.sa(y))))
              if(any(pminusind)) {
                # Making sure that periods are in mantissa...
                init <- substring(x[pminusind], 1, pminus[
                  pminusind, 1] - 1)
                #...preceeding exponent
                pper <- myRegMatchPos(init, "\\.")
                pperind <- apply(pper, 1, function(y)
                  (all(in.sa(y))))
                if(any(pperind))
                  init[pperind] <- paste(init[pperind],
                                         ".0", sep = "")
                # Changing the format of the exponent...
                x[pminusind] <- paste(init, "E-", substring(
                  x[pminusind], pminus[pminusind, 2] +
                    1), sep = "")
              }
              x
            }
          if(!is.list(DATA))
            stop("DATA must be a named list or data frame.")
          dlnames <- names(DATA)
          if(is.data.frame(DATA))
            DATA <- as.list(DATA)
          #
          # Checking for lists in DATA....
          lind <- sapply(DATA, is.list)
          # Checking for data frames in DATA....
          dfind <- sapply(DATA, is.data.frame)
          # Any lists that are not data frames?...
          if(any(lind & !dfind)) stop("DATA may not contain lists.")
          # Checking for unnamed elements of list that are not data frames....
          if(any(dlnames[!dfind] == "")) stop(
            "When DATA is a list, all its elements that are not data frames must be named."
          )
          if(any(dfind)) {
            dataold <- DATA
            DATA <- vector("list", 0)
            for(i in seq(along = dataold)) {
              if(dfind[i])
                DATA <- c(DATA, as.list(dataold[[i]]))
              else DATA <- c(DATA, dataold[i])
            }
            dataold <- NULL
          }
          dlnames <- names(DATA)
          # Making sure all names are valid ...
          testWinbugsNames(dlnames)
          # Checking for factors....
          factorind <- sapply(DATA, is.factor)
          if(any(factorind))
            stop(paste(
              "DATA may not include factors. One or more factor variables were detected:",
              paste(dlnames[factorind], collapse = ", ")))
          # Checking for character vectors....
          charind <- sapply(DATA, is.character)
          if(any(charind))
            stop(paste(
              "WinBUGS does not handle character data. One or more character variables were detected:",
              paste(dlnames[charind], collapse = ", ")))
          # Checking for complex vectors....
          complexind <- sapply(DATA, is.complex)
          if(any(complexind))
            stop(paste(
              "WinBUGS does not handle complex data. One or more complex variables were detected:",
              paste(dlnames[complexind], collapse = ", ")))
          # Checking for values farther from zero than 1E+38 (which is limit of single precision)....
          toobigind <- sapply(DATA, function(x)
          {
            y <- abs(x[!in.sa(x)])
            any(y[y > 0] > 1.0e+038)
          }
          )
          if(any(toobigind))
            stop(paste(
              "WinBUGS works in single precision. The following variables contain data outside the range +/-1.0E+38: ",
              paste(dlnames[toobigind], collapse = ", "),
              ".\n", sep = ""))
          # Checking for values in range +/-1.0E-38 (which is limit of single precision)....
          toosmallind <- sapply(DATA, function(x)
          {
            y <- abs(x[!in.sa(x)])
            any(y[y > 0] < 1.0e-038)
          }
          )
          n <- length(dlnames)
          data.string <- as.list(rep(NA, n))
          for(i in 1:n) {
            ldi <- length(DATA[[i]])
            if(ldi == 1) {
              ac <- toSingle(DATA[[i]])
              data.string[[i]] <- c(
                names(DATA)[i],
                "=",
                paste(ac),
                "," )
              next
            }
            if(is.vector(DATA[[i]]) & ldi > 1) {
              ac <- toSingle(DATA[[i]])
              data.string[[i]] <- c(
                names(DATA)[i],
                "= c(",
                paste(ac[-ldi], ",", sep=""),
                paste(ac[ldi], ")", sep=""),
                "," )
              next
            }
            if(is.array(DATA[[i]])) {
              ac <- toSingle(aperm(DATA[[i]]))
              data.string[[i]] <- c(
                names(DATA)[i],
                "= structure(.Data = c(",
                paste(ac[-ldi], ",", sep=""),
                paste(ac[ldi], "),", sep=""),
                ".Dim=c(",
                paste(as.character(dim(DATA[[i]])),collapse = ", "),
                "))",
                "," )
            }
          }
          data.string <- unlist(data.string)
          data.tofile <- c(
            "list(",
            data.string[-length(data.string)],
            ")" )
          if(any(toosmallind))
            warning(paste(
              "WinBUGS works in single precision. The following variables contained nonzero data",
              "\ninside the range +/-1.0E-38 that were set to zero: ",
              paste(dlnames[toosmallind], collapse = ", "),
              ".\n", sep = ""))
          return(data.tofile)
        }
      cat(formatDataR(DATA), file = towhere, fill = fill)
      formatDataR(DATA)
      invisible(0)
    }





pasalo<-data1[["n = 250, baseline = NB, mu = 1000"]][[1]]
library('dplyr')
# dplyr - Select columns by label name & gender
names(pasalo)
pasalo1<-pasalo  %>% select('Mids1','TestR1', 'count1')
names(pasalo1)<-c("Mids",  "TestR", "count")
writeDatafileR(data1[["n = 250, baseline = NB, mu = 1000"]][[1]])
pasalo1$Mids<- trunc(pasalo1$Mids)*24
writeDatafileR(pasalo1)
}