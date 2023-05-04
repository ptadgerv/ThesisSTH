rm(list=ls())
setwd("~/Master EPI/SEM6/Thesis/DATA_LEVECKE/Leta et al., 2018/")
ind <- read.csv("GT_Ind_19JAN2016.csv")
pool <- read.csv("GT_Pool_19JAN2016.csv")

## Organizing Data
library(plyr)
library(xlsx)
dim(ind)
ind$n <- 1

# EPG
ind$EPG_AL <- ind$Eggs.of.Ascaris*24
ind$EPG_HK <- ind$Eggs.of.Hookworm*24
ind$EPG_TT <- ind$Eggs.of.Trichuris*24
ind$EPG_SM <- ind$Eggs.of.Schistosoma*24

pool$EPG_AL_P <- pool$Eggs.of.Ascaris*24
pool$EPG_HK_P <- pool$Eggs.of.Hookworm*24
pool$EPG_TT_P <- pool$Eggs.of.Trichuris*24
pool$EPG_SM_P <- pool$Eggs.of.Schistosoma*24

# positive/negative
ind$P_AL <- ifelse(ind$Eggs.of.Ascaris==0,0,1)
ind$P_HK <- ifelse(ind$Eggs.of.Hookworm==0,0,1)
ind$P_TT <- ifelse(ind$Eggs.of.Trichuris==0,0,1)
ind$P_SM <- ifelse(ind$Eggs.of.Schistosoma==0,0,1)
ind$mix <- ind$P_AL + ind$P_HK + ind$P_TT + ind$P_SM
ind$sth <- ind$P_AL + ind$P_HK + ind$P_TT
ind$P_sth <- ifelse(ind$sth==0,0,1)
ind$P_ntd <- ifelse(ind$mix==0,0,1)

pool$P_AL_P <- ifelse(pool$Eggs.of.Ascaris==0,0,1)
pool$P_HK_P <- ifelse(pool$Eggs.of.Hookworm==0,0,1)
pool$P_TT_P <- ifelse(pool$Eggs.of.Trichuris==0,0,1)
pool$P_SM_P <- ifelse(pool$Eggs.of.Schistosoma==0,0,1)

## intensity of infection
ind$Int_AL <- ifelse(ind$EPG_AL == 0, 0, ifelse(ind$EPG_AL < 5000, 1, ifelse(ind$EPG_AL>49999,3,2)))
ind$Int_TT <- ifelse(ind$EPG_TT == 0, 0, ifelse(ind$EPG_TT < 1000, 1, ifelse(ind$EPG_TT>9999,3,2)))
ind$Int_HK<- ifelse(ind$EPG_HK == 0, 0, ifelse(ind$EPG_HK < 2000, 1, ifelse(ind$EPG_TT>3999,3,2)))
ind$Int_SM<- ifelse(ind$EPG_SM == 0, 0, ifelse(ind$EPG_SM < 100, 1, ifelse(ind$EPG_TT>399,3,2)))

### Descriptive statistics
sum(ind$P_ntd)
round(100*mean(ind$P_ntd),1)
round(100*mean(ind$P_AL),1)
round(100*mean(ind$P_TT),1)
round(100*mean(ind$P_HK),1)
round(100*mean(ind$P_SM),1)

sum(ind$P_ntd)
sum(ind$P_AL)
sum(ind$P_TT)
sum(ind$P_HK)
sum(ind$P_SM)

round(mean(ind$EPG_AL),1)
round(mean(ind$EPG_TT),1)
round(mean(ind$EPG_HK),1)
round(mean(ind$EPG_SM),1)
table(ind$Int_AL)
table(ind$Int_TT)
table(ind$Int_HK)
table(ind$Int_SM)



sum <- ddply(ind, .(Woreda.code,School.Code), summarize,
             n = sum(n),schi_p = round(100*mean(P_SM),1),asc_p= round(100*mean(P_AL),1),hoo_p = round(100*mean(P_HK),1),tri_p= round(100*mean(P_TT),1),
             schi = round(mean(EPG_SM),1),asc= round(mean(EPG_AL),1),hoo = round(mean(EPG_HK),1),tri= round(mean(EPG_TT),1))

sum$ntd <- sum$schi_p + sum$asc_p + sum$hoo_p + sum$tri_p

min(sum$schi_p[sum$schi_p>0])
max(sum$schi_p[sum$schi_p>0])
min(sum$schi[sum$schi_p>0])
max(sum$schi[sum$schi_p>0])

min(sum$asc_p[sum$asc_p>0])
max(sum$asc_p[sum$asc_p>0])
min(sum$asc[sum$asc_p>0])
max(sum$asc[sum$asc_p>0])

min(sum$tri_p[sum$tri_p>0])
max(sum$tri_p[sum$tri_p>0])
min(sum$tri[sum$tri_p>0])
max(sum$tri[sum$tri_p>0])

min(sum$hoo_p[sum$hoo_p>0])
max(sum$hoo_p[sum$hoo_p>0])
min(sum$hoo[sum$hoo_p>0])
max(sum$hoo[sum$hoo_p>0])

ntd <- subset(sum, sum$ntd==0)
write.xlsx(sum, 'sum.xlsx')


### Merging individual and pooled samples
sum2 <- ddply(ind, .(Woreda.code,School.Code,Pool.ID), summarize,
              n = sum(n),schi = mean(EPG_SM),asc= mean(EPG_AL),hoo= mean(EPG_HK),tri= mean(EPG_TT),
              min = sum(Reading.time..MIN.),  sec= sum(Reading.time..SEC., na.rm=T))

pool[1:10,]
sum2[1:10,]

TOT <- merge(sum2,pool,by=c("Woreda.code","School.Code",'Pool.ID'), all.Y = T)
TOT[1:10,]
TOT$N <- 1


library(Hmisc)
var.labels = c(Woreda.code="Woreda code",
               School.Code="School Code",
               Pool.ID="Pool ID",
               n="n",
               schi="schi",
               asc="asc",
               hoo="hoo",
               tri="tri",
               min="Time min",
               sec="Time sec",
               Woreda.name="Woreda name",
               School.Name="School Name",
               Eggs.of.Schistosoma="Eggs of Schistosoma",
               Eggs.of.Ascaris="Eggs of Ascaris",
               Eggs.of.Hookworm="Eggs of Hookworm",
               Eggs.of.Trichuris="Eggs of Trichuris",
               other="other",
               Reading.time..MIN.="Reading time MIN",
               Reading.time..SEC.="Reading time SEC",
               X="X",
               X.1="X.1",
               X.2="X.2",
               X.3="X.3",
               X.4="X.4",
               X.5="X.5",
               EPG_AL_P="EPG_AL_P",
               EPG_HK_P="EPG_HK_P",
               EPG_TT_P="EPG_TT_P",
               EPG_SM_P="EPG_SM_P",
               P_AL_P="P_AL_P",
               P_HK_P="P_HK_P",
               P_TT_P="P_TT_P",
               P_SM_P="P_SM_P",
               N="N")

label(TOT) = as.list(var.labels[match(names(TOT), names(var.labels))])

label(TOT)

library(prevalence)

truePrev(x = sum(TOT$P_AL), n = nrow(TOT),SE = ~dunif(0.60, 1.00), SP = ~dunif(0.9, 1.00))
truePrev(x = sum(TOT$P_HK), n = nrow(TOT),SE = ~dunif(0.60, 1.00), SP = ~dunif(0.9, 1.00))
truePrev(x = sum(TOT$P_TT), n = nrow(TOT),SE = ~dunif(0.60, 1.00), SP = ~dunif(0.9, 1.00))
truePrev(x = sum(TOT$P_SM), n = nrow(TOT),SE = ~dunif(0.60, 1.00), SP = ~dunif(0.9, 1.00))

prevalence::truePrevPools(x=as.integer(TOT$P_AL_P), n=as.integer(TOT$n),  prior = c(1, 1),nchains = 2, burnin = 10000, update = 10000,verbose = FALSE,SE = ~dunif(0.60, 1.00), SP = ~dunif(0.9, 1.00))
prevalence::truePrevPools(x=as.integer(TOT$P_HK_P), n=as.integer(TOT$n),SE = ~dunif(0.60, 1.00), SP = ~dunif(0.9, 1.00), prior = c(1, 1),nchains = 2, burnin = 10000, update = 10000,verbose = FALSE)
prevalence::truePrevPools(x=as.integer(TOT$P_TT_P), n=as.integer(TOT$n),SE = ~dunif(0.60, 1.00), SP = ~dunif(0.9, 1.00), prior = c(1, 1),nchains = 2, burnin = 10000, update = 10000,verbose = FALSE)
prevalence::truePrevPools(x=as.integer(TOT$P_SM_P), n=as.integer(TOT$n),SE = ~dunif(0.60, 1.00), SP = ~dunif(0.9, 1.00), prior = c(1, 1),nchains = 2, burnin = 10000, update = 10000,verbose = FALSE)

