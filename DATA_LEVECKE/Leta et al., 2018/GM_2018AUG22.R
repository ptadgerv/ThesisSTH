### Importing of data
###
###
#GT_Ind <- read.csv("~/Master EPI/SEM6/Thesis/DATA_LEVECKE/Leta et al., 2018/GT_Ind_19JAN2016.csv")
#GT_Pool <- read.csv("~/Master EPI/SEM6/Thesis/DATA_LEVECKE/Leta et al., 2018/GT_Pool_19JAN2016.csv")

ind <- read.csv("~/Master EPI/SEM6/Thesis/DATA_LEVECKE/Leta et al., 2018/GT_Ind_19JAN2016.csv")
pool <- read.csv("~/Master EPI/SEM6/Thesis/DATA_LEVECKE/Leta et al., 2018/GT_Pool_19JAN2016.csv")
#GT_Ind$Sub.ID %in% GT_Pool$Pool.ID
#sum(GT_Ind$Pool.ID %in% GT_Pool$Pool.ID)
#unique(GT_Ind$Pool.ID)
#prepool <- read.csv('/Users/blevecke/Dropbox/PhDstudents/Gemechu/Gemechu_data/Final data/GT_Prep_Pool_19JAN2016.csv', header = TRUE, sep = ";")
#prepkk <- read.csv('/Users/blevecke/Dropbox/PhDstudents/Gemechu/Gemechu_data/Final data/GT_Prep_KK_19JAN2016.csv', header = TRUE, sep = ";")
#team <- read.csv('/Users/blevecke/Dropbox/Communications/Publicaties/Drafts/Pooling/Human/Ethiopia/Gemechu/Final Data/Teams.csv', header = TRUE, sep = ";")

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



### timing
team <- read.csv('/Users/blevecke/Dropbox/Communications/Publicaties/Drafts/Pooling/Human/Ethiopia/Gemechu/Final Data/teams.csv', header = TRUE, sep = ";")
Team <- team$Team.numbert
Woreda.code <- team$Woreda.code
team <- data.frame(Team,Woreda.code)

team2 <- ddply(team, .(Woreda.code,Team), summarize,
                  n = sum(n, na.rm=T))

## Individual reading
readind <- merge(ind, team2, by=c("Woreda.code"))
dim(readind)
readind$tot <- readind$Reading.time..MIN.*60 + readind$Reading.time..SEC.
readind$n <- 1
dim(readind)
readsum1 <- ddply(readind, .(Team), summarize, n = sum(n),sum_ind = sum(tot), mean_ind = round(mean(tot)/60,1),sd_ind = round(sd(tot)/60,1))
mean(readind$tot)/60
readind55 <- subset(readind, !Team==12)
readsum55i <- ddply(readind55, .(Team), summarize, n = sum(n),sum_ind = sum(tot), mean_ind = round(mean(tot)/60,1),sd_ind = round(sd(tot)/60,1))


## Pooled reading
readpool <- merge(pool, team2, by=c("Woreda.code"))
dim(readpool)
readpool$tot <- readpool$Reading.time..MIN.*60 + readpool$Reading.time..SEC.
readpool$n <- 1
dim(readpool)
readsum2 <- ddply(readpool, .(Team), summarize, n = sum(n),sum_pool=sum(tot),mean_pool = round(mean(tot)/60,1),sd_pool = round(sd(tot)/60,1))
mean(readpool$tot)/60
readpool55 <- subset(readpool, !Team==12)
readsum55p <- ddply(readpool55, .(Team), summarize, n = sum(n),sum_pool=sum(tot),mean_pool = round(mean(tot)/60,1),sd_pool = round(sd(tot)/60,1))

## Preparing KK
prep <- merge(prepkk, team2, by=c("Woreda.code"))
dim(prep)
prep$tot <- prep$Preparation.time..MIN.*60 + prep$Preparation.time..SEC.
prep$n <- 1
dim(prep)
prepsum1 <- ddply(prep, .(Team), summarize, n = sum(n),sum_prepkk = sum(tot),mean_prepkk = round(mean(tot)/60,1),sd_prepkk = round(sd(tot)/60,1))
prep55 <- subset(prep, !Team==12)
prepsum55 <- ddply(prep55, .(Team), summarize, n = sum(n),sum_prepkk = sum(tot),mean_prepkk = round(mean(tot)/60,1),sd_prepkk = round(sd(tot)/60,1))


## Preparing pools
prepPool <- merge(prepool, team2, by=c("Woreda.code"))
dim(prepPool)
prepPool$tot <- prepPool$Preparation.time..MIN.*60 + prepPool$Preparation.time..SEC.
prepPool$n <- 1
dim(prepPool)
prepsum2 <- ddply(prepPool, .(Team), summarize, n = sum(n),sum_preppool=sum(tot),mean_prepool = round(mean(tot)/60,1),sd_preppool = round(sd(tot)/60,1))
mean(prepPool$tot)/60
prepPool55 <- subset(prepPool, !Team==12)
prepsum55p <- ddply(prepPool55, .(Team), summarize, n = sum(n),sum_preppool=sum(tot),mean_prepool = round(mean(tot)/60,1),sd_preppool = round(sd(tot)/60,1))

timing <- cbind(readsum1,readsum2,prepsum1,prepsum2)
timing55 <- cbind(readsum55i,readsum55p,prepsum55,prepsum55p)

write.xlsx(timing55, '/Users/blevecke/Desktop/timing55.xlsx')


### Sensitivity at pooled level
## Ascaris
AL <- subset(TOT, TOT$asc > 0 |TOT$EPG_AL_P > 0)
dim(AL)
AL$inP <- ifelse(AL$asc>0,1,0)
AL$poolP <- ifelse(AL$EPG_AL_P>0,1,0)
round(100*mean(AL$inP),1)
round(100*mean(AL$poolP),1)


## CI
inPR <- rep(NA,5000)
poolPR <- rep(NA,5000)
for (j in 1:5000)
{
  inPR[j] <- mean(sample(AL$inP, length(AL[,1]), replace=T))
  poolPR[j] <- mean(sample(AL$poolP, length(AL[,1]), replace=T))
}
round(100*quantile(inPR, prob=c(0.025,0.975)),1)
round(100*quantile(poolPR, prob=c(0.025,0.975)),1)

## differences across examination strategy
diff<- mean(AL$inP) - mean(AL$poolP)
AL$diff <- AL$inP - AL$poolP
diffAL <- rep(NA,5000)
for (j in 1:5000)
{
  t <- c(rep(1,length(AL[,1])/2),rep(0,length(AL[,1])/2))
  toss <- sample(t, length(AL[,1]), replace=FALSE)
  for (i in 1:length(AL[,1]))
  {
    ifelse(toss[i] == 1, AL$diffR[i] <- AL$diff[i], AL$diffR[i] <- -AL$diff[i])
  }

  diffAL[j] <- abs(mean(AL$diffR))
}


hist(diffAL)
sum(ifelse(diffAL > abs(diff), 1, 0))/5000

## Variation across infection intensity
quantile(AL$asc, prob=c(0.33,0.66))
AL$int <- ifelse(AL$asc<=7.368,1,ifelse(AL$asc<=34.272,2,3))
alint <- tapply(AL$poolP,AL$int,mean)
tapply(AL$N,AL$int,sum)

## CI
poolPR1 <- rep(NA,5000)
poolPR2 <- rep(NA,5000)
poolPR3 <- rep(NA,5000)

for (j in 1:5000)
{
  poolPR1[j] <- mean(sample(AL$poolP[AL$int==1], length(AL$int[AL$int==1]), replace=T))
  poolPR2[j] <- mean(sample(AL$poolP[AL$int==2], length(AL$int[AL$int==2]), replace=T))
  poolPR3[j] <- mean(sample(AL$poolP[AL$int==3], length(AL$int[AL$int==3]), replace=T))
}
round(100*quantile(poolPR1, prob=c(0.025,0.975)),1)
round(100*quantile(poolPR2, prob=c(0.025,0.975)),1)
round(100*quantile(poolPR3, prob=c(0.025,0.975)),1)
x <- seq(0,100,1)
y <- seq(0,100,1)

par(mfrow=c(1,1),ps = 12, cex = 1, cex.main = 1)
plot(x,y, type='n', ylab='Sensitivity of a pooled examination strategy (%)', xlab="Level of egg excretion", xaxt = 'n')
axis(side = 1, at=c(33/2,33/2+33,33/2+66), lab=c('Mean FEC =< q33', 'q33 < mean FEC =< q66', 'mean FEC > q66'))
points(33/2-5,40.7, pch = 19)
lines(c(33/2-5,33/2-5),c(22.2,53.8))
points(33+33/2-5,42.3, pch = 19)
lines(c(33+33/2-5,33+33/2-5),c(23.1,61.5))
points(66+33/2-5,81.5, pch = 19)
lines(c(66+33/2-5,66+33/2-5),c(66.7,96.3))
lines(c(33/2-5,33+33/2-5,66+33/2-5), c(40.7,42.3,81.5))

## significant differences
alT <- rep(NA,5000)

for (j in 1:5000)
{
  AL$t <- sample(AL$poolP,length(AL[,1]), replace = FALSE)
  d13 <- abs(mean(AL$t[AL$int==1]) - mean(AL$t[AL$int==3]))
  d12 <- abs(mean(AL$t[AL$int==1]) - mean(AL$t[AL$int==2]))
  d23 <- abs(mean(AL$t[AL$int==2]) - mean(AL$t[AL$int==3]))
  alT[j] <- max(d13,d12,d23)
}

#hist(alT)
diff12 <- abs(alint[1]-alint[2])
sum(ifelse(alT>diff12,1,0))/5000

diff13 <- abs(alint[1]-alint[3])
sum(ifelse(alT>diff13,1,0))/5000

diff23 <- abs(alint[3]-alint[2])
sum(ifelse(alT>diff23,1,0))/5000


# Trichuris
TT <- subset(TOT, TOT$tri > 0 |TOT$EPG_TT_P > 0)
dim(TT)
TT$inP <- ifelse(TT$tri>0,1,0)
TT$poolP <- ifelse(TT$EPG_TT_P>0,1,0)
round(100*mean(TT$inP),1)
round(100*mean(TT$poolP),1)


## CI
inPR <- rep(NA,5000)
poolPR <- rep(NA,5000)
for (j in 1:5000)
{
  inPR[j] <- mean(sample(TT$inP, length(TT[,1]), replace=T))
  poolPR[j] <- mean(sample(TT$poolP, length(TT[,1]), replace=T))
}
round(100*quantile(inPR, prob=c(0.025,0.975)),1)
round(100*quantile(poolPR, prob=c(0.025,0.975)),1)

## differences across examination strategy
diff<- mean(TT$inP) - mean(TT$poolP)
TT$diff <- TT$inP - TT$poolP
diffTT <- rep(NA,5000)
for (j in 1:5000)
{
  t <- c(rep(1,length(TT[,1])/2),rep(0,length(TT[,1])/2))
  toss <- sample(t, length(TT[,1]), replace=FALSE)
  for (i in 1:length(TT[,1]))
  {
    ifelse(toss[i] == 1, TT$diffR[i] <- TT$diff[i], TT$diffR[i] <- -TT$diff[i])
  }

  diffTT[j] <- abs(mean(TT$diffR))
}

#hist(diffTT)
sum(ifelse(diffTT > abs(diff), 1, 0))/5000

# hookworm
HW <- subset(TOT, TOT$hoo > 0 |TOT$EPG_HK_P > 0)
dim(HW)
HW$inP <- ifelse(HW$hoo>0,1,0)
HW$poolP <- ifelse(HW$EPG_HK_P>0,1,0)
round(100*mean(HW$inP),1)
round(100*mean(HW$poolP),1)


## CI
inPR <- rep(NA,5000)
poolPR <- rep(NA,5000)
for (j in 1:5000)
{
  inPR[j] <- mean(sample(HW$inP, length(HW[,1]), replace=T))
  poolPR[j] <- mean(sample(HW$poolP, length(HW[,1]), replace=T))
}
round(100*quantile(inPR, prob=c(0.025,0.975)),1)
round(100*quantile(poolPR, prob=c(0.025,0.975)),1)

## differences across examination strategy
diff<- mean(HW$inP) - mean(HW$poolP)
HW$diff <- HW$inP - HW$poolP
diffHW <- rep(NA,5000)
for (j in 1:5000)
{
  t <- c(rep(1,length(HW[,1])/2),rep(0,length(HW[,1])/2))
  toss <- sample(t, length(HW[,1]), replace=FALSE)
  for (i in 1:length(HW[,1]))
  {
    ifelse(toss[i] == 1, HW$diffR[i] <- HW$diff[i], HW$diffR[i] <- -HW$diff[i])
  }

  diffHW[j] <- abs(mean(HW$diffR))
}

#hist(diffHW)
sum(ifelse(diffHW > abs(diff), 1, 0))/5000


## Variation across infection intensity
quantile(HW$hoo, prob=c(0.33,0.66))
HW$int <- ifelse(HW$hoo<=4.800,1,ifelse(HW$hoo<26.928,2,3))
HWint<- tapply(HW$poolP,HW$int,mean)
tapply(HW$N,HW$int,sum)

## CI
poolPR1 <- rep(NA,5000)
poolPR2 <- rep(NA,5000)
poolPR3 <- rep(NA,5000)

for (j in 1:5000)
{
  poolPR1[j] <- mean(sample(HW$poolP[HW$int==1], length(HW$int[HW$int==1]), replace=T))
  poolPR2[j] <- mean(sample(HW$poolP[HW$int==2], length(HW$int[HW$int==2]), replace=T))
  poolPR3[j] <- mean(sample(HW$poolP[HW$int==3], length(HW$int[HW$int==3]), replace=T))
}
round(100*quantile(poolPR1, prob=c(0.025,0.975)),1)
round(100*quantile(poolPR2, prob=c(0.025,0.975)),1)
round(100*quantile(poolPR3, prob=c(0.025,0.975)),1)

## significant differences
HWT <- rep(NA,5000)

for (j in 1:5000)
{
  HW$t <- sample(HW$poolP,length(HW[,1]), replace = FALSE)
  d13 <- abs(mean(HW$t[HW$int==1]) - mean(HW$t[HW$int==3]))
  d12 <- abs(mean(HW$t[HW$int==1]) - mean(HW$t[HW$int==2]))
  d23 <- abs(mean(HW$t[HW$int==2]) - mean(HW$t[HW$int==3]))
  HWT[j] <- max(d13,d12,d23)
}

#hist(HWT)
diff12 <- abs(HWint[1]-HWint[2])
sum(ifelse(HWT>diff12,1,0))/5000

diff13 <- abs(HWint[1]-HWint[3])
sum(ifelse(HWT>diff13,1,0))/5000

diff23 <- abs(HWint[3]-HWint[2])
sum(ifelse(HWT>diff23,1,0))/5000
#plot(x,y, type='n', ylab='Sensitivity of a pooled examination strategy (%)', xaxt="n", xlab = '')
#axis(side = 1, at=c(33/2,33/2+33,33/2+66), lab=c('Mean FEC =< q33', 'q33 < mean FEC =< q66', 'mean FEC > q66'))
points(33/2,52, pch = 19, col = 2)
lines(c(33/2,33/2),c(32,72), col = 2)
points(33+33/2,50, pch = 19, col = 2)
lines(c(33+33/2,33+33/2),c(30,70), col = 2)
points(66+33/2,82.6, pch = 19, col = 2)
lines(c(66+33/2,66+33/2),c(65.2,95.7), col = 2)
lines(c(33/2,33+33/2,66+33/2), c(52,50,82.6), col = 2)


### Schistosoma mansoni
SM <- subset(TOT, TOT$schi > 0 |TOT$EPG_SM_P > 0)
dim(SM)
SM$inP <- ifelse(SM$schi>0,1,0)
SM$poolP <- ifelse(SM$EPG_SM_P>0,1,0)
round(100*mean(SM$inP),1)
round(100*mean(SM$poolP),1)


## CI
inPR <- rep(NA,5000)
poolPR <- rep(NA,5000)
for (j in 1:5000)
{
  inPR[j] <- mean(sample(SM$inP, length(SM[,1]), replace=T))
  poolPR[j] <- mean(sample(SM$poolP, length(SM[,1]), replace=T))
}
round(100*quantile(inPR, prob=c(0.025,0.975)),1)
round(100*quantile(poolPR, prob=c(0.025,0.975)),1)

## differences across examination strategy
diff<- mean(SM$inP) - mean(SM$poolP)
SM$diff <- SM$inP - SM$poolP
diffSM <- rep(NA,5000)
for (j in 1:5000)
{
  t <- c(rep(1,((length(SM[,1])-1)/2)+1),rep(0,(length(SM[,1])-1)/2))
  toss <- sample(t, length(SM[,1]), replace=FALSE)
  for (i in 1:length(SM[,1]))
  {
    ifelse(toss[i] == 1, SM$diffR[i] <- SM$diff[i], SM$diffR[i] <- -SM$diff[i])
  }

  diffSM[j] <- abs(mean(SM$diffR))
}

#hist(diffSM)
sum(ifelse(diffSM > abs(diff), 1, 0))/5000


## Variation across infection intensity
quantile(SM$schi, prob=c(0.33,0.66))
SM$int <- ifelse(SM$schi<=7.2,1,ifelse(SM$schi<=15.552,2,3))
SMint<- tapply(SM$poolP,SM$int,mean)
tapply(SM$N,SM$int,sum)

## CI
poolPR1 <- rep(NA,5000)
poolPR2 <- rep(NA,5000)
poolPR3 <- rep(NA,5000)

for (j in 1:5000)
{
  poolPR1[j] <- mean(sample(SM$poolP[SM$int==1], length(SM$int[SM$int==1]), replace=T))
  poolPR2[j] <- mean(sample(SM$poolP[SM$int==2], length(SM$int[SM$int==2]), replace=T))
  poolPR3[j] <- mean(sample(SM$poolP[SM$int==3], length(SM$int[SM$int==3]), replace=T))
}
round(100*quantile(poolPR1, prob=c(0.025,0.975)),1)
round(100*quantile(poolPR2, prob=c(0.025,0.975)),1)
round(100*quantile(poolPR3, prob=c(0.025,0.975)),1)

## significant differences
SMT <- rep(NA,5000)

for (j in 1:5000)
{
  SM$t <- sample(SM$poolP,length(SM[,1]), replace = FALSE)
  d13 <- abs(mean(SM$t[SM$int==1]) - mean(SM$t[SM$int==3]))
  d12 <- abs(mean(SM$t[SM$int==1]) - mean(SM$t[SM$int==2]))
  d23 <- abs(mean(SM$t[SM$int==2]) - mean(SM$t[SM$int==3]))
  SMT[j] <- max(d13,d12,d23)
}

#hist(SMT)
diff12 <- abs(SMint[1]-SMint[2])
sum(ifelse(SMT>diff12,1,0))/5000

diff13 <- abs(SMint[1]-SMint[3])
sum(ifelse(SMT>diff13,1,0))/5000

diff23 <- abs(SMint[3]-SMint[2])
sum(ifelse(SMT>diff23,1,0))/5000

## making graph
points(33/2+5,30.8, pch = 19, col = 3)
lines(c(33/2+5,33/2+5),c(7.7,53.8), col = 3)
points(33+33/2+5,50, pch = 19, col = 3)
lines(c(33+33/2+5,33+33/2+5),c(16.7,83.3), col = 3)
points(66+33/2+5,80.0, pch = 19, col = 3)
lines(c(66+33/2+5,66+33/2+5),c(50.0,100), col = 3)
lines(c(33/2+5,33+33/2+5,66+33/2+5), c(30.8,50,80.0), col = 3)

### Sensitivity at school level
## Ascaris
AL <- subset(TOT, TOT$asc > 0 |TOT$EPG_AL_P > 0)
dim(AL)
AL$inP <- ifelse(AL$asc>0,1,0)
AL$poolP <- ifelse(AL$EPG_AL_P>0,1,0)
sum3 <- ddply(AL, .(Woreda.code,School.Code), summarize,
              ind = mean(inP),pool = mean(poolP))
dim(sum3)
sum3$inP<- ifelse(sum3$ind>0,1,0)
sum3$poolP <- ifelse(sum3$pool>0,1,0)
round(100*mean(sum3$inP),1)
round(100*mean(sum3$poolP),1)


## CI
inPR <- rep(NA,5000)
poolPR <- rep(NA,5000)
for (j in 1:5000)
{
  inPR[j] <- mean(sample(sum3$inP, length(sum3[,1]), replace=T))
  poolPR[j] <- mean(sample(sum3$poolP, length(sum3[,1]), replace=T))
}
round(100*quantile(inPR, prob=c(0.025,0.975)),1)
round(100*quantile(poolPR, prob=c(0.025,0.975)),1)

## differences across examination strategy
diff<- mean(sum3$inP) - mean(sum3$poolP)
sum3$diff <- sum3$inP - sum3$poolP
diffAL <- rep(NA,5000)
for (j in 1:5000)
{
  t <- c(rep(1,length(sum3[,1])/2),rep(0,length(sum3[,1])/2))
  toss <- sample(t, length(sum3[,1]), replace=F)
  for (i in 1:length(sum3[,1]))
  {
    ifelse(toss[i] == 1, sum3$diffR[i] <- sum3$diff[i], sum3$diffR[i] <- -sum3$diff[i])
  }

  diffAL[j] <- abs(mean(sum3$diffR))
}


hist(diffAL)
sum(ifelse(diffAL > abs(diff), 1, 0))/5000


# Trichuris
sum3 <- ddply(TT, .(Woreda.code,School.Code), summarize,
              ind = mean(inP),pool = mean(poolP))
dim(sum3)
sum3$inP<- ifelse(sum3$ind>0,1,0)
sum3$poolP <- ifelse(sum3$pool>0,1,0)
round(100*mean(sum3$inP),1)
round(100*mean(sum3$poolP),1)


## CI
inPR <- rep(NA,5000)
poolPR <- rep(NA,5000)
for (j in 1:5000)
{
  inPR[j] <- mean(sample(sum3$inP, length(sum3[,1]), replace=T))
  poolPR[j] <- mean(sample(sum3$poolP, length(sum3[,1]), replace=T))
}
round(100*quantile(inPR, prob=c(0.025,0.975)),1)
round(100*quantile(poolPR, prob=c(0.025,0.975)),1)

## differences across examination strategy
diff<- mean(sum3$inP) - mean(sum3$poolP)
sum3$diff <- sum3$inP - sum3$poolP
diffTT <- rep(NA,5000)
for (j in 1:5000)
{
  t <- c(rep(1,(length(sum3[,1])-1)/2+1),rep(0,(length(sum3[,1])/2)))
  toss <- sample(t, length(sum3[,1]), replace=F)
  for (i in 1:length(sum3[,1]))
  {
    ifelse(toss[i] == 1, sum3$diffR[i] <- sum3$diff[i], sum3$diffR[i] <- -sum3$diff[i])
  }

  diffTT[j] <- abs(mean(sum3$diffR))
}


hist(diffTT)
sum(ifelse(diffTT > abs(diff), 1, 0))/5000


# hookworm
sum3 <- ddply(HW, .(Woreda.code,School.Code), summarize,
              ind = mean(inP),pool = mean(poolP))
dim(sum3)
sum3$inP<- ifelse(sum3$ind>0,1,0)
sum3$poolP <- ifelse(sum3$pool>0,1,0)
round(100*mean(sum3$inP),1)
round(100*mean(sum3$poolP),1)


## CI
inPR <- rep(NA,5000)
poolPR <- rep(NA,5000)
for (j in 1:5000)
{
  inPR[j] <- mean(sample(sum3$inP, length(sum3[,1]), replace=T))
  poolPR[j] <- mean(sample(sum3$poolP, length(sum3[,1]), replace=T))
}
round(100*quantile(inPR, prob=c(0.025,0.975)),1)
round(100*quantile(poolPR, prob=c(0.025,0.975)),1)

## differences across examination strategy
diff<- mean(sum3$inP) - mean(sum3$poolP)
sum3$diff <- sum3$inP - sum3$poolP
diffHW <- rep(NA,5000)
for (j in 1:5000)
{
  t <- c(rep(1,(length(sum3[,1])-1)/2+1),rep(0,(length(sum3[,1])-1)/2))
  toss <- sample(t, length(sum3[,1]), replace=F)
  for (i in 1:length(sum3[,1]))
  {
    ifelse(toss[i] == 1, sum3$diffR[i] <- sum3$diff[i], sum3$diffR[i] <- -sum3$diff[i])
  }

  diffHW[j] <- abs(mean(sum3$diffR))
}


hist(diffHW)
sum(ifelse(diffHW > abs(diff), 1, 0))/5000

### Schistosoma mansoni
sum3 <- ddply(SM, .(Woreda.code,School.Code), summarize,
              ind = mean(inP),pool = mean(poolP))
dim(sum3)
sum3$inP<- ifelse(sum3$ind>0,1,0)
sum3$poolP <- ifelse(sum3$pool>0,1,0)
round(100*mean(sum3$inP),1)
round(100*mean(sum3$poolP),1)


## CI
inPR <- rep(NA,5000)
poolPR <- rep(NA,5000)
for (j in 1:5000)
{
  inPR[j] <- mean(sample(sum3$inP, length(sum3[,1]), replace=T))
  poolPR[j] <- mean(sample(sum3$poolP, length(sum3[,1]), replace=T))
}
round(100*quantile(inPR, prob=c(0.025,0.975)),1)
round(100*quantile(poolPR, prob=c(0.025,0.975)),1)

## differences across examination strategy
diff<- mean(sum3$inP) - mean(sum3$poolP)
sum3$diff <- sum3$inP - sum3$poolP
diffSM <- rep(NA,5000)
for (j in 1:5000)
{
  t <- c(rep(1,length(sum3[,1])/2),rep(0,length(sum3[,1])/2))
  toss <- sample(t, length(sum3[,1]), replace=FALSE)
  for (i in 1:length(sum3[,1]))
  {
    ifelse(toss[i] == 1, sum3$diffR[i] <- sum3$diff[i], sum3$diffR[i] <- -sum3$diff[i])
  }

  diffSM[j] <- abs(mean(sum3$diffR))
}


hist(diffSM)
sum(ifelse(diffSM > abs(diff), 1, 0))/5000

### FECs
## Scatterplots
par(mfrow=c(1,3),ps = 8, cex = 1, cex.main = 1)
TOT$AL <- TOT$asc + TOT$EPG_AL_P
TOT$TT <- TOT$tri + TOT$EPG_TT_P
TOT$HK <- TOT$hoo + TOT$EPG_HK_P

plot(c(0,9),c(0,9), type='n', xaxt='n', yaxt='n',main =  'A', xlab = 'Mean FEC of ind. samples (EPG)', ylab='FEC of pooled sample (EPG)')
points(log(TOT$asc+1),log(TOT$EPG_AL_P+1), pch=19)
axis(side = 1, at=c(log(1), log(6), log(51), log(501), log(5001)), lab=c(0,5,50,500,5000))
axis(side = 2, at=c(log(1), log(6), log(51), log(501), log(5001)), lab=c(0,5,50,500,5000))
abline(a=0,b=1, lty = 2)
alr <- cor(log(TOT$asc+1),log(TOT$EPG_AL_P+1))
text(log(500),log(5), 'R = 0.68')

plot(c(0,9),c(0,9), type='n', main =  'B', xaxt='n', yaxt='n',xlab = 'Mean FEC of ind. samples (EPG)', ylab='FEC of pooled sample (EPG)')
points(log(TOT$hoo+1),log(TOT$EPG_HK_P+1), pch=19)
axis(side = 1, at=c(log(1), log(6), log(51), log(501), log(5001)), lab=c(0,5,50,500,5000))
axis(side = 2, at=c(log(1), log(6), log(51), log(501), log(5001)), lab=c(0,5,50,500,5000))
abline(a=0,b=1, lty = 2)
hkr <- cor(log(TOT$hoo+1),log(TOT$EPG_HK_P+1))
text(log(500),log(5), 'R = 0.65')

plot(c(0,6),c(0,6), type='n', main =  'C', xaxt='n', yaxt='n', xlab = 'Mean FEC of ind. samples (EPG)', ylab='FEC of pooled sample (EPG)')
points(log(TOT$schi+1),log(TOT$EPG_SM_P+1), pch=19)
axis(side = 1, at=c(log(1), log(6), log(51), log(501), log(5001)), lab=c(0,5,50,500,5000))
axis(side = 2, at=c(log(1), log(6), log(51), log(501), log(5001)), lab=c(0,5,50,500,5000))
abline(a=0,b=1, lty = 2)
text(log(50),log(3), 'R = 0.75')
smr <- cor(log(TOT$schi+1),log(TOT$EPG_SM_P+1))

ttr <- cor(log(TOT$tri+1),log(TOT$EPG_TT_P+1))

## Level of significance
corrAL <- rep(NA,5000)
corrTT <- rep(NA,5000)
corrHK <- rep(NA,5000)
corrSM <- rep(NA,5000)

for (j in 1:5000)
{
  TOT$ALR <- sample(TOT$asc,length(TOT[,1]),replace=F)
  TOT$HKR <- sample(TOT$hoo,length(TOT[,1]),replace=F)
  TOT$TTR <- sample(TOT$tri,length(TOT[,1]),replace=F)
  TOT$SMR <- sample(TOT$schi,length(TOT[,1]),replace=F)
  corrAL[j] <- abs(cor(log(TOT$ALR+1),log(TOT$EPG_AL_P+1)))
  corrHK[j] <- abs(cor(log(TOT$HKR+1),log(TOT$EPG_HK_P+1)))
  corrTT[j] <- abs(cor(log(TOT$TTR+1),log(TOT$EPG_TT_P+1)))
  corrSM[j] <- abs(cor(log(TOT$SMR+1),log(TOT$EPG_SM_P+1)))
}

sum(ifelse(corrAL > abs(alr), 1, 0))/5000
sum(ifelse(corrTT > abs(ttr), 1, 0))/5000
sum(ifelse(corrHK > abs(hkr), 1, 0))/5000
sum(ifelse(corrSM > abs(smr), 1, 0))/5000

## Mean EPG
mean(TOT$asc)
mean(TOT$EPG_AL_P)

mean(TOT$tri)
mean(TOT$EPG_TT_P)

mean(TOT$hoo)
mean(TOT$EPG_HK_P)

mean(TOT$schi)
mean(TOT$EPG_SM_P)

## CI
TOT$diffAL <- TOT$asc-TOT$EPG_AL_P
TOT$diffTT <- TOT$tri-TOT$EPG_TT_P
TOT$diffHK <- TOT$hoo-TOT$EPG_HK_P
TOT$diffSM <- TOT$schi-TOT$EPG_SM_P

alind <- rep(NA,5000)
alpool <- rep(NA,5000)
diffal <- rep(NA,5000)
ttind <- rep(NA,5000)
ttpool <- rep(NA,5000)
difftt <- rep(NA,5000)
hkind <- rep(NA,5000)
hkpool <- rep(NA,5000)
diffhk <- rep(NA,5000)
smind <- rep(NA,5000)
smpool <- rep(NA,5000)
diffsm <- rep(NA,5000)

for (j in 1:5000)
{
  alind[j] <- mean(sample(TOT$asc,length(TOT[,1]),replace=T))
  alpool[j] <- mean(sample(TOT$EPG_AL,length(TOT[,1]),replace=T))
  diffal[j] <- mean(sample(TOT$diffAL,length(TOT[,1]),replace=T))

  ttind[j] <- mean(sample(TOT$tri,length(TOT[,1]),replace=T))
  ttpool[j] <- mean(sample(TOT$EPG_TT,length(TOT[,1]),replace=T))
  difftt[j] <- mean(sample(TOT$diffTT,length(TOT[,1]),replace=T))

  hkind[j] <- mean(sample(TOT$hoo,length(TOT[,1]),replace=T))
  hkpool[j] <- mean(sample(TOT$EPG_HK,length(TOT[,1]),replace=T))
  diffhk[j] <- mean(sample(TOT$diffHK,length(TOT[,1]),replace=T))

  smind[j] <- mean(sample(TOT$schi,length(TOT[,1]),replace=T))
  smpool[j] <- mean(sample(TOT$EPG_SM,length(TOT[,1]),replace=T))
  diffsm[j] <- mean(sample(TOT$diffSM,length(TOT[,1]),replace=T))
}
quantile(alind, prob=c(0.025,0.975))
quantile(alpool, prob=c(0.025,0.975))
quantile(diffal, prob=c(0.025,0.975))

quantile(ttind, prob=c(0.025,0.975))
quantile(ttpool, prob=c(0.025,0.975))
quantile(difftt, prob=c(0.025,0.975))

quantile(hkind, prob=c(0.025,0.975))
quantile(hkpool, prob=c(0.025,0.975))
quantile(diffhk, prob=c(0.025,0.975))

quantile(smind, prob=c(0.025,0.975))
quantile(smpool, prob=c(0.025,0.975))
quantile(diffsm, prob=c(0.025,0.975))


## Difference
diffal <- rep(NA,5000)
difftt <- rep(NA,5000)
diffhk <- rep(NA,5000)
diffsm <- rep(NA,5000)

for (j in 1:5000)
{
  t <- c(rep(1,(length(TOT[,1])-1)/2+1),rep(0,(length(TOT[,1])-1)/2))
  toss <- sample(t, length(TOT[,1]), replace=F)
  for (i in 1:length(TOT[,1]))
  {
    ifelse(toss[i] == 1, TOT$diffalR[i] <- TOT$diffAL[i], TOT$diffalR[i] <- -TOT$diffAL[i])
    ifelse(toss[i] == 1, TOT$diffttR[i] <- TOT$diffTT[i], TOT$diffttR[i] <- -TOT$diffTT[i])
    ifelse(toss[i] == 1, TOT$diffhkR[i] <- TOT$diffHK[i], TOT$diffhkR[i] <- -TOT$diffHK[i])
    ifelse(toss[i] == 1, TOT$diffsmR[i] <- TOT$diffSM[i], TOT$diffsmR[i] <- -TOT$diffSM[i])
  }

  diffal[j] <- abs(mean(TOT$diffalR))
  difftt[j] <- abs(mean(TOT$diffttR))
  diffhk[j] <- abs(mean(TOT$diffhkR))
  diffsm[j] <- abs(mean(TOT$diffsmR))
}

sum(ifelse(diffal>abs(mean(TOT$diffAL)), 1,0))/5000
sum(ifelse(difftt>abs(mean(TOT$diffTT)), 1,0))/5000
sum(ifelse(diffhk>abs(mean(TOT$diffHK)), 1,0))/5000
sum(ifelse(diffsm>abs(mean(TOT$diffSM)), 1,0))/5000

### Time to process samples
pool[1:10,]
ind[1:10,]
prepkk[1:10,]
prepool[1:10,]


### Preperatio Kato-Katz
prepkk$tot <- prepkk$Preparation.time..MIN.*60 + prepkk$Preparation.time..SEC.
mean(prepkk$tot)/60
sqrt(var(prepkk$tot))/60


### Tornado plots
par(mfrow=c(2,3),ps = 8, cex = 1, cex.main = 1)

## Impact on total cost
## Individual
# Poor accessibility
tplot <- read.csv('/Users/blevecke/Dropbox/Communications/Publicaties/Drafts/Pooling/Human/Ethiopia/Gemechu/Let et al., 2017/tornplots.csv', header = TRUE, sep = ";")
tplot[1:10,]
barplot(tplot$Pdiff10.[tplot$Acces == 'poor'],xlab='Relative change in total costs (%)',main='Poor accessibility',names.arg=c('Data entry','Material','Fee teacher','Fuel', 'Salary','Rent car'), horiz = T, las=1, xlim = c(-5,5), ylab = '', col=c('white'))
barplot(tplot$Mdiff10.[tplot$Acces == 'poor'], las = 1,horiz = T, xlim = c(-6,6), ylab ='',beside=T, col=c('black'), add = TRUE)

## Moderate accessibility
tplot <- read.csv('/Users/blevecke/Desktop/tornplots.csv', header = TRUE, sep = ";")
tplot[1:10,]
barplot(tplot$Pdiff10.[tplot$Acces == 'Moderate'],xlab='Relative change in total costs (%)',main='Moderate accessibility',names.arg=c('Data entry','Material','Fee teacher','Fuel', 'Salary','Rent car'), horiz = T, las=1, xlim = c(-5,5), ylab = '', col=c('white'))
barplot(tplot$Mdiff10.[tplot$Acces == 'Moderate'], las = 1,horiz = T, xlim = c(-6,6), ylab ='',beside=T, col=c('black'), add = TRUE)

## High accessibility
tplot <- read.csv('/Users/blevecke/Desktop/tornplots.csv', header = TRUE, sep = ";")
tplot[1:10,]
barplot(tplot$Pdiff10.[tplot$Acces == 'High'],xlab='Relative change in total costs (%)',main='High accessibility',names.arg=c('Data entry','Material','Fee teacher','Fuel', 'Salary','Rent car'), horiz = T, las=1, xlim = c(-5,5), ylab = '', col=c('white'))
barplot(tplot$Mdiff10.[tplot$Acces == 'High'], las = 1,horiz = T, xlim = c(-6,6), ylab ='',beside=T, col=c('black'), add = TRUE)


## Pooled
# Poor accessibility
tplot <- read.csv('/Users/blevecke/Desktop/tornplots.csv', header = TRUE, sep = ";")
tplot[1:10,]
barplot(tplot$Pdiff10._P[tplot$Acces == 'poor'],xlab='Relative change in total costs (%)',main='',names.arg=c('Data entry','Material','Fee teacher','Fuel', 'Salary','Rent car'), horiz = T, las=1, xlim = c(-5,5), ylab = '', col=c('white'))
barplot(tplot$Mdiff10._P[tplot$Acces == 'poor'], las = 1,horiz = T, xlim = c(-6,6), ylab ='',beside=T, col=c('black'), add = TRUE)

## Moderate accessibility
tplot <- read.csv('/Users/blevecke/Desktop/tornplots.csv', header = TRUE, sep = ";")
tplot[1:10,]
barplot(tplot$Pdiff10._P[tplot$Acces == 'Moderate'],xlab='Relative change in total costs (%)',main='',names.arg=c('Data entry','Material','Fee teacher','Fuel', 'Salary','Rent car'), horiz = T, las=1, xlim = c(-5,5), ylab = '', col=c('white'))
barplot(tplot$Mdiff10._P[tplot$Acces == 'Moderate'], las = 1,horiz = T, xlim = c(-6,6), ylab ='',beside=T, col=c('black'), add = TRUE)

## High accessibility
tplot <- read.csv('/Users/blevecke/Desktop/tornplots.csv', header = TRUE, sep = ";")
tplot[1:10,]
barplot(tplot$Pdiff10._P[tplot$Acces == 'High'],main='',xlab='Relative change in total costs (%)',names.arg=c('Data entry','Material','Fee teacher','Fuel', 'Salary','Rent car'), horiz = T, las=1, xlim = c(-5,5), ylab = '', col=c('white'))
barplot(tplot$Mdiff10._P[tplot$Acces == 'High'], las = 1,horiz = T, xlim = c(-6,6), ylab ='',beside=T, col=c('black'), add = TRUE)



### Impact on reduction
par(mfrow=c(1,3),ps = 8, cex = 1, cex.main = 1)

# Poor accessibility
#tplot <- read.csv('/Users/blevecke/Desktop/tornplots.csv', header = TRUE, sep = ";")
tplot[1:10,]
tplot$diffRL <- tplot$ReferenceR-tplot$X.10.D
tplot$diffRU <- tplot$ReferenceR-tplot$X10.D
tplot <- tplot[order(tplot$Acces,abs(tplot$diffRL)),]

barplot(tplot$diffRL[tplot$Acces == 'poor'],xlab='Point percent difference in cost savings',main='Poor accessibility',names.arg=c('Data entry','Material','Fee teacher','Fuel', 'Rent car', 'Salary'), horiz = T, las=1, xlim = c(-1,1), ylab = '', col=c('white'))
barplot(tplot$diffRU[tplot$Acces == 'poor'], las = 1,horiz = T, xlim = c(-1,1), ylab ='',beside=T, col=c('black'), add = TRUE)

## Moderate accessibility
barplot(tplot$diffRL[tplot$Acces == 'Moderate'],xlab='Point percent difference in cost savings',main='Moderate accessibility', horiz = T, las=1, xlim = c(-1,1), ylab = '', col=c('white'))
barplot(tplot$diffRU[tplot$Acces == 'Moderate'], las = 1,horiz = T, xlim = c(-1,1), ylab ='',beside=T, col=c('black'), add = TRUE)

## High accessibility
barplot(tplot$diffRL[tplot$Acces == 'High'],xlab='Point percent difference in cost savings',main='High accessibility', horiz = T, las=1, xlim = c(-1,1), ylab = '', col=c('white'))
barplot(tplot$diffRU[tplot$Acces == 'High'], las = 1,horiz = T, xlim = c(-1,1), ylab ='',beside=T, col=c('black'), add = TRUE)

