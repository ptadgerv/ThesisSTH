expand.grid
library(dplyr)
tablita<-expand.grid(Day = c(2, 3),Sample = c(1,2,3),NxSchool = 50,TP = c(1,2,3),LikeR = c(1.5,3))
tablita$Scenario<-rownames(tablita)
tablita<-tablita[,c(6,1:5)]
tablita%>%flextable::flextable()
