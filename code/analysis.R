#set working directory
setwd("~/Desktop/KU/Projects/GeographicTreelinePatterns")
#load in packages
library(dplyr)
# read in data
regs<-read.csv("data/Regressions_3Jan.csv", header = TRUE)

#Keep linear models only
regs<-regs[regs$Regression_Type == "Linear",]

#remove extraneous columns
regs<-regs[,!(names(regs)%in%c("Regression_Type","X"))]

## should we only keep significant models?
## should I subtract 8488 from 1317? do a range from 8488 to 9802 and then 1317? Calculate rate?
## what should be included in final model? Final model type? Any random effect?

##slope is ndvi~elevation...so the relationship between vegetation density and elevation

##what should the response variable be? 