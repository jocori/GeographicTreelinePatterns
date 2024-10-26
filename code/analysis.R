#set working directory
setwd("~/Desktop/KU/Projects/GeographicTreelinePatterns")
#load in packages
library(dplyr)
# read in data
regs<-read.csv("data/Regressions_3Jan.csv", header = TRUE)
count(regs[regs$Significance >= 0.05,])

#Keep linear models only
regs<-regs[regs$Regression_Type == "Linear",]

#remove extraneous columns
regs<-regs[,!(names(regs)%in%c("Regression_Type","X"))]

##percent linear models that are significant
count(regs[regs$Significance >= 0.05,])
## should we only keep significant models?
## should I subtract 8488 from 1317? do a range from 8488 to 9802 and then 1317? Calculate rate?
## what should be included in final model? Final model type? Any random effect?

##slope is ndvi~elevation...so the relationship between vegetation density and elevation

##what should the response variable be? 

##############################
#keep best model for each:
#regs<-regs %>%
 # filter(regs$Regression_Type == regs$Best_Model)
#unique(regs$Regression_Type)
#dim(regs[regs$Regression_Type == "Reciprocal_Quadratic",])

### solve y = mx + b
###midpoint btw min and max NDVI, plug into equation and see what elevation
### range divide by 2 add to minimum
###elevation of treeline
### difference betweeen present and 80s
