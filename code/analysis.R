#set working directory
setwd("~/Desktop/KU/Projects/GeographicTreelinePatterns")
#load in packages
library(dplyr)
library(lme4)
library(AICcmodavg)
# read in data
regs_full<-read.csv("data/Regressions_28Oct.csv", header = TRUE)
count(regs_full[regs_full$Significance >= 0.05,])

#Keep linear models only
regs<-regs_full[regs_full$Regression_Type == "Linear",]

#Keep only significant models
regs<-regs[regs$Significance<=0.05,]

#Keep only models with negative slope
regs<-regs[regs$Slope < 0,]

#Keep only average NDVI NOT Max (relevant in year calculation)
regs <- regs %>% filter(Year %in% c("Avg8488", "Avg9802", "Avg1317"))

#Which mountains are excluded now?
x1<-unique(regs$Peak)
x2<-unique(regs_full$Peak)
# no mountains are completely excluded! 

#remove extraneous columns
regs<-regs[,!(names(regs)%in%c("Regression_Type","X"))]


## what should be included in final model? Final model type? Any random effect?



##what should the response variable be? 


### solve y = mx + b
# NDVI = slope(elevation) + intercept....
###midpoint btw min and max NDVI, plug into equation and see what elevation
### range divide by 2 add to minimum
###elevation of treeline
### difference betweeen present and 80s

## calculating response variable

#elevation at treeline = (NDVI- Intercept)/Slope
#regs$elevation_treeline<-(regs$treeline_NDVI-regs$Intercept)/regs$Slope
#exclude rows where elevation at treeline is negative...only 8/1748
#regs<-regs[regs$elevation_treeline >0,]
head(regs)

#now...subtract 80s from 2017 to get final response variable
# Ensure that rows are grouped by Peak and Direction so that calculations are done within these groups
#y = mx+ b...NDVI = slope*treeline_elevation + intercept
regs$treeline_NDVI <- (regs$Slope*regs$Treeline_elevation) + regs$Intercept
# Filter out peak-direction combinations that do not have both years
regs_filtered <- regs %>%
  group_by(Peak, Direction) %>%
  filter(all(c("Avg8488", "Avg1317") %in% Year)) %>%
  ungroup()

# Calculate the change in treeline elevation
regs <- regs_filtered %>%
  group_by(Peak, Direction) %>%
  mutate(
    change_in_treeline_NDVI = treeline_NDVI[Year == "Avg1317"] - treeline_NDVI[Year == "Avg8488"]
  ) %>%
  ungroup()

# View the updated dataset
head(regs)

#remove unnecessary columns
data<-regs[,-c(3:8,14,15)]
# Remove duplicate rows based on all columns
data <- data %>%
  distinct()
#Which mountains are excluded now?
x1<-unique(data$Peak)
x2<-unique(regs_full$Peak)
setdiff(x2, x1)

###### mountains excluded from analysis ######
## "Mt. Harrison"   "Mt. Ovington"   "Mt. Timpanogos" "Star Peak"  
head(data)

## model selection
hist(data$change_in_treeline_NDVI)
plot(data$change_in_treeline_NDVI~data$Lat, data = data)
plot(data$change_in_treeline_NDVI~data$Long, data = data)


m1<-lmer(change_in_treeline_NDVI~Lat +Long + (1|Peak_ID), data = data)
summary(m1)
AIC(m1)
m2<-lmer(change_in_treeline_NDVI~Lat *Long + (1|Peak_ID), data = data)
summary(m2)
AIC(m2)
m3<-lmer(change_in_treeline_NDVI~Lat +Long + Stations_After_Treeline +(1|Peak_ID), data = data)
summary(m3)
AIC(m3)
m4<-lmer(change_in_treeline_NDVI~  Stations_After_Treeline +Lat*Long+(1|Peak_ID), data = data)
summary(m4)
AIC(m4)
m5<-lmer(change_in_treeline_NDVI~Lat +Long + Direction +(1|Peak_ID), data = data)
summary(m5)
AIC(m5)
m6<-lmer(change_in_treeline_NDVI~ Direction +Lat*Long+(1|Peak_ID), data = data)
summary(m6)
AIC(m6)
m7<-lmer(change_in_treeline_NDVI~Lat +Long + Stations_After_Treeline +Direction+(1|Peak_ID), data = data)
summary(m7)
AIC(m7)
m8<-lmer(change_in_treeline_NDVI~Stations_After_Treeline +Direction+Lat*Long+(1|Peak_ID), data = data)
summary(m8)
AIC(m8)

mods<-list(m1,m2,m3,m4,m5,m6,m7,m8)
aictab(cand.set = mods)

##Model 1 is the best model
m1<-lmer(change_in_treeline_NDVI~Lat +Long + (1|Peak_ID), data = data)
summary(m1)
confint(m1)

###extract elevation from treesbegin  for present day keyed to NDVI from y = mx + b  , use that elevatino to get NDVI at treeline in the 1980s

#### y axis should be change in elevation at treeline.....to get that, 
#elevation at trees begin -> present day elevation at treeline, 
#put in present day regression equation > 
#get out of that the NDVI at present at treesbegin elevation > 
#plug that NDVI into 1980s relationship from regressions y = mx +b > 
#solve for elevation at that NDVI value

#split present day and 1980s into seperate dataframes 
present<-regs_full[regs_full$Year == "Avg1317",]
past<-regs_full[regs_full$Year == "Avg8488",]
regs_full$Treeline_elevation

#present day regression equation
present$present_NDVI<-(present$Treeline_elevation*present$Slope) + present$Intercept

#plug present day NDVI at treesbegin value into past regression equation to calculate past elevation
#y = mx + b....x = y-b/m...(present ndvi - intercept)/slope
past$past_elevation <- (present$present_NDVI -past$Intercept)/past$Slope
colnames(present$Treeline_elevation) <- "present_elevation" # fix column name
#get difference