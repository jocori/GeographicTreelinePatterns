#set working directory
setwd("~/Desktop/KU/Projects/GeographicTreelinePatterns")
#load in packages
library(dplyr)
library(lme4)
library(AICcmodavg)
# read in data
regs_full<-read.csv("data/Regressions_28Oct.csv", header = TRUE)
count(regs_full[regs_full$Significance >= 0.05,])
#remove mt. Washington
regs <-regs_full[regs_full$Peak != "Mt. Washington",]
#Keep linear models only
regs<-regs[regs$Regression_Type == "Linear",]
unique(regs$Best_Model)
regs_full[regs_full$Best_Model=="Reciprocal_Quadratic",]

#Keep only average NDVI NOT Max (relevant in year calculation)
regs <- regs %>% filter(Year %in% c("Avg8488", "Avg1317"))
#Keep only significant models
regs<-regs[regs$Significance<=0.05,]

#Keep only models with negative slope
regs<-regs[regs$Slope < 0,]

# Filter out peak-direction combinations that do not have both years
regs <- regs %>%
  group_by(Peak, Direction) %>%
  filter(all(c("Avg8488", "Avg1317") %in% Year)) %>%
  ungroup()
unique(regs$Peak)
#remove extraneous columns
regs<-regs[,!(names(regs)%in%c("Regression_Type","X","Best_Model"))]
####
##All models with change in NDVI as response variable
##have treeline_elevation...this is present day
## equations are ndvi~elevation...so ndvi = m(elevation) + b
## present day NDVI = m (elevation 2017) + b
## solve for past elevation....NDVI from present day calculation = m (elevation past) + b


####NDVI at 2017 treeline elevation
#separate past and present

# solve for 80s and 2017 ndvi using 2017 treeline elevation
regs$treeline_ndvi<-(regs$Slope*regs$Treeline_elevation) + regs$Intercept


# Calculate the change in treeline elevation
regs <- regs %>%
  group_by(Peak, Direction) %>%
  mutate(
    change_in_treeline_NDVI = treeline_ndvi[Year == "Avg1317"] - treeline_ndvi[Year == "Avg8488"]
  ) %>%
  ungroup()

# NDVI should get higher...point should be bare rock in 1980s and green in 2017

##other response variable...change in treeline elevation
past<-regs[regs$Year == "Avg8488",]
present<-regs[regs$Year == "Avg1317",]
#present day regression equation
present$present_NDVI<-(present$Treeline_elevation*present$Slope) + present$Intercept

#plug present day NDVI at treesbegin value into past regression equation to calculate past elevation
#y = mx + b....x = y-b/m...(present ndvi - intercept)/slope
past$past_elevation <- (present$present_NDVI -past$Intercept)/past$Slope
cols<-names(present)
colnames(present)[13] <- "present_elevation" # fix column name
colnames(past)[13] <- "present_elevation" # fix column name

#merge past and present
regs <- merge(present, past, by = intersect(names(present),names(past)), all.y = TRUE, all.x = TRUE)
#remove unnecessary columns
regs<-regs[,-c(4:7)]

# Calculate the change in treeline elevation
regs <- regs %>%
  group_by(Peak, Direction) %>%
  mutate(
    change_in_treeline_elevation = present_elevation[Year == "Avg1317"] - past_elevation[Year == "Avg8488"]
  ) %>%
  ungroup()

#remove present elevation, past elevation, present NDVI, and year columns
regs<-regs[,-c(3,10,12:13)]

#remove duplicate rows
regs <- regs %>%
  distinct()

#Which mountains are excluded now?
x1<-unique(regs$Peak)
x2<-unique(regs_full$Peak)
setdiff(x2, x1) #"Mt. Harrison"   "Mt. Ovington"   "Mt. Timpanogos" "Star Peak"  

###### mountains excluded from analysis ######
## "Mt. Harrison"   "Mt. Ovington"   "Mt. Timpanogos" "Star Peak"  

unique(regs$Peak)
unique(regs_full$Peak)

write.csv(regs, "data/regsJan2025.csv")
