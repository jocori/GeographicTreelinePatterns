#set working directory
setwd("~/Desktop/KU/Projects/GeographicTreelinePatterns")
#load in packages
library(dplyr)
library(lme4)
library(AICcmodavg)
library(ggplot2)
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




###########################################
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

# Filter out peak-direction combinations that do not have both years
regs <- regs %>%
  group_by(Peak, Direction) %>%
  filter(all(c("Avg8488", "Avg1317") %in% Year)) %>%
  ungroup()

present<-regs[regs$Year == "Avg1317",]
past<-regs[regs$Year == "Avg8488",]

#present day regression equation
present$present_NDVI<-(present$Treeline_elevation*present$Slope) + present$Intercept

#plug present day NDVI at treesbegin value into past regression equation to calculate past elevation
#y = mx + b....x = y-b/m...(present ndvi - intercept)/slope
past$past_elevation <- (present$present_NDVI -past$Intercept)/past$Slope
cols<-names(present)
colnames(present)[14] <- "present_elevation" # fix column name

#merge past and present
regs <- merge(present, past, by = intersect(names(present),names(past)), all.y = TRUE, all.x = TRUE)
#remove unnecessary columns
regs<-regs[,-c(4:8,15:16)]

# Calculate the change in treeline elevation
regs <- regs %>%
  group_by(Peak, Direction) %>%
  mutate(
    change_in_treeline_elevation = present_elevation[Year == "Avg1317"] - past_elevation[Year == "Avg8488"]
  ) %>%
  ungroup()

#remove present and past elevation columns
regs<-regs[,-c(9:10)]

#remove year value
regs<-regs[,-c(3)]
#remove duplicate rows
regs <- regs %>%
  distinct()

#Which mountains are excluded now?
x1<-unique(regs$Peak)
x2<-unique(regs_full$Peak)
setdiff(x2, x1) #"Mt. Harrison"   "Mt. Ovington"   "Mt. Timpanogos" "Star Peak"  

###### mountains excluded from analysis ######
## "Mt. Harrison"   "Mt. Ovington"   "Mt. Timpanogos" "Star Peak"  


## model selection
hist(regs$change_in_treeline_elevation)
plot(regs$change_in_treeline_elevation~regs$Lat, data = regs)
plot(regs$change_in_treeline_elevation~regs$Long, data = regs)


m1<-lmer(change_in_treeline_elevation~Lat +Long + (1|Peak_ID), data = regs)
summary(m1)
AIC(m1)
m2<-lmer(change_in_treeline_elevation~Lat *Long + (1|Peak_ID), data = regs)
summary(m2)
AIC(m2)
m3<-lmer(change_in_treeline_elevation~Lat +Long + Stations_After_Treeline +(1|Peak_ID), data = regs)
summary(m3)
AIC(m3)
m4<-lmer(change_in_treeline_elevation~  Stations_After_Treeline +Lat*Long+(1|Peak_ID), data = regs)
summary(m4)
AIC(m4)
m5<-lmer(change_in_treeline_elevation~Lat +Long + Direction +(1|Peak_ID), data = regs)
summary(m5)
AIC(m5)
m6<-lmer(change_in_treeline_elevation~ Direction +Lat*Long+(1|Peak_ID), data = regs)
summary(m6)
AIC(m6)
m7<-lmer(change_in_treeline_elevation~Lat +Long + Stations_After_Treeline +Direction+(1|Peak_ID), data = regs)
summary(m7)
AIC(m7)
m8<-lmer(change_in_treeline_elevation~Stations_After_Treeline +Direction+Lat*Long+(1|Peak_ID), data = regs)
summary(m8)
AIC(m8)

mods<-list(m1,m2,m3,m4,m5,m6,m7,m8)
aictab(cand.set = mods)

##Model 8 is the best model
m8<-lmer(change_in_treeline_elevation~Stations_After_Treeline +Direction+Lat*Long+(1|Peak_ID), data = regs)
summary(m8)
confint(m8)

regs_plot<-regs[,-c(13:24)]
regs_plot %>%
  ggplot(aes(y = change_in_treeline_elevation, x = Lat)) + 
  geom_smooth() + 
  xlab("Latitude") + 
  ylab("Change in Treeline Elevation (2017 - 1984)")

regs_plot %>%
  ggplot(aes(y = change_in_treeline_elevation, x = Long)) + 
  geom_smooth() + 
  xlab("Longitude") + 
  ylab("Change in Treeline Elevation (2017 - 1984)")

regs_plot %>%
  ggplot(aes(x = Lat, y = change_in_treeline_elevation, color = Long)) + 
  geom_point(size = 3, alpha = 0.7) + 
  scale_color_gradient(name = "Longitude", low = "blue", high = "red") + 
  xlab("Latitude") + 
  ylab("Change in Treeline Elevation (2017 - 1984)") + theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = 
          element_line(colour = "black", linewidth = 0.25), 
        axis.text = element_text(family = "sans"), 
        axis.ticks = element_line(),
        axis.title = element_text(size = 13, family = "sans"))

regs_full %>%
  filter(Year %in% c("Avg8488", "Avg1317")) %>% 
  mutate(Region = ifelse(Lat > 23.5, "Tropical Mountains", "Temperate Mountains")) %>%
  ggplot(aes(y = Treeline_elevation, x = Year)) + geom_boxplot() + ylim(1900,2200)+
  facet_wrap(~ Region)

regs_full %>%
  filter(Year %in% c("Avg8488", "Avg1317")) %>%
  mutate(Region = ifelse(Lat > 23.5, "Tropical Mountains", "Temperate Mountains")) %>%
  ggplot(aes(x = Treeline_elevation, fill = Year)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Region) +
  xlab("Treeline Elevation") +
  ylab("Density") +
  ggtitle("Distribution of Treeline Elevations by Region and Year")

regs_full %>%
  filter(Year %in% c("Avg8488", "Avg1317")) %>%
  mutate(Region = ifelse(Lat > 23.5, "Tropical Mountains", "Temperate Mountains")) %>%
  group_by(Region) %>%
  summarize(diff_elevation = diff(Treeline_elevation[order(Year)]),
            .groups = "drop") %>%
  ggplot(aes(x = Region, y = diff_elevation, fill = Region)) +
  geom_bar(stat = "identity") +
  ylab("Change in Treeline Elevation (2017 - 1984)") +
  ggtitle("Change in Treeline Elevation by Region")


#should we throw in distance to coast as another predictor??

### genetics and molecular biology 32:203 Diniz-Filho et al 2009...don't necessarily do this, just check it out





