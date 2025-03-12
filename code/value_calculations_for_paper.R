### calculating summary values for results
#set working directory
setwd("~/Desktop/KU/Projects/GeographicTreelinePatterns")
#read in data
regs<-read.csv("data/regs_final22jan25.csv")
head(regs)

#On average, across all mountain peaks in our analyses, treeline was located at XXXX m. 
mean(regs$present_elevation)
#tropical treelines averaged XXXX m, 
#tropic of cancer is 23.5 degrees n
trop<-regs[regs$Lat <=23.5,]
mean(trop$present_elevation)
#whereas temperate-zone treelines were lower, at YYYY m.
temp<-regs[regs$Lat >23.5,]
mean(temp$present_elevation)

#treeline shifts ranged from upward expansion of XXXX m to 
#2017 - 1984, if 2017 is bigger, positive...upward shift
pos<-regs[regs$change_in_treeline_elevation>=0,]
max(pos$change_in_treeline_elevation)
mean(pos$change_in_treeline_elevation)
#downward shifts of YYYY m. (negative values of change in treeline elevation)
neg<-regs[regs$change_in_treeline_elevation<0,]
min(neg$change_in_treeline_elevation)
#The average change was ZZZZ m, 
mean(regs$change_in_treeline_elevation)

#Average with absolute values only is 240 m
mean(abs(regs$change_in_treeline_elevation))

#with a modal value of AAAA m.
max(regs$change_in_treeline_elevation)

#pearson correlation between distance to coast and latitude
cor.test(regs$Long, regs$dist_coast, method="pearson")
cor.test( regs$dist_coast,regs$change_in_treeline_elevation, method="pearson")
cor.test(regs$change_in_treeline_NDVI, regs$dist_coast, method="pearson")

##past and present treeline NDVI for Malinche
#ndvi = elevation*slope+intercept
elevation_pres<-regs$Treeline_elevation[regs$Peak=="Malinche"&regs$Direction=="N"&regs$Year == "Avg1317"]
#present elevation = 3960
slope_pres<-regs$Slope[regs$Peak=="Malinche"&regs$Direction=="N"&regs$Year == "Avg1317"]
#slope = -0.001058405
intercept_pres<-regs$Intercept[regs$Peak=="Malinche"&regs$Direction=="N"&regs$Year == "Avg1317"]
#intercept = 4.504737
ndvi_pres<-(elevation_pres*slope_pres) +intercept_pres
#ndvi = 0.3134539
slope_past<-regs$Slope[regs$Peak=="Malinche"&regs$Direction=="N"&regs$Year == "Avg8488"]
intercept_past<-regs$Intercept[regs$Peak=="Malinche"&regs$Direction=="N"&regs$Year == "Avg8488"]
#past elevation = (present ndvi-past intercept)/past slope
past_elevation<-(ndvi_pres-intercept_past)/slope_past
elevation_pres-past_elevation


###troubleshooting why extreme values
regs<-read.csv("testing_data_includes_slope")

subset_data <- regs %>% 
  filter(present_NDVI > Intercept, Year == "Avg8488")
print(subset_data)

# Calculate the difference between present_NDVI and the 1980s Intercept
past<-regs[regs$Year =="Avg8488",]
present<-regs[regs$Year == "Avg1317",]
dNDVI <- present$present_NDVI - past$Intercept

# Summarize the differences
summary(regs$dNDVI)
summary(past$Slope)
# Plot to visualize the differences
library(ggplot2)
ggplot(regs, aes(x = dNDVI, y = Slope)) +
  geom_point() +
  labs(x = "present_NDVI - 1980s Intercept", y = "Slope (1980s)", 
       title = "Difference between present_NDVI and 1980s Intercept vs Slope")

##calculating median and wuartiles as it is a better estimate of the treeline movement
median(regs$present_elevation[regs$Year == "Avg1317"]) #2114
median(regs$past_elevation[regs$Year == "Avg8488"]) #2210.469

median(regs$treeline_ndvi[regs$Year == "Avg1317"]) #0.1030599
median(regs$treeline_ndvi[regs$Year == "Avg8488"]) #0.1349535

summary(regs$present_elevation[regs$Year == "Avg1317"]) #quartiles: 1754, 3278
summary(regs$past_elevation[regs$Year == "Avg8488"]) # quartiles: 1806.8, 3246.4
summary(regs$treeline_ndvi[regs$Year == "Avg1317"]) #quartiles: 0.0640, 0.1714
summary(regs$treeline_ndvi[regs$Year == "Avg8488"]) #quartiles: 0.08467, 0.19705

median(regs$change_in_treeline_elevation) #-36.12678
summary(regs$change_in_treeline_elevation) #quartiles: -165.09, 126.61

median(regs$change_in_treeline_NDVI) #-0.009955125
summary(regs$change_in_treeline_NDVI) #quartiles: -0.055891, 0.029747

