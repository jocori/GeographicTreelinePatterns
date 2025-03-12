#create data for distance to coast calculation in qgis
setwd("~/Desktop/KU/Projects/GeographicTreelinePatterns")
# read in data
regs_full<-read.csv("data/Regressions_28Oct.csv", header = TRUE)

regs <- regs_full[regs_full$Year %in% c("Avg8488", "Avg1317"),]

unique(regs$Year)
#Keep linear models only
regs<-regs[regs$Regression_Type == "Linear",]
unique(regs$Year) #gets rid of 1317 here....why????

#Keep only significant models
regs<-regs[regs$Significance<=0.05,]

#Keep only models with negative slope
regs<-regs[regs$Slope < 0,]

unique(regs$Year)

#Which mountains are excluded now?
x1<-unique(regs$Peak)
x2<-unique(regs_full$Peak)
setdiff(x1,x2)
# no mountains are completely excluded! 

#remove extraneous columns
regs<-regs[,!(names(regs)%in%c("Regression_Type","X"))]
head(regs)
# Filter out peak-direction combinations that do not have both years
regs<-regs %>%
  group_by(Peak, Direction) %>%
  filter(n_distinct(Year) == 2 & all(c("Avg8488", "Avg1317") %in% Year)) %>%
  ungroup()

present<-regs[regs$Year == "Avg1317",]
past<-regs[regs$Year == "Avg8488",]

#present day regression equation
present$present_NDVI<-(present$Treeline_elevation*present$Slope) + present$Intercept
unique(regs$Peak)


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

coords<-regs[,-c(1:4,7,8)]
write.csv(coords, "data/coords.csv")

#create regs including distance to coast
coast<-read.csv("data/coastdist_peaks.csv")
colnames(coast)[5] <- "dist_coast"
coast<- coast[,-c(1:2)]
library(dplyr)
dist_coast<- read.csv("data/Regressions12Nov24.csv")

# Add only the 'dist_coast' column to 'regs'
regs <- regs %>%
  left_join(dist_coast %>% select(Peak, Direction, dist_coast), by = c("Peak", "Direction"))

regs <- regs %>%
  left_join(coast, by = c("Lat" = "Lat", "Long" = "Long"))

write.csv(regs, "data/Regressions12Nov24.csv")
