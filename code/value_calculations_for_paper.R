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
#with a modal value of AAAA m.
max(regs$change_in_treeline_elevation)

#pearson correlation between distance to coast and latitude
cor.test(regs$Long, regs$dist_coast, method="pearson")
cor.test( regs$dist_coast,regs$change_in_treeline_elevation, method="pearson")
cor.test(regs$change_in_treeline_NDVI, regs$dist_coast, method="pearson")
