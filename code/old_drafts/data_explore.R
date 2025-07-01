######data exploration
setwd("C:/Users/j364g325/Documents/mountain peaks")
regs<-read.csv("peak_regressions_atLeastLast10Stations.csv")
plot(regs$Slope,regs$Significance)
head(regs)
range(na.omit(regs$Significance))
regs[regs$Significance >=0.05,]
lin<-regs[regs$Regression_Type == "Linear",]
non.sig<-lin[lin$Significance >=0.05,]
