#####plots
setwd("C:/Users/j364g325/Documents/mountain peaks/data")
regressions<-read.csv("Regressions_8Dec.csv") #read in data
sig.88.lin.neg<-regressions[regressions$Significance <= 0.05 & regressions$Year == "Avg8488"&
              regressions$Regression_Type == "Linear"
              & regressions$Slope <= 0,] #subset to significant linear models with negative slopes from 1984-1988 only
#How many models show what we expect?
x<-regressions[regressions$Significance <= 0.05&
              regressions$Regression_Type == "Linear"
            & regressions$Slope <= 0,]
plot(Slope~Lat,data = regressions,
     ylim = c(-1,1)) #74530/108378 linear models show what we expect
74530/108378 #68% of models


#Create a plot of one model to show example: Abercrombie Mountain
#subset to significant linear models with negative slopes from 1984-1988 only
sig.88.lin.neg<-regressions[regressions$Significance <= 0.05 
                            & regressions$Year == "Avg8488"&
                              regressions$Regression_Type == "Linear"
                            & regressions$Slope <= 0,] 
#get slope and intercept from 2013-2017
#subset to significant linear models with negative slopes from 2013-2017 only
ab.17<-regressions[regressions$Peak == "Abercrombie Mtn." & regressions$Year == "Avg1317"&
                              regressions$Regression_Type == "Linear"
                            & regressions$Direction == "E",] 
intercept<-sig.88.lin.neg[1,7] 
slope<-sig.88.lin.neg[1,8] 
intercept17<-ab.17[1,7] 
slope17<-ab.17[1,8] 
peaks <- read.csv("20231208_Averages.csv", header = TRUE) # to get elevation for graph
range(peaks$elevation)
png("fig.1.png")
curve(intercept +slope*x, from =261,
      to = 5581, col = "darkseagreen4",
      xlab = "Elevation",
      ylab = "NDVI",
      main = "Eastern Abercrombie Mountain Treeline Change",
      lwd = 3)
#compare with 2017 (add to same plot)
curve(intercept17 +slope17*x, from =261,
to = 5581, col = "peru", add = TRUE,
lwd = 3)
#add legend
legend("topright", legend = c("1984-1988", "2013-2017"), col = c("darkseagreen4", "peru"), lty = 1, lwd = 3)
dev.off()
##why are there ten elevation models for 1 peak, direction and year?