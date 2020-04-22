# 13 week trial 

temp.brood <- read.csv(file=here::here("data","Experimental-temperatures","Avtech-Temps_2017-11-20_2018-04-27.csv"))
colnames(temp.brood) <- c("Date", "A1", "A2", "B1", "B2", "C1", "C2", "D1", "D2")
temp.brood$Date <- as.POSIXct(temp.brood$Date, format="%m/%d/%y %H:%M")
temp.brood <- melt(data = temp.brood, id.vars = "Date")
colnames(temp.brood) <- c("Date", "trt.rep", "Temperature")
pdf(file = "results/treatment-temperatures.pdf", width = 8, height = 4)
ggplotly(ggplot(data = temp.brood, aes(x=Date, y=Temperature, col=trt.rep)) + geom_line() + theme_minimal() + theme(legend.position = "none") + 
  scale_color_manual(values=c('#92c5de', '#92c5de', '#f4a582', '#f4a582', '#0571b0', '#0571b0', '#ca0020', '#ca0020')))
dev.off()

aggregate(Temperature ~ trt.rep, data=subset(temp.brood, Date>"2017-12-15 09:30:00" & Date<"2018-02-27 10:00:00"), mean) 
aggregate(Temperature ~ trt.rep, data=subset(temp.brood, Date>"2017-12-15 09:30:00" & Date<"2018-02-27 10:00:00"), sd) 
 
# Mud Bay 

# shallow 
temp.mud.shallow1 <- read.csv(file = "data/Mud-Bay-Temperatures/2018-02-12-MudBay-Shallow-temp.csv")[2:3]
colnames(temp.mud.shallow1) <- c("Date", "Temperature")
temp.mud.shallow2 <- read.csv(file="data/Mud-Bay-Temperatures/2018-02-28-MudBay-Shallow-temp-redeployment.csv")[2:3]
colnames(temp.mud.shallow2) <- c("Date", "Temperature")
temp.mud.shallow <- rbind(temp.mud.shallow1, temp.mud.shallow2)
temp.mud.shallow$Date <- as.POSIXct(temp.mud.shallow$Date, format="%m/%d/%y %H:%M")
pdf(file = "results/mudbay-shallow.pdf", width = 7, height = 4)
subset(temp.mud.shallow, Date>"2017-07-11 15:45:00" & Date<"2018-02-09 01:30:00") %>%
ggplot(aes(x=Date, y=Temperature)) + geom_line() + theme_minimal()+ ggtitle("Shallow Deployment") + ylim(0,15)
dev.off()

mean(subset(temp.mud.shallow, Date>"2017-07-11 15:45:00" & Date<"2018-02-09 01:30:00")[,"Temperature"]) #9.31C
sd(subset(temp.mud.shallow, Date>"2017-07-11 15:45:00" & Date<"2018-02-09 01:30:00")[,"Temperature"]) #2.58

# deep
temp.mud.deep <- read.csv(file = "data/Mud-Bay-Temperatures/2018-02-28_MudBay-Deep-temp.csv")[,2:3]
colnames(temp.mud.deep) <- c("Date", "Temperature")
temp.mud.deep$Date <- as.POSIXct(temp.mud.deep$Date, format="%m/%d/%y %H:%M")
pdf(file = "results/mudbay-deep.pdf", width = 7, height = 4)
subset(temp.mud.deep, Date>"2017-07-11 15:45:00" & Date<"2018-02-09 01:30:00") %>%
  ggplot(aes(x=Date, y=Temperature)) + geom_line() + theme_minimal() + ggtitle("Deep Deployment") + ylim(0,15)
dev.off()
mean(subset(temp.mud.deep, Date>"2017-07-11 15:45:00" & Date<"2018-02-09 01:30:00")[,"Temperature"]) #9.33C
sd(subset(temp.mud.deep, Date>"2017-07-11 15:45:00" & Date<"2018-02-09 01:30:00")[,"Temperature"]) #2.65C


# very deep
temp.mud.very.deep <- read.csv(file = "data/Mud-Bay-Temperatures/2018-06-3_MudBay-Very-Deep-temp.csv")[,2:3]
colnames(temp.mud.very.deep) <- c("Date", "Temperature")
temp.mud.very.deep$Date <- as.POSIXct(temp.mud.very.deep$Date, format="%m/%d/%y %H:%M")
pdf(file = "results/mudbay-very-deep.pdf", width = 7, height = 4)
subset(temp.mud.very.deep, Date>"2017-07-11 15:45:00" & Date<"2018-02-09 01:30:00") %>%
  ggplot(aes(x=Date, y=Temperature)) + geom_line() + theme_minimal()+ ggtitle("Very Deep, Mid-Channel Deployment") + ylim(0,15)
dev.off()
mean(subset(temp.mud.very.deep, Date>"2017-07-11 15:45:00" & Date<"2018-02-09 01:30:00")[,"Temperature"]) #9.13C
sd(subset(temp.mud.very.deep, Date>"2017-07-11 15:45:00" & Date<"2018-02-09 01:30:00")[,"Temperature"]) #2.63C



