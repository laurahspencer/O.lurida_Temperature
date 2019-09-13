# 13 week trial 

temp.brood <- read.csv(file="data/Avtech-Temps_2017-11-20_2018-04-27.csv")
colnames(temp.brood) <- c("Timestamp", "A1", "A2", "B1", "B2", "C1", "C2", "D1", "D2")
temp.brood$Timestamp <- as.POSIXct(temp.brood$Timestamp, format="%m/%d/%y %H:%M")
temp.brood <- melt(data = temp.brood, id.vars = "Timestamp")
colnames(temp.brood) <- c("Timestamp", "trt.rep", "temperature")
ggplot(data = temp.brood, aes(x=Timestamp, y=temperature, col=trt.rep)) + geom_line()
