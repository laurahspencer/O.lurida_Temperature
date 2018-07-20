library(car)

Collection <- read.csv("data/Larvae-collection.csv", header=T, stringsAsFactors = F)
Collection <- Collection[c(1,2,3,4,6,7,8,9,10,11,12,13)]
Collection$Date.collected <- as.Date(Collection$Date.collected, format = "%m/%d/%y")
Collection$Treatment <- as.factor(Collection$Treatment)
Collection$Bag.. <- as.factor(Collection$Bag..)
Collection$Broodstock <- as.numeric(Collection$Broodstock)
Collection$Vol..for.count..mL. <- as.numeric(Collection$Vol..for.count..mL.)
Collection$Count.A.LIVE <- as.numeric(Collection$Count.A.LIVE)
Collection$Count.A.DEAD <- as.numeric(Collection$Count.A.DEAD)
Collection$Count.B.LIVE <- as.numeric(Collection$Count.B.LIVE)
Collection$Count.B.DEAD <- as.numeric(Collection$Count.B.DEAD)
Collection$Count.C.LIVE <- as.numeric(Collection$Count.C.LIVE)
Collection$Count.C.DEAD <- as.numeric(Collection$Count.C.DEAD)
Collection$Live.Larvae <- (rowMeans(subset(Collection, select = c("Count.A.LIVE","Count.C.LIVE","Count.C.LIVE")), na.rm = TRUE))/Collection$Vol..for.count..mL.*Collection$Total.Vol..mL.
Collection$Dead.Larvae <- (rowMeans(subset(Collection, select = c("Count.A.DEAD","Count.C.DEAD","Count.C.DEAD")), na.rm = TRUE))/Collection$Vol..for.count..mL.*Collection$Total.Vol..mL.
Collection$Live.Larvae.norm <- Collection$Live.Larvae/Collection$Broodstock
Collection <- merge(x=Collection, y=unique(survival[,c("TRT.REP", "TEMP", "FOOD")]), by.x = "Treatment", by.y = "TRT.REP", all.x=T, all.y=FALSE)
Collection$TREAT <- as.factor(paste(Collection$TEMP, "-", Collection$FOOD))
colnames(Collection) <- c("Rep", "Date", "Bag", "Broodstock", "Vol.sampled", "Vol.total", "Live.A", "Dead.A", "Live.B", "Dead.B", "Live.C", "Dead.C", "Live.Larvae", "Dead.Larvae", "Live.Larvae.norm" ,"TEMP", "FOOD", "TREAT")

summary(Collection$Date)
treat_total <- aggregate(cbind(Live.Larvae, Live.Larvae.norm) ~ Date + TREAT + TEMP + FOOD, data = Collection, sum, na.rm = TRUE)

#Calculate cumulative larvae released through time
Collection.cumul <- group_by(treat_total, TREAT) %>% 
  mutate(cum.total=cumsum(Live.Larvae),cum.percap = cumsum(Live.Larvae.norm),CalDay = format(Date,"%j")) %>% 
  arrange(Date) %>% dplyr::select(Date,CalDay,TREAT,TEMP,FOOD,Live.Larvae,Live.Larvae.norm,cum.total,cum.percap)

library(ggplot2)
# Larvae released (all groups )
p1 <- ggplot(data=Collection.cumul, aes(x=Date, y=Live.Larvae, fill=TREAT)) + 
  geom_bar(stat="identity",width=.5, position = position_dodge(width=2)) + ylab("Number of Larvae Released") +
  ggtitle("O. lurida 2018 Pre-Conditioning Experiment\nLarvae Release by Treatment (not normalized by broodstock)") + theme_bw(base_size = 16) +   
  theme(plot.title = element_text(face = 'bold',size = 20, hjust = 0)) + 
  scale_x_date(date_breaks = "1 week",date_labels ="%b-%d") + 
  theme(legend.position = c(0.15, 0.85)) + scale_fill_manual(values=c("skyblue3", "seagreen3", "indianred2",  "orange1"))
Larvae.release <- p1+ geom_line(data=Collection.cumul, aes(x=Date, y=cum.total/5, group=TREAT, color=TREAT),size=.75) +
  scale_color_manual(values=c("royalblue4", "seagreen4", "indianred4",  "orange3")) + scale_y_continuous(sec.axis = sec_axis(~.*5,name="Cumulative Larvae Released"))

# Larvae released, normalized by # broodstock (all groups)
p2 <- ggplot(data=Collection.cumul, aes(x=Date, y=Live.Larvae.norm, fill=TREAT)) + 
  geom_bar(stat="identity",width=.5, position = position_dodge(width=2)) + ylab("Larvae Released, per broodstock") +
  ggtitle("O. lurida 2018 Pre-Conditioning Experiment\nLarvae Release by Treatment (normalized by broodstock)") + theme_bw(base_size = 16) +   
  theme(plot.title = element_text(face = 'bold',size = 20, hjust = 0)) + 
  scale_x_date(date_breaks = "1 week",date_labels ="%b-%d") + 
  theme(legend.position = c(0.15, 0.85)) + scale_fill_manual(values=c("skyblue3", "seagreen3", "indianred2",  "orange1"))
Larvae.release.norm <- p2+ geom_line(data=Collection.cumul, aes(x=Date, y=cum.percap/5, group=TREAT, color=TREAT),size=.75) +
  scale_color_manual(values=c("royalblue4", "seagreen4", "indianred4",  "orange3")) + scale_y_continuous(sec.axis = sec_axis(~.*5,name="Cumulative Larvae Released, per broodstock"))
levels(Collection.cumul$TREAT)

# Larvae released (COLD - HIGH)
p3 <- ggplot(data=subset(Collection.cumul, TREAT=="COLD - HIGH"), aes(x=Date, y=Live.Larvae, fill=TREAT)) + 
  geom_bar(stat="identity",width=.5, position = position_dodge(width=2)) + ylab("Number of Larvae Released") +
  ggtitle("7C/High Food Treatment, \nLarvae Release") + theme_bw(base_size = 16) +   
  theme(plot.title = element_text(face = 'bold',size = 20, hjust = 0)) + 
  scale_x_date(date_breaks = "1 week",date_labels ="%b-%d", limits=c(as.Date("2018-03-30"), as.Date("2018-04-27"))) + 
  theme(legend.position = "none") + scale_fill_manual(values=c("skyblue3", "seagreen3", "indianred2",  "orange1"))
Larvae.release.COLD.HIGH <- p3+ geom_line(data=subset(Collection.cumul, TREAT=="COLD - HIGH"), aes(x=Date, y=cum.total/5, group=TREAT, color=TREAT),size=.75) + scale_color_manual(values="royalblue4") + scale_y_continuous(sec.axis = sec_axis(~.*5,name="Cumulative Larvae Released"))

# Larvae released (COLD - LOW)
p4 <- ggplot(data=subset(Collection.cumul, TREAT=="COLD - LOW"), aes(x=Date, y=Live.Larvae, fill=TREAT)) + 
  geom_bar(stat="identity",width=.5, position = position_dodge(width=2)) + ylab("Number of Larvae Released") +
  ggtitle("7C/Low Food Treatment, \nLarvae Release") + theme_bw(base_size = 16) +   
  theme(plot.title = element_text(face = 'bold',size = 20, hjust = 0)) + 
  scale_x_date(date_breaks = "1 week",date_labels ="%b-%d", limits=c(as.Date("2018-03-30"), as.Date("2018-04-27"))) + 
  theme(legend.position = "none") + scale_fill_manual(values=c("seagreen3"))
Larvae.release.COLD.LOW <- p4+ geom_line(data=subset(Collection.cumul, TREAT=="COLD - LOW"), aes(x=Date, y=cum.total/5, group=TREAT, color=TREAT),size=.75) + scale_color_manual(values="seagreen4") + scale_y_continuous(sec.axis = sec_axis(~.*5,name="Cumulative Larvae Released"))

summary(Collection.cumul$Date)

# Larvae released (WARM - HIGH)
p5 <- ggplot(data=subset(Collection.cumul, TREAT=="WARM - HIGH"), aes(x=Date, y=Live.Larvae, fill=TREAT)) + 
  geom_bar(stat="identity",width=.5, position = position_dodge(width=2)) + ylab("Number of Larvae Released") +
  ggtitle("10C/High Food Treatment, \nLarvae Release") + theme_bw(base_size = 16) +   
  theme(plot.title = element_text(face = 'bold',size = 20, hjust = 0)) + 
  scale_x_date(date_breaks = "1 week",date_labels ="%b-%d", limits=c(as.Date("2018-03-30"), as.Date("2018-04-27"))) + 
  theme(legend.position = "none") + scale_fill_manual(values=c("indianred2"))
Larvae.release.WARM.HIGH <- p5+ geom_line(data=subset(Collection.cumul, TREAT=="WARM - HIGH"), aes(x=Date, y=cum.total/5, group=TREAT, color=TREAT),size=.75) + scale_color_manual(values="indianred2") + scale_y_continuous(sec.axis = sec_axis(~.*5,name="Cumulative Larvae Released"))

# Larvae released (WARM - LOW)
p6 <- ggplot(data=subset(Collection.cumul, TREAT=="WARM - LOW"), aes(x=Date, y=Live.Larvae, fill=TREAT)) +  geom_bar(stat="identity",width=.5, position = position_dodge(width=2)) + ylab("Number of Larvae Released") +
  ggtitle("10C/Low Food Treatment, \nLarvae Release") + theme_bw(base_size = 16) +   
  theme(plot.title = element_text(face = 'bold',size = 20, hjust = 0)) + 
  scale_x_date(date_breaks = "1 week",date_labels ="%b-%d", limits=c(as.Date("2018-03-30"), as.Date("2018-04-27"))) + 
  theme(legend.position = "none") + scale_fill_manual(values=c("orange1"))
Larvae.release.WARM.LOW <- p6+ geom_line(data=subset(Collection.cumul, TREAT=="WARM - LOW"), aes(x=Date, y=cum.total/5, group=TREAT, color=TREAT),size=.75) + scale_color_manual(values="orange1") + scale_y_continuous(sec.axis = sec_axis(~.*5,name="Cumulative Larvae Released"))

summary(na.omit(Collection$Live.Larvae.norm))
hist((na.omit(Collection$Live.Larvae.norm))^(0.25))
qqnorm(na.omit(Collection$Live.Larvae.norm)^(0.25))
qqline(na.omit(Collection$Live.Larvae.norm)^(0.25))
shapiro.test(na.omit(Collection$Live.Larvae.norm)^(0.25)) #normal distributed after transformation
leveneTest(na.omit(Collection$Live.Larvae.norm)^(0.25), group=Collection$TREAT) #variance OK 

larvae.temp.aov <- aov(Live.Larvae.norm^(0.25) ~ TEMP, data= Collection)
summary(larvae.temp.aov) #no diff, but higher F-value 

larvae.food.aov <- aov(Live.Larvae.norm^(0.25) ~ FOOD, data= Collection)
summary(larvae.food.aov) #no diff. 

# Nested anova on normalized larval release data 
larvae.lm <- lm(Live.Larvae.norm^(0.25) ~ TREAT/Rep/Bag, data= Collection)
anova(larvae.lm)
larvae.lm <- lm(Live.Larvae.norm^(0.25) ~ TREAT/Rep, data= Collection)
anova(larvae.lm)
larvae.lm <- lm(Live.Larvae.norm^(0.25) ~ TREAT, data= Collection)
anova(larvae.lm)
summary(larvae.lm)
larvae.aov <- aov(Live.Larvae.norm^(0.25) ~ TREAT, data= Collection)
summary(larvae.aov)

Collection.total <- aggregate(cbind(Live.Larvae, Live.Larvae.norm) ~  TREAT + TEMP + FOOD + Rep + Bag, data = Collection, sum, na.rm = TRUE)
library(reshape2)
chisq.test(x=Collection.total$TEMP, y=Collection.total$Live.Larvae.norm, simulate.p.value = T)
chisq.test(x=Collection.total$FOOD, y=Collection.total$Live.Larvae.norm, simulate.p.value = T)

glm.cum.larvae <- glm(Live.Larvae.norm ~ TEMP*FOOD, data=Collection.total)
summary(glm.cum.larvae)
anodev.fn(glm.cum.larvae)

glm.cum.larvae2 <- update(glm.cum.larvae, ~.- TEMP:LOW) 
anova(glm.cum.larvae,glm.cum.larvae2,test="Chi") #interaction not significant 

glm.cum.larvae3 <- glm(Live.Larvae.norm ~ TEMP+FOOD, data=Collection.total)
summary(glm.cum.larvae3)
anodev.fn(glm.cum.larvae)


plot(Live.Larvae.norm ~ TREAT, data=Collection.total, col=c("skyblue3", "seagreen3",  "indianred2","orange1"), main="Total larvae released per treatment\n(4 reps, normalized by broodstock)", xlab="Treatment", ylab="Cumulative larvae, normalized", cex.lab=1.5, cex.main=1.5, par(mar=c(5,5,4.1,2.1)))

plot(Live.Larvae.norm ~ Rep, data=Collection.total, col=c("skyblue3", "seagreen3",  "indianred2","orange1"), main="Total larvae released per treatment\n(4 reps, normalized by broodstock)", xlab="Treatment", ylab="Cumulative larvae, normalized", cex.lab=1.5, cex.main=1.5, par(mar=c(5,5,4.1,2.1)))

plot(Live.Larvae ~ Bag, data=Collection.total,  main="Total larvae released per treatment across reps", xlab="Treatment", ylab="Cumulative larvae, not normalized", cex.lab=1.5, cex.main=1.5, par(mar=c(5,5,4.1,2.1)))

ggplot(Collection.total, aes(x=TREAT, y=Live.Larvae.norm)) + geom_point(size=6, aes(color=TREAT)) + labs(title="Total larvae released per treatment (points = reps)", y=("Cumulative release per broodstock")) + scale_color_manual(values=c("skyblue3", "seagreen3", "indianred2",  "orange1"))+ scale_shape_manual(values=c(16, 8)) + theme(text = element_text(size = 20)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

# Summary stats on larval collection 
mean(Collection$Live.Larvae)

aggregate(Live.Larvae ~ TEMP + FOOD, data = Collection, mean, na.rm = TRUE)
aggregate(Live.Larvae ~ TEMP + FOOD, data = Collection, median, na.rm = TRUE)
aggregate(Live.Larvae ~ TEMP + FOOD, data = Collection, max, na.rm = TRUE)
aggregate(Live.Larvae ~ TEMP + FOOD, data = Collection, sum, na.rm = TRUE)

aggregate(Broodstock ~ TEMP + FOOD, data = Collection, median, na.rm = TRUE)
