Collection <- read.csv("data/Larvae-collection.csv", header=T, stringsAsFactors = F)
Collection <- Collection[c(1,2,3,4,5,6,7,8,9,10,11,12,13,22)]
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
Collection <- merge(x=Collection, y=unique(survival[,c("TRT.REP", "TEMP", "FOOD")]), by.x = "Treatment", by.y = "TRT.REP", all.x=T, all.y=T)

Collection$TREAT <- as.factor(paste(Collection$TEMP, "-", Collection$FOOD))
colnames(Collection) <- c("Rep", "Date", "Bag", "Broodstock", "Group", "Vol.sampled", "Vol.total", "Live.A", "Dead.A", "Live.B", "Dead.B", "Live.C", "Dead.C", "Sample.number", "Live.Larvae", "Dead.Larvae", "Live.Larvae.norm" ,"TEMP", "FOOD", "TREAT")

Collection$Perc.live <- (Collection$Live.Larvae/(Collection$Dead.Larvae + Collection$Live.Larvae))*100
survival

saveRDS(Collection, file = "data/larvae-collection-data.rds") #save for later use 

nrow(subset(Collection, Live.Larvae>10000 & TREAT=="Cold - Low")) == (13+9+7+8)
nrow(subset(Collection, Live.Larvae>10000 & TREAT=="Cold - High")) == (8+4+7+9)
nrow(subset(Collection, Live.Larvae>10000 & TREAT=="Warm - Low")) == (13+12+9+12)
nrow(subset(Collection, Live.Larvae>10000 & TREAT=="Warm - High")) == (9+10+11+12)

# Calculate summary stats by treatment including reps 
treat_total.rep <- aggregate(cbind(Live.Larvae, Live.Larvae.norm) ~ Date + TREAT + TEMP + FOOD + Bag, data = Collection, sum, na.rm = TRUE) %>%
  group_by(TEMP, FOOD, TREAT, Bag) %>% 
  mutate(cum.total=cumsum(Live.Larvae),cum.percap = cumsum(Live.Larvae.norm),CalDay = format(Date,"%j")) %>% 
  group_by(TEMP, FOOD, TREAT, Bag) %>% dplyr::summarize(overall_Total = sum(Live.Larvae, na.rm = T), mean.Live.Larvae = mean(Live.Larvae, na.rm=T), mean.percap = mean(Live.Larvae.norm, na.rm=T), total.percap = sum(Live.Larvae.norm, na.rm=T), maxday = as.numeric(CalDay[which.max(Live.Larvae)]), max = max(Live.Larvae), max.percap = max(Live.Larvae.norm), first.big = as.numeric(CalDay[which(Live.Larvae > 10000)[1]]), release.days = as.numeric(length(CalDay[Live.Larvae > 10000])))
treat_total.rep <- merge(x=treat_total.rep, y=aggregate(Broodstock ~ Bag, data=Collection, median), by.x = "Bag", by.y = "Bag")[,-1]
treat_total.rep$trial <- c("13") 

# Estimate # and % broodstock that reproduced as females 
treat_total.rep$noFemales <- treat_total.rep$overall_Total/215000 #estimate of # females 
treat_total.rep$perc.spawn <- 100*(treat_total.rep$noFemales/treat_total.rep$Broodstock) #estimate of # females 
ggplot(data=treat_total.rep, aes(x=TEMP:FOOD, y=perc.spawn)) + geom_boxplot() + geom_point()
summary(treat_total.rep$perc.spawn)

# Treatment total (all reps combined), by date 
treat_total <- aggregate(cbind(Live.Larvae, Live.Larvae.norm) ~ Date + TREAT + TEMP + FOOD, data = Collection, sum, na.rm = TRUE)

# Stats 

# Daily larvae released - did it differ between treatments?
hist((na.omit(subset(Collection, Live.Larvae>0)$Live.Larvae))^(1/3))
qqnorm(na.omit(subset(Collection, Live.Larvae>0)$Live.Larvae)^(1/3))
qqline(na.omit(subset(Collection, Live.Larvae>0)$Live.Larvae)^(1/3))
shapiro.test(na.omit(subset(Collection, Live.Larvae>0)$Live.Larvae)^(1/3))  #hmmmm 
bartlett.test(Live.Larvae ~ TREAT, data=na.omit(subset(Collection, Live.Larvae>0))) #variance OK
leveneTest(na.omit(subset(Collection, Live.Larvae>0)$Live.Larvae)^(1/3), group=subset(Collection, Live.Larvae>0)$TREAT) #variance OK 
anova(lm(I(Live.Larvae^(1/3)) ~ TEMP*FOOD, data=subset(Collection, Live.Larvae>0))) #no effect

# Total release normalized by broodstock 
hist(treat_total.rep$total.percap)
shapiro.test(treat_total.rep$total.percap)
bartlett.test(total.percap ~ TREAT, data=treat_total.rep) #variances differ 
leveneTest(treat_total.rep$total.percap, group=treat_total.rep$TREAT) #borderline 
anova(lm(total.percap ~ TEMP*FOOD, data=treat_total.rep)) #no diff in total released per brood 
kruskal.test(total.percap ~ TREAT, data=treat_total.rep) #no diff in total released per brood (trying non-parametric)

?bartlett.test
?leveneTest
# Total release not normalized
shapiro.test(treat_total.rep$overall_Total)
leveneTest(treat_total.rep$overall_Total, group=treat_total.rep$TREAT)
summary(aov(overall_Total ~ TEMP*FOOD, data=treat_total.rep)) #no diff in total release
?aov()

# First big day 
shapiro.test(treat_total.rep$first.big) #can't easily convert to normal 
kruskal.test(first.big ~ TREAT, data = treat_total.rep) #Diff 
kruskal.test(first.big ~ TEMP, data = treat_total.rep) #no diff temp only 
kruskal.test(first.big ~ FOOD, data = treat_total.rep) #Diff 
plot(x=treat_total.rep$TREAT, treat_total.rep$first.big) # earlier release in Low food group plot(x=treat_total.rep$FOOD, treat_total.rep$first.big) # earlier release in Low food group 

# Max day  
shapiro.test(treat_total.rep$maxday) #can't easily convert to normal 
kruskal.test(maxday ~ TREAT, data = treat_total.rep) #no diff all grps
kruskal.test(maxday ~ TEMP, data = treat_total.rep) #no diff temp
kruskal.test(maxday ~ FOOD, data = treat_total.rep) #no diff all FOOD

# Estimated % spawn as females   
shapiro.test(treat_total.rep$perc.spawn) # normal 
hist(treat_total.rep$perc.spawn) # normal 
summary(aov(perc.spawn ~ TEMP*FOOD, data = treat_total.rep)) #No diff 
aggregate(perc.spawn ~ TEMP+FOOD, data=treat_total.rep, mean) #mean percent by treatments  
mean(treat_total.rep$perc.spawn)
sd(treat_total.rep$perc.spawn)

Collection.total <- aggregate(cbind(Live.Larvae, Live.Larvae.norm) ~ TREAT + TEMP + FOOD + Rep + Bag, data = Collection, sum, na.rm = TRUE)

plot(Live.Larvae.norm ~ TREAT, data=Collection.total, col=c("skyblue3", "seagreen3",  "indianred2","orange1"), main="Total larvae released per treatment\n13-week exposure\n4 reps, normalized by broodstock", xlab="Treatment", ylab="Cumulative larvae, normalized", cex.lab=1.5, cex.main=1.5, par(mar=c(5,5,4.1,2.1)))

plot(Live.Larvae.norm ~ Rep, data=Collection.total, col=c("skyblue3", "seagreen3",  "indianred2","orange1"), main="Total larvae released per treatment\n(4 reps, normalized by broodstock)", xlab="Treatment", ylab="Cumulative larvae, normalized", cex.lab=1.5, cex.main=1.5, par(mar=c(5,5,4.1,2.1)))

plot(Live.Larvae ~ Bag, data=Collection.total,  main="Total larvae released per treatment across reps", xlab="Treatment", ylab="Cumulative larvae, not normalized", cex.lab=1.5, cex.main=1.5, par(mar=c(5,5,4.1,2.1)))

png(filename = "results/total-larvae-released-points-12wk.png", width=400, height = 300)
ggplot(Collection.total, aes(x=TREAT, y=Live.Larvae.norm)) + geom_point(size=5, aes(color=TREAT)) + theme_bw() + labs(title="Total larvae released per treatment\n13-week exposure\n4 replicates per treatment", y=("Cumulative release per broodstock")) + scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food")) + theme(text = element_text(size = 12)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) + scale_y_continuous(lim=c(30000,225000))
dev.off()

# Summary stats on larval collection 
mean(Collection$Live.Larvae)

aggregate(Live.Larvae ~ TEMP + FOOD, data = Collection, mean, na.rm = TRUE)
aggregate(Live.Larvae ~ TEMP + FOOD, data = Collection, median, na.rm = TRUE)
aggregate(Live.Larvae ~ TEMP + FOOD, data = Collection, max, na.rm = TRUE)
aggregate(Live.Larvae ~ TEMP + FOOD, data = Collection, sum, na.rm = TRUE)
aggregate(Broodstock ~ TEMP + FOOD, data = Collection, median, na.rm = TRUE)
aggregate(Live.Larvae ~ TEMP + FOOD, data = Collection, sd, na.rm = TRUE)
aggregate(Live.Larvae ~ TEMP + FOOD + Bag, data = Collection, sum, na.rm = TRUE)
aggregate(Live.Larvae ~ TEMP + FOOD + Bag, data = Collection, sum, na.rm = TRUE)


#Calculate cumulative larvae released through time for plots 
Collection.cumul <- aggregate(cbind(Live.Larvae, Live.Larvae.norm) ~ Date + TREAT + TEMP + FOOD, data = Collection, sum, na.rm = TRUE) %>%
  group_by(TEMP, FOOD, TREAT) %>% 
  mutate(cum.total=cumsum(Live.Larvae),cum.percap = cumsum(Live.Larvae.norm),CalDay = format(Date,"%j")) %>% 
  arrange(Date) %>% dplyr::select(Date,CalDay,TREAT,TEMP,FOOD,Live.Larvae,Live.Larvae.norm, cum.total,cum.percap)

# Relevel food factors such that Low comes before High
Collection.cumul$FOOD <- factor(Collection.cumul$FOOD, levels = rev(levels(Collection.cumul$FOOD)))

TREATS <- levels(Collection.cumul$FOOD:Collection.cumul$TEMP)
TREATS.col <- c("#92c5de","#f4a582","#0571b0","#ca0020")

p <- list()
p[[1]] <- ggplot(data=Collection.cumul, aes(x=Date, y=cum.percap, group=FOOD:TEMP, color=FOOD:TEMP)) + 
    geom_line(size=.75) + 
  scale_color_manual(values=TREATS.col, name=element_blank()) + 
    theme_classic(base_size = 14) + 
    #ggtitle("Cumulative larvae released (normalized by no. broodstock)\n13-week exposure") + 
    scale_x_date(date_breaks = "1 week",date_labels ="%b-%d", 
                 limits=c(as.Date("2018-03-30"), as.Date("2018-04-27"))) + 
    theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), 
          plot.title = element_text(size = 14, hjust = 0), title = element_blank(), 
          axis.line = element_line(size = .5, colour = "gray50")) + 
  scale_y_continuous(limits = c(0, 600000), position = "right") 

for (i in 1:length(TREATS)) {
  p[[i+1]] <- ggplot(data=subset(Collection.cumul, FOOD:TEMP==TREATS[i]), aes(x=Date, y=Live.Larvae, fill=FOOD:TEMP)) + geom_bar(fill=TREATS.col[i], stat="identity",width=.5, position = position_dodge(width=2)) + 
    theme_classic(base_size = 14) +  
    theme(plot.title = element_text(size = 14, hjust = 0), axis.title.y = element_blank(), 
          axis.title.x = element_blank(),  axis.title = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none",
          axis.line = element_line(size = .5, colour = "gray50")) + 
    scale_x_date(date_breaks = "1 week",date_labels ="%b-%d", 
                 limits=c(as.Date("2018-03-30"), as.Date("2018-04-27"))) + 
    scale_y_continuous(labels = scales::scientific, 
                       limits = c(0,5780000), breaks=c(2500000, 5000000), position = "right") +
    geom_smooth(method="loess", color="gray60", size=0.6, linetype="dashed", span = 0.25, method.args = list(degree=1), se = FALSE)
}
pdf(file = "results/13wk-larvae.pdf", width = 5, height = 8)
plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], align = "v", nrow = 5, rel_heights = c(2/6, 1/6, 1/6, 1/6, 1/6))
dev.off()











### ========================= old crap 

# Number of substantial larval release days (>10k larvae)
shapiro.test(treat_total.rep$release.days)  
leveneTest(treat_total.rep$release.days, group=treat_total.rep$TREAT)
summary(aov(release.days ~ TEMP*FOOD, data=treat_total.rep)) # Temperature effect 
TukeyHSD(aov(release.days ~ TEMP*FOOD, data=treat_total.rep)) # Temperature effect
plot(x=treat_total.rep$TEMP, treat_total.rep$release.days) # more release days in Warm treatment 
plot(x=treat_total.rep$TREAT, y=treat_total.rep$release.days)

# Larvae released, normalized by # broodstock (all groups)
p2 <- ggplot(data=Collection.cumul, aes(x=Date, y=Live.Larvae.norm, fill=TREAT)) + 
  geom_bar(stat="identity",width=.5) + ylab("Daily larvae released") +
  ggtitle("Larvae release (normalized by broodstock)\n13-week exposure") + theme_bw(base_size = 16) +   
  theme(plot.title = element_text(size = 16, hjust = 0)) + 
  scale_x_date(date_breaks = "1 week",date_labels ="%b-%d") + 
  theme(legend.position = "bottom") + scale_fill_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food"))
Larvae.release.norm <- p2 + geom_line(data=Collection.cumul, aes(x=Date, y=cum.percap/5, group=TREAT, color=TREAT),size=.75) +
  scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food"))  + scale_y_continuous(limits=c(0,120000), labels = scales::comma, sec.axis = sec_axis(~.*5,name="Total larvae released, per broodstock"))
png(filename = "results/larval-release-chart-12wk.png", width = 700, height = 450)
print(Larvae.release.norm)
dev.off()

# Larvae released (Cold - High)
p3 <- ggplot(data=subset(Collection.cumul, TREAT=="Cold - High"), aes(x=Date, y=Live.Larvae, fill=TREAT)) + geom_bar(fill="#0571b0", stat="identity",width=.5, position = position_dodge(width=2)) + ylab("Number of Larvae Released") +
  ggtitle("Cold / High Food, 13-week exposure \nLarvae Release") + theme_bw(base_size = 14) +   
  theme(plot.title = element_text(size = 14, hjust = 0)) + 
  scale_x_date(date_breaks = "1 week",date_labels ="%b-%d", limits=c(as.Date("2018-03-30"), as.Date("2018-04-27"))) + theme(legend.position = "none")

png(filename = "results/larval-release-chart-Cold-High-12wk.png", width = 500, height = 400)
print(Larvae.release.Cold.High <- p3+ geom_line(data=subset(Collection.cumul, TREAT=="Cold - High"), aes(x=Date, y=cum.total/5, group=TREAT, color=TREAT),size=.75) + scale_color_manual(values="#0571b0") + scale_y_continuous(limits = c(0, 6e6), sec.axis = sec_axis(~.*6.5,name="Cumulative Larvae Released")))
dev.off()

# Larvae released (Cold - Low)
p4 <- ggplot(data=subset(Collection.cumul, TREAT=="Cold - Low"), aes(x=Date, y=Live.Larvae, fill=TREAT)) + geom_bar(fill="#92c5de", stat="identity",width=.5, position = position_dodge(width=2)) + ylab("Number of Larvae Released") +
  ggtitle("Cold / Low Food, 13-week exposure \nLarvae Release") + theme_bw(base_size = 14) +   
  theme(plot.title = element_text(size = 14, hjust = 0)) + 
  scale_x_date(date_breaks = "1 week",date_labels ="%b-%d", limits=c(as.Date("2018-03-30"), as.Date("2018-04-27"))) + 
  theme(legend.position = "none") + scale_fill_manual(values=c("seagreen3"))

png(filename = "results/larval-release-chart-Cold-Low-12wk.png", width = 500, height = 400)
print(Larvae.release.Cold.Low <- p4+ geom_line(data=subset(Collection.cumul, TREAT=="Cold - Low"), aes(x=Date, y=cum.total/5, group=TREAT, color=TREAT),size=.75) + scale_color_manual(values="#92c5de") + scale_y_continuous(limits = c(0, 6e6), sec.axis = sec_axis(~.*6.5,name="Cumulative Larvae Released")))
dev.off()

# Larvae released (Warm - High)
p5 <- ggplot(data=subset(Collection.cumul, TREAT=="Warm - High"), aes(x=Date, y=Live.Larvae, fill=TREAT)) + geom_bar(fill="#ca0020", stat="identity",width=.5, position = position_dodge(width=2)) + ylab("Number of Larvae Released") +
  ggtitle("Warm / High Food, 13-week exposure \nLarvae Release") + theme_bw(base_size = 14) +   
  theme(plot.title = element_text(size = 14, hjust = 0)) + 
  scale_x_date(date_breaks = "1 week",date_labels ="%b-%d", limits=c(as.Date("2018-03-30"), as.Date("2018-04-27"))) + theme(legend.position = "none")

png(filename = "results/larval-release-chart-Warm-High-12wk.png", width = 500, height = 400)
print(Larvae.release.Warm.High <- p5+ geom_line(data=subset(Collection.cumul, TREAT=="Warm - High"), aes(x=Date, y=cum.total/5, group=TREAT, color=TREAT),size=.75) + scale_color_manual(values="#ca0020") + scale_y_continuous(limits = c(0, 6e6), sec.axis = sec_axis(~.*6.5,name="Cumulative Larvae Released")))
dev.off()

# Larvae released (Warm - Low)
p6 <- ggplot(data=subset(Collection.cumul, TREAT=="Warm - Low"), aes(x=Date, y=Live.Larvae, fill=TREAT)) +  geom_bar(fill="#f4a582", stat="identity",width=.5, position = position_dodge(width=2)) + ylab("Number of Larvae Released") +
  ggtitle("Warm / Low Food, 13-week exposure \nLarvae Release") + theme_bw(base_size = 14) +   
  theme(plot.title = element_text(size = 14, hjust = 0)) + 
  scale_x_date(date_breaks = "1 week",date_labels ="%b-%d", limits=c(as.Date("2018-03-30"), as.Date("2018-04-27"))) + theme(legend.position = "none")

png(filename = "results/larval-release-chart-Warm-Low-12wk.png", width = 500, height = 400)
print(Larvae.release.Warm.Low <- p6+ geom_line(data=subset(Collection.cumul, TREAT=="Warm - Low"), aes(x=Date, y=cum.total/5, group=TREAT, color=TREAT),size=.75) + scale_color_manual(values="#f4a582") + scale_y_continuous(limits = c(0, 6e6), sec.axis = sec_axis(~.*6.5,name="Cumulative Larvae Released")))
dev.off()