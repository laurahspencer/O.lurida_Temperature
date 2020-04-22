
# Mortality after 6 week treatment 
brood.mortality.6wk <- read.csv("data/Broodstock-mortality-6wk-trt.csv", stringsAsFactors = F)
brood.mortality.6wk$Date <- as.Date(brood.mortality.6wk$Date, format = "%m/%d/%y")
brood.mortality.6wk$TRT <- factor(brood.mortality.6wk$TRT, levels=c("C", "A", "D", "B"))
str(brood.mortality.6wk)

brood.mortality.6wk.sum <- aggregate(Alive ~ Date + TRT, data = brood.mortality.6wk, min, na.rm = TRUE)

ggplot(data=brood.mortality.6wk.sum, aes(x=Date, y=100*(Alive/165), group=TRT, col=TRT)) + theme_bw(base_size = 13) + ggtitle("Broodstock survival over time\n6 week treatment") + geom_line()+ geom_point() + xlab("Date") + ylab("% Alive") + scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food"))

# Larvae collected after 6 week treatment 
Collection.6wks <- read.csv("data/Larvae-collection-6wk-trt.csv", header=T, stringsAsFactors = F)
Collection.6wks$Date <- as.Date(Collection.6wks$Date, format = "%m/%d/%y")
Collection.6wks$Bucket <- as.factor(Collection.6wks$Bucket)
Collection.6wks$Temperature <- as.factor(Collection.6wks$Temperature)
Collection.6wks$Food <- as.factor(Collection.6wks$Food)
Collection.6wks$position <- as.factor(Collection.6wks$position)
Collection.6wks$TREAT <- as.factor(paste(Collection.6wks$Temperature, "-", Collection.6wks$Food))
Collection.6wks$Larvae.norm <- Collection.6wks$Larvae/Collection.6wks$Alive.buck
Collection.6wks <- subset(Collection.6wks, Date<"2018-03-22") #IMPORTANT - only retain first 4 weeks of larval coll. data

plot(x=Collection.6wks$position, y=Collection.6wks$Larvae)
hist(Collection.6wks$Larvae)
hist(subset(Collection.6wks, Larvae!=0)$Larvae)

# Calculate summary stats by treatment including reps 
treat_total.6wks.rep <- aggregate(cbind(Larvae, Larvae.norm) ~ Date + TREAT + Temperature + Food + Bucket, data = Collection.6wks, sum, na.rm = TRUE) %>%
  group_by(Temperature, Food, TREAT, Bucket) %>% 
  mutate(cum.total=cumsum(Larvae),cum.percap = cumsum(Larvae.norm),CalDay = format(Date,"%j")) %>% 
  group_by(Temperature, Food, TREAT, Bucket) %>% dplyr::summarize(overall_Total = sum(Larvae, na.rm = T), mean.larvae = mean(Larvae, na.rm=T), mean.percap = mean(Larvae.norm, na.rm=T), total.percap = sum(Larvae.norm, na.rm=T), maxday = as.numeric(CalDay[which.max(Larvae)]), max = max(Larvae), max.percap = max(Larvae.norm), first.big = as.numeric(CalDay[which(Larvae > 10000)[1]]), release.days = as.numeric(length(CalDay[Larvae > 10000])))
treat_total.6wks.rep <- merge(x=treat_total.6wks.rep, y=aggregate(Alive.buck ~ Bucket, data=Collection.6wks, mean), by.x = "Bucket", by.y = "Bucket")[,-1]
treat_total.6wks.rep$trial <- c("6")

# Estimate # and % broodstock that reproduced as females 
treat_total.6wks.rep$noFemales <- round(treat_total.6wks.rep$overall_Total/215000, digits = 0) #estimate of # females 
treat_total.6wks.rep$perc.spawn <- 100*(treat_total.6wks.rep$noFemales/treat_total.6wks.rep$Alive.buck) #estimate of # females 
ggplot(data=treat_total.6wks.rep, aes(x=Temperature:Food, y=perc.spawn)) + geom_boxplot() + geom_point()
summary(treat_total.6wks.rep$perc.spawn)

#Calculate cumulative larvae released through time without reps 
Collection.6wks.cumul <- aggregate(cbind(Larvae, Larvae.norm) ~ Date + TREAT + Temperature + Food, data = Collection.6wks, sum, na.rm = TRUE) %>%
  group_by(Temperature, Food, TREAT) %>% 
  mutate(cum.total=cumsum(Larvae),cum.percap = cumsum(Larvae.norm),CalDay = format(Date,"%j")) %>% 
  arrange(Date) %>% dplyr::select(Date,CalDay,TREAT,Temperature,Food,Larvae,Larvae.norm, cum.total,cum.percap)

# Stats 

# Ave. no. larvae released per day - did it differ between treatments?
hist((na.omit(subset(Collection.6wks.cumul, Larvae>0)$Larvae))^(1/3))
qqnorm(na.omit(subset(Collection.6wks.cumul, Larvae>0)$Larvae)^(1/3))
qqline(na.omit(subset(Collection.6wks.cumul, Larvae>0)$Larvae)^(1/3))
shapiro.test(na.omit(subset(Collection.6wks.cumul, Larvae>0)$Larvae)^(1/3))  #hmmmm 
leveneTest(na.omit(subset(Collection.6wks.cumul, Larvae>0)$Larvae)^(1/3), group=subset(Collection.6wks.cumul, Larvae>0)$TREAT) #variance OK 
anova(lm(I(Larvae^(1/3)) ~ Temperature*Food, data=subset(Collection.6wks.cumul, Larvae>0))) #no effect

# Total release normalized by broodstock 
shapiro.test(treat_total.6wks.rep$total.percap)
leveneTest(treat_total.6wks.rep$total.percap, group=treat_total.6wks.rep$TREAT)
anova(lm(total.percap ~ Temperature*Food, data=treat_total.6wks.rep)) #no diff in total released per brood 

# Total release not normalized
shapiro.test(treat_total.6wks.rep$overall_Total)
leveneTest(treat_total.6wks.rep$overall_Total, group=treat_total.6wks.rep$TREAT)
anova(lm(overall_Total ~ Temperature*Food, data=treat_total.6wks.rep)) #no diff in total release

# First big day 
shapiro.test(treat_total.6wks.rep$first.big) #can't easily convert to normal 
kruskal.test(first.big ~ TREAT, data = treat_total.6wks.rep) #no diff all grps
kruskal.test(first.big ~ Temperature, data = treat_total.6wks.rep) #no diff temp
kruskal.test(first.big ~ Food, data = treat_total.6wks.rep) #no diff all food

# Max day  
shapiro.test(treat_total.6wks.rep$maxday) #can't easily convert to normal 
kruskal.test(maxday ~ TREAT, data = treat_total.6wks.rep) #no diff all grps
kruskal.test(maxday ~ Temperature, data = treat_total.6wks.rep) #no diff temp
kruskal.test(maxday ~ Food, data = treat_total.6wks.rep) #no diff all food

# Estimated % spawn as females   
shapiro.test(treat_total.6wks.rep$perc.spawn) # normal 
hist(treat_total.6wks.rep$perc.spawn) # normal 
summary(aov(perc.spawn ~ Temperature*Food, data = treat_total.6wks.rep)) #No diff 
aggregate(perc.spawn ~ Temperature+Food, data=treat_total.6wks.rep, mean) #mean percent by treatments  
mean(treat_total.6wks.rep$perc.spawn)
sd(treat_total.6wks.rep$perc.spawn)

# Total larvae released (normalized by # broodstock)
png(filename = "results/total-larvae-released-points-6wks.png", width=400, height = 300)
ggplot(treat_total.6wks.rep, aes(x=TREAT, y=total.percap)) + geom_point(size=5, aes(color=TREAT)) + theme_bw() + labs(title="Total larvae released per treatment\n6-week exposure\n3 replicates per treatment", y=("Total released (normalized/broodstock)")) + scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food")) + theme(text = element_text(size = 12)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) + scale_y_continuous(lim=c(30000,225000))
dev.off()


# Summary info on larval Collection.6wks 
mean(Collection.6wks$Larvae)
aggregate(Larvae ~ Temperature + Food, data = Collection.6wks, mean, na.rm = TRUE)
aggregate(Larvae ~ Temperature + Food, data = Collection.6wks, median, na.rm = TRUE)
aggregate(Larvae ~ Temperature + Food, data = Collection.6wks, max, na.rm = TRUE)
aggregate(Larvae ~ Temperature + Food, data = Collection.6wks, sum, na.rm = TRUE)
aggregate(Larvae ~ Temperature + Food, data = Collection.6wks, sd, na.rm = TRUE)


### ==========
# New plots for paper 8/26/2019

p <- list()
p[[1]] <- ggplot(data=Collection.6wks.cumul, aes(x=Date, y=cum.percap, group=TREAT, color=TREAT)) + 
  geom_line(size=.75) + scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food")) + 
  theme_classic(base_size = 14) + 
  #ggtitle("Cumulative larvae released (normalized by no. broodstock)\n13-week exposure") + 
  scale_x_date(date_breaks = "1 week",date_labels ="%b-%d") + 
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), 
        plot.title = element_text(size = 14, hjust = 0), title = element_blank(), 
        axis.line = element_line(size = .5, colour = "gray50"), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()) + 
  scale_y_continuous(labels = scales::comma, position = "right") 
max(Collection.6wks.cumul$Larvae)
for (i in 1:length(TREATS)) {
  p[[i+1]] <- ggplot(data=subset(Collection.6wks.cumul, TREAT==TREATS[i]), aes(x=Date, y=Larvae, fill=TREAT)) + geom_bar(fill=TREATS.col[i], stat="identity",width=.5, position = position_dodge(width=2)) + 
    theme_classic(base_size = 14) +  
    theme(plot.title = element_text(size = 14, hjust = 0), axis.title.y = element_blank(), 
          axis.title.x = element_blank(),  axis.title = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none",
          axis.line = element_line(size = .5, colour = "gray50"), axis.text.y = element_blank(), 
          axis.ticks.y = element_blank()) + 
    scale_x_date(date_breaks = "1 week",date_labels ="%b-%d", 
                 limits=c(as.Date("2018-02-19"), as.Date("2018-03-21"))) + 
    scale_y_continuous(labels = scales::scientific, limits = c(0,5780000), 
                       breaks=c(2500000, 5000000),position = "right") +
    geom_smooth(color="gray60", size=0.6, linetype="dashed", 
                span = 0.25, method.args = list(degree=1), se = FALSE)
}
pdf(file = "results/6wk-larvae.pdf", width = 4.5, height = 8)
plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], align = "v", nrow = 5, rel_heights = c(2/6, 1/6, 1/6, 1/6, 1/6))
dev.off()














#===
# Old crap, not to be used 

# Number of substantial larval release days (>10k larvae)
shapiro.test(treat_total.6wks.rep$release.days) #can't easily convert to normal 
kruskal.test(release.days ~ TREAT, data = treat_total.6wks.rep) #no diff all grps
kruskal.test(release.days ~ Temperature, data = treat_total.6wks.rep) #no diff temp
kruskal.test(release.days ~ Food, data = treat_total.6wks.rep) #no diff all food


# Larvae released (all groups, normalized)
p1.6wks <- ggplot(data=Collection.6wks.cumul, aes(x=Date, y=Larvae.norm, fill=TREAT)) + 
  geom_bar(stat="identity",width=.5) + ylab("Daily larvae released") +
  ggtitle("Larvae release (normalized by broodstock)\n6-week exposure") + theme_bw(base_size = 16) +   
  theme(plot.title = element_text(size = 16, hjust = 0)) + 
  scale_x_date(date_breaks = "1 week",date_labels ="%b-%d", limits=c(as.Date("2018-02-19"), as.Date("2018-03-21"))) + 
  theme(legend.position = "none") + scale_fill_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food"))
png(filename = "results/larval-release-chart-6wk.png", width = 700, height = 400)
print(Larvae.release.6wks <- p1.6wks + geom_line(data=Collection.6wks.cumul, aes(x=Date, y=cum.percap/5, group=TREAT, color=TREAT),size=.75) +
  scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food")) + scale_y_continuous(limits=c(0,120000), labels = scales::comma, sec.axis = sec_axis(~.*5,name="Total larvae released, per broodstock")))
dev.off()

# Larvae released (COLD - HIGH)
p3.6wks <- ggplot(data=subset(Collection.6wks.cumul, TREAT=="Cold - High"), aes(x=Date, y=Larvae, fill=TREAT)) + geom_bar(fill="#0571b0", stat="identity",width=.5, position = position_dodge(width=2)) + ylab("Number of Larvae Released") +
  ggtitle("Cold / High Food, 6-week exposure \nLarvae Release") + theme_bw(base_size = 14) +   
  theme(plot.title = element_text(size = 14, hjust = 0)) + 
  scale_x_date(date_breaks = "1 week",date_labels ="%b-%d", limits=c(as.Date("2018-02-19"), as.Date("2018-03-21"))) + theme(legend.position = "none")
png(filename = "results/larval-release-chart-Cold-High-6wk.png", width = 500, height = 400)
print(Larvae.release.COLD.HIGH.6wks <- p3.6wks + geom_line(data=subset(Collection.6wks.cumul, TREAT=="Cold - High"), aes(x=Date, y=cum.total/6.5, group=TREAT, color=TREAT),size=.75) + scale_color_manual(values="#0571b0") + scale_y_continuous(limits = c(0, 6e6), sec.axis = sec_axis(~.*6.5,name="Cumulative Larvae Released")))
dev.off()

# Larvae released (COLD - LOW)
p4 <- ggplot(data=subset(Collection.6wks.cumul, TREAT=="Cold - Low"), aes(x=Date, y=Larvae, fill=TREAT)) + geom_bar(fill="#92c5de", stat="identity",width=.5, position = position_dodge(width=2)) + ylab("Number of Larvae Released") +
  ggtitle("Cold / Low Food, 6-week exposure \nLarvae Release") + theme_bw(base_size = 14) +   
  theme(plot.title = element_text(size = 14, hjust = 0)) + 
  scale_x_date(date_breaks = "1 week",date_labels ="%b-%d", limits=c(as.Date("2018-02-19"), as.Date("2018-03-21"))) + 
  theme(legend.position = "none") + scale_fill_manual(values=c("seagreen3"))
png(filename = "results/larval-release-chart-Cold-Low-6wk.png", width = 500, height = 400)
print(Larvae.release.COLD.LOW <- p4+ geom_line(data=subset(Collection.6wks.cumul, TREAT=="Cold - Low"), aes(x=Date, y=cum.total/6.5, group=TREAT, color=TREAT),size=.75) + scale_color_manual(values="#92c5de") + scale_y_continuous(limits = c(0, 6e6), sec.axis = sec_axis(~.*6.5,name="Cumulative Larvae Released")))
dev.off()

# Larvae released (WARM - HIGH)
p5 <- ggplot(data=subset(Collection.6wks.cumul, TREAT=="Warm - High"), aes(x=Date, y=Larvae, fill=TREAT)) + geom_bar(fill="#ca0020", stat="identity",width=.5, position = position_dodge(width=2)) + ylab("Number of Larvae Released") +
  ggtitle("Warm / High Food, 6-week exposure \nLarvae Release") + theme_bw(base_size = 14) +   
  theme(plot.title = element_text(size = 14, hjust = 0)) + 
  scale_x_date(date_breaks = "1 week",date_labels ="%b-%d", limits=c(as.Date("2018-02-19"), as.Date("2018-03-21"))) + theme(legend.position = "none")
png(filename = "results/larval-release-chart-Warm-High-6wk.png", width = 500, height = 400)
print(Larvae.release.WARM.HIGH <- p5+ geom_line(data=subset(Collection.6wks.cumul, TREAT=="Warm - High"), aes(x=Date, y=cum.total/6.5, group=TREAT, color=TREAT),size=.75) + scale_color_manual(values="#ca0020") + scale_y_continuous(limits = c(0, 6e6), sec.axis = sec_axis(~.*6.5,name="Cumulative Larvae Released")))
dev.off()

# Larvae released (WARM - LOW)
p6 <- ggplot(data=subset(Collection.6wks.cumul, TREAT=="Warm - Low"), aes(x=Date, y=Larvae, fill=TREAT)) +  geom_bar(fill="#f4a582", stat="identity",width=.5, position = position_dodge(width=2)) + ylab("Number of Larvae Released") +
  ggtitle("Warm / Low Food, 6-week exposure \nLarvae Release") + theme_bw(base_size = 14) +   
  theme(plot.title = element_text(size = 14, hjust = 0)) + 
  scale_x_date(date_breaks = "1 week",date_labels ="%b-%d", limits=c(as.Date("2018-02-19"), as.Date("2018-03-21"))) + 
  theme(legend.position = "none")
png(filename = "results/larval-release-chart-Warm-Low-6wk.png", width = 500, height = 400)
print(Larvae.release.WARM.LOW <- p6+ geom_line(data=subset(Collection.6wks.cumul, TREAT=="Warm - Low"), aes(x=Date, y=cum.total/6.5, group=TREAT, color=TREAT),size=.75) + scale_color_manual(values="#f4a582") + scale_y_continuous(limits = c(0, 6e6), sec.axis = sec_axis(~.*6.5,name="Cumulative Larvae Released")))
dev.off()

