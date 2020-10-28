
# Merge size data here 

#survival.collect$Sample.number <- as.factor(survival.collect$Sample.number) #first convert sample # to factor 
#master <- merge(x=size.surv, y=survival.collect, by.x="Sample", by.y="Sample.number") #merge 

#merge survival & collection data (survival.collect) with larval size data (average.all)
# survival.collect object created in the "Survival" script 
# average.all created in the "Larvae-length.Rmd" file

master <- merge(x=survival.collect, 
                y=average.all[,c("MaxFeret", "MinFeret", "Length", "Width", "Sample.number")], 
                by="Sample.number")
  
# Add larval length variance for each family 
master <- merge(x=master, y=larvalsize %>%
                  #mutate(sheet=as.factor(sheet)) %>%
                  group_by(Sample.number) %>%
                  summarise(length.cv = var(Length, na.rm=TRUE), 
                            width.cv = var(Width, na.rm=TRUE)), by="Sample.number")

# Test each factor separately 
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ Date.stocked, data=master, quasibinomial)) #sign.
summary(glm(cbind(Live.35.days, Dead.35.days) ~ Date.stocked, data=master, quasibinomial)) #positive eatimate
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ Live.Larvae, data=master, quasibinomial)) #tended to be sign.
summary(glm(cbind(Live.35.days, Dead.35.days) ~ Live.Larvae, data=master, quasibinomial)) #positive eatimate

Anova(glm(cbind(Live.35.days, Dead.35.days) ~ Length, data=master, quasibinomial)) #not
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ Width, data=master, quasibinomial)) #not
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ I(Length*Width), data=master, quasibinomial)) #not 
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ TEMP, data=master, quasibinomial)) #not
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ FOOD, data=master, quasibinomial)) #not

# Test ALL factors in one model, no interaction
Anova(glm.all <- glm(cbind(Live.35.days, Dead.35.days) ~ Live.Larvae + Length + Width + Date.stocked + TEMP + FOOD, data=master, quasibinomial))

# Date stocked against treatments 
Anova(glm.all <- glm(cbind(Live.35.days, Dead.35.days) ~ Date.stocked*TEMP + Date.stocked*FOOD, data=master, quasibinomial))

summary(glm.all <- glm(cbind(Live.35.days, Dead.35.days) ~ as.numeric(Date.stocked) + Live.Larvae, data=master, quasibinomial))

# size & survival plots 

ggplot(master, aes(x=1000*Length, y=100*Live.35.days/(800))) + geom_point(size=3.5, aes(color=TREAT.x)) + labs(title="% survival ~ shell length upon release", y=("% survival"), x="length upon release") + scale_color_manual(values=c("skyblue3", "seagreen3", "indianred2",  "orange1")) + theme(text = element_text(size = 16))

ggplot(master, aes(x=1000*Width, y=100*Live.35.days/(800))) + geom_point(size=3.5, aes(color=TREAT.x)) + labs(title="% survival ~ shell width upon release", y=("% survival"), x="width upon release") + scale_color_manual(values=c("skyblue3", "seagreen3", "indianred2",  "orange1")) + theme(text = element_text(size = 16))

ggplot(master, aes(x=TREAT.x, y=Length)) + geom_boxplot(aes(fill=TREAT.x)) + labs(title="shell length ~ treatment", y=("shell length"), x="treatment") + scale_fill_manual(values=c("skyblue3", "seagreen3", "indianred2",  "orange1")) + theme(text = element_text(size = 16)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

ggplot(master, aes(x=TREAT.x, y=Width)) + geom_boxplot(aes(fill=TREAT.x)) + labs(title="shell width ~ treatment", y=("shell width"), x="treatment") + scale_fill_manual(values=c("skyblue3", "seagreen3", "indianred2",  "orange1")) + theme(text = element_text(size = 16)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

# Merge total fecundity for 7 wk and 12 wk 

# add standardized calendar day columns to each df 
treat_total.6wks.rep[17:18] <- data.frame(treat_total.6wks.rep[c("maxday", "first.big")], 
           lapply(treat_total.6wks.rep[c("maxday", "first.big")], function(x) x - min(treat_total.6wks.rep$first.big)))[3:4]

treat_total.rep[17:18] <- data.frame(treat_total.rep[c("maxday", "first.big")], 
           lapply(treat_total.rep[c("maxday", "first.big")], function(x) x - min(treat_total.rep$first.big)))[3:4]

colnames(treat_total.rep) <- colnames(treat_total.6wks.rep)
master.release <- rbind(treat_total.rep, treat_total.6wks.rep)
master.release$trial <- as.factor(master.release$trial)

plot(x=master.release$trial, y=master.release$mean.larvae)
plot(x=master.release$trial, y=master.release$total.percap)
plot(x=master.release$trial, y=master.release$maxday.1)
plot(x=master.release$trial, y=master.release$perc.spawn)
plot(y=master.release$maxday.1, x=master.release$first.big.1,col=master.release$TREAT)
plot(x=master.release$release.days, y=master.release$total.percap,col=master.release$TREAT)

# Total release norm. differ by treatment & trial? 
shapiro.test(master.release$total.percap)
summary(aov(total.percap ~ trial*Temperature*Food, data=master.release)) #test trial + treats
aggregate(total.percap ~ trial, data=master.release, mean)
aggregate(total.percap ~ trial, data=master.release, sd)
164025.23/92243.35

# Total release (not norm.) differ by treatment & trial? 
shapiro.test(master.release$overall_Total^(1/2))
summary(aov(overall_Total^(1/2) ~ trial*Temperature*Food, data=master.release))

# Mean daily larvae released differ by treatment & trial?
shapiro.test(master.release$mean.larvae)
summary(aov(mean.larvae ~ trial*Temperature*Food, data=master.release))

# Date of max release (after onset) difer by trial? 
shapiro.test(master.release$maxday.1^(1/3))
summary(aov(maxday.1^(1/3) ~ trial*Temperature*Food, data=master.release))
TukeyHSD(aov(maxday.1^(1/3) ~ trial*Temperature*Food, data=master.release))
plot(x=master.release$trial, y=master.release$maxday.1)
plot(maxday.1~temp_food_trial, data=mutate(master.release, temp_food_trial=as.factor(paste(Temperature, Food, trial, sep="_"))), main="Peak release ~ Temp:Food:Trial")


aggregate(maxday.1 ~ trial, data=master.release, mean)
aggregate(maxday.1 ~ trial, data=master.release, sd)

# Date of onset differ by trial? 
shapiro.test(master.release$first.big.1) #can't make normal
hist(master.release$first.big.1)
kruskal.test(first.big.1 ~ trial, data=master.release)
kruskal.test(first.big.1 ~ Temperature, data=master.release)
kruskal.test(first.big.1 ~ Food, data=master.release)
plot(first.big.1 ~ Food, data=master.release)
aggregate(first.big ~ Food, data=master.release, mean)
72.85714-74.92857
kruskal.test(first.big.1 ~ TREAT, data=master.release)
kruskal.test(first.big.1 ~ temp_trial, data=mutate(master.release, temp_trial=as.factor(paste(Temperature, trial, sep="_"))))
kruskal.test(first.big.1 ~ food_trial, data=mutate(master.release, food_trial=as.factor(paste(Food, trial, sep="_"))))

# Plot to assess pairwise differences 
plot(first.big.1~food_trial, data=mutate(master.release, food_trial=as.factor(paste(Food, trial, sep="_"))), main="Release onset ~ Food:Trial")
aggregate(first.big.1 ~ food_trial, data=mutate(master.release, food_trial=as.factor(paste(Food, trial, sep="_"))), sd)
4.375-0.625

kruskal.test(first.big.1 ~ temp_food_trial, data=mutate(master.release, temp_food_trial=as.factor(paste(Temperature, Food, trial, sep="_"))))
plot(first.big.1~temp_food_trial, data=mutate(master.release, temp_food_trial=as.factor(paste(Temperature, Food, trial, sep="_"))), main="Release onset ~ Temp:Food:Trial")


# Perc. spawn as female differ by trial? 
shapiro.test(master.release$perc.spawn)
hist(master.release$perc.spawn)
summary(aov(perc.spawn ~ trial*Temperature*Food, data=master.release))
plot(x=master.release$trial, y=master.release$perc.spawn)

master.release$trial <- factor(master.release$trial, levels = rev(levels(master.release$trial)))

pdf("results/total-larvae-released-points-both-trials.pdf", width=5, height = 5)
ggplot(master.release, aes(x=TREAT, y=total.percap)) + geom_jitter(width=0.2, size=5, aes(color=TREAT, shape=trial)) + theme_bw() + labs(title="Total larvae released\nby treatment and exposure time", y=("Total released (per broodstock)")) + scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("7°C+high-food", "7°C+low-food", "10°C+high-food", "10°C+low-food")) + theme(text = element_text(size = 12)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) + scale_y_continuous(labels=scales::comma,lim=c(30000,225000)) + scale_shape_manual(values=c(8,16), name="Days in treatment", labels=c("7-week", "12-week")) 
dev.off()

#larval.length.release <- merge(x=survival.collect[,c("Sample.number", "Live.Larvae")], y=do.call(data.frame, aggregate(mm ~ Sample + TREAT + TEMP + FOOD, subset(new.length.ann, Length.Width=="length"), FUN = function(x) c(mean = mean(x), sd = sd(x), cv=100*sd(x)/mean(x)))), by.x="Sample.number", by.y="Sample")


ggplot(master, aes(x=Live.Larvae, y=length.cv)) +geom_jitter(width=0.2, size=5, aes(color=TREAT.x)) + theme_bw() + labs(title="Larval length variance within sample ~ Total larvae released", y=("Larval length variance (in sample)"), x="# Larvae collected from spawning tank") + scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food")) + theme(text = element_text(size = 12))

ggplot(master, aes(x=TREAT.x, y=length.cv)) +geom_boxplot(aes(color=TREAT.x)) + theme_bw() + labs(title="CV (%) Larval length within sample", y=("Larval length CV (in sample)"), x=element_blank()) + scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food")) + theme(text = element_text(size = 12))+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

ggplot(master, aes(x=length.cv, y=100*(Live.35.days/800))) +geom_jitter(width=0.2, size=5, aes(color=TREAT.x)) + theme_bw() + labs(title="Larval survival ~ Variance in larval length within sample", y=("% Survival to post-settlement"), x="Larval length variance (in sample)") + scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food")) + theme(text = element_text(size = 12))

larval.width.release <- merge(x=survival.collect[,c("Sample.number", "Live.Larvae")], y=do.call(data.frame, aggregate(mm ~ Sample + TREAT + TEMP + FOOD, subset(new.length.ann, Length.Width=="width"), FUN = function(x) c(mean = mean(x), sd = sd(x), cv=100*sd(x)/mean(x)))), by.x="Sample.number", by.y="Sample")

ggplot(master, aes(x=Live.Larvae, y=width.cv)) +geom_jitter(width=0.2, size=5, aes(color=TREAT.x)) + theme_bw() + labs(title="CV (%) Larval width within sample ~ Total larvae released", y=("Larval width CV (in sample)"), x="# Larvae collected from spawning tank") + scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food")) + theme(text = element_text(size = 12)) 

ggplot(master, aes(x=Live.Larvae, y=100*(Live.35.days/800))) +geom_jitter(width=0.2, size=5, aes(color=TREAT.x)) + theme_bw() + labs(title="% Survival by treatment and no. larvae released in group", y=("% Survival through metamorphosis"), x="No. larvae released") + scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food")) + theme(text = element_text(size = 12)) 

ggplot(master, aes(x=Date.stocked, y=100*(Live.35.days/800))) +geom_jitter(width=0.2, size=5, aes(color=TREAT.x)) + theme_bw() + labs(title="% Survival by treatment and date larvae were released", y=("% Survival through metamorphosis"), x="Date released") + scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food")) + theme(text = element_text(size = 12)) 

ggplot(master, aes(x=width.cv, y=100*(Live.35.days/800))) +geom_jitter(width=0.2, size=5, aes(color=TREAT.x)) + theme_bw() + labs(title="Larval survival ~ CV (%) Larval width within sample", y=("% Survival to post-settlement"), x="Larval wkidth CV (in sample)") + scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food")) + theme(text = element_text(size = 12))

# Broodstock mortality, both trials combined

brood.mortality.all <- read.csv("data/broodstock-mortality-all-trials.csv", header=T, stringsAsFactors = F, na.strings = "NA")
brood.mortality.all$Date <- as.Date(brood.mortality.all$Date, format = "%m/%d/%y")
brood.mortality.all <- melt(data = brood.mortality.all, id.vars = "Date", value.name = "Alive", variable.name = "TREAT.trial")
brood.mortality.all$TREAT <- brood.mortality.all$TREAT.trial
brood.mortality.all$TREAT <- factor(gsub(".6|.13", "", brood.mortality.all$TREAT), levels = c("C", "A", "D", "B"))
brood.mortality.all$TREAT.trial <- as.factor(gsub("A.|B.|C.|D.", "", brood.mortality.all$TREAT.trial))
brood.mortality.all$Alive <- as.numeric(brood.mortality.all$Alive)

brood.mort.6 <-  ggplot(data=subset(brood.mortality.all, TREAT.trial=="6" & Alive!="NA"), aes(x=Date, y=100*Alive, group=TREAT, col=TREAT)) + 
  theme_bw(base_size = 11) + 
  geom_step() + geom_point(shape=16) + 
  ggtitle("Broodstock survival over time") + 
  scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank()) + 
  scale_x_date(date_breaks = "3 week",date_labels ="%b %d", limits = c(min(brood.mortality.all$Date), max(brood.mortality.all$Date))) +
  scale_y_continuous(limits=c(50, 100)) +
  geom_vline(xintercept = as.numeric(as.Date("2018-01-24")), linetype="dashed", color = "gray50", size=.5) + 
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

brood.mort.13 <- ggplot(data=subset(brood.mortality.all, TREAT.trial=="13" & Alive!="NA"), aes(x=Date, y=100*Alive, group=TREAT, col=TREAT)) + 
  theme_bw(base_size = 12) + 
  geom_step() + geom_point(shape=16) + xlab("Date") + ylab("% survival") + 
  scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank()) + 
  scale_x_date(date_breaks = "3 week",date_labels ="%b %d", limits = c(min(brood.mortality.all$Date), max(brood.mortality.all$Date))) +
  scale_y_continuous(limits=c(50, 100)) +
  geom_vline(xintercept = as.numeric(as.Date("2018-02-28")), linetype="dashed", color = "gray50", size=.5) +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank())

pdf(file = "results/broodstock-survival.pdf", width = 6, height = 5)
grid.arrange(brood.mort.6, brood.mort.13)
dev.off()

