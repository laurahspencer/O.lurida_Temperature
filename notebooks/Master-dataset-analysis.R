
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
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ FOOD*TEMP, data=master, quasibinomial)) #not
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ Length, data=master, quasibinomial)) #not
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ Width, data=master, quasibinomial)) #not
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ I(Length*Width), data=master, quasibinomial)) #not 

# Test ALL factors in one model
Anova(glm.all <- glm(cbind(Live.35.days, Dead.35.days) ~ Live.Larvae + Length + Width + Date.stocked + TEMP + FOOD + TEMP:FOOD, data=master, quasibinomial))

# Date stocked against treatments 
Anova(glm.all <- glm(cbind(Live.35.days, Dead.35.days) ~ Date.stocked*TEMP + Date.stocked*FOOD, data=master, quasibinomial))

summary(glm.all <- glm(cbind(Live.35.days, Dead.35.days) ~ as.numeric(Date.stocked) + Live.Larvae, data=master, quasibinomial))

# size & survival plots 

ggplot(master, aes(x=1000*Length, y=100*Live.35.days/(800))) + geom_point(size=3.5, aes(color=TREAT.x)) + labs(title="% survival ~ shell length upon release", y=("% survival"), x="length upon release") + scale_color_manual(values=c("skyblue3", "seagreen3", "indianred2",  "orange1")) + theme(text = element_text(size = 16))

ggplot(master, aes(x=Date.stocked, y=100*Live.35.days/(800))) + geom_point(size=3.5, aes(color=TREAT.x)) + labs(title="% survival ~ date stocked", y=("% survival"), x="length upon release") + scale_color_manual(values=c("skyblue3", "seagreen3", "indianred2",  "orange1")) + theme(text = element_text(size = 16))

ggplotly(ggplot(master, aes(x=1000*Width, y=100*Live.35.days/(800))) + geom_point(size=3.5, aes(color=TREAT.x)) + labs(title="% survival ~ shell width upon release", y=("% survival"), x="width upon release") + scale_color_manual(values=c("skyblue3", "seagreen3", "indianred2",  "orange1")) + theme(text = element_text(size = 16)))

ggplot(master, aes(x=TREAT.x, y=Length)) + geom_boxplot(aes(fill=TREAT.x)) + labs(title="shell length ~ treatment", y=("shell length"), x="treatment") + scale_fill_manual(values=c("skyblue3", "seagreen3", "indianred2",  "orange1")) + theme(text = element_text(size = 16)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

ggplot(master, aes(x=TREAT.x, y=Width)) + geom_boxplot(aes(fill=TREAT.x)) + labs(title="shell width ~ treatment", y=("shell width"), x="treatment") + scale_fill_manual(values=c("skyblue3", "seagreen3", "indianred2",  "orange1")) + theme(text = element_text(size = 16)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

# Merge total fecundity for 7 wk and 12 wk 

# add standardized calendar day columns to each df 
treat_total.6wks.rep[17:18] <- data.frame(treat_total.6wks.rep[c("maxday", "first.big")], 
           lapply(treat_total.6wks.rep[c("maxday", "first.big")], function(x) x - min(treat_total.6wks.rep$first.big)))[3:4]

treat_total.rep[17:18] <- data.frame(treat_total.rep[c("maxday", "first.big")], 
           lapply(treat_total.rep[c("maxday", "first.big")], function(x) x - min(treat_total.rep$first.big)))[3:4]

colnames(treat_total.rep) <- colnames(treat_total.6wks.rep)
subset(master.release, trial==13) <- rbind(treat_total.rep, treat_total.6wks.rep)
subset(master.release, trial==13)$trial <- as.factor(subset(master.release, trial==13)$trial)

plot(x=subset(master.release, trial==13)$trial, y=subset(master.release, trial==13)$mean.larvae)
plot(x=subset(master.release, trial==13)$trial, y=subset(master.release, trial==13)$total.percap)
plot(x=subset(master.release, trial==13)$trial, y=subset(master.release, trial==13)$maxday.1)
plot(x=subset(master.release, trial==13)$trial, y=subset(master.release, trial==13)$perc.spawn)
plot(y=subset(master.release, trial==13)$maxday.1, x=subset(master.release, trial==13)$first.big.1,col=subset(master.release, trial==13)$TREAT)
plot(x=subset(master.release, trial==13)$release.days, y=subset(master.release, trial==13)$total.percap,col=subset(master.release, trial==13)$TREAT)

# IMPORTANT  
# Decision 10/28/2020 - don't include any data from the 6 week trial in the manuscript. Remove all 6 week data, and omit that factor ("trial") from stats 

# Date of onset differ? 
shapiro.test(subset(master.release, trial==13)$first.big.1) #can't make normal
hist(subset(master.release, trial==13)$first.big.1)
kruskal.test(first.big.1 ~ TREAT, data=subset(master.release, trial==13))
kruskal.test(first.big.1 ~ Temperature, data=subset(master.release, trial==13))
kruskal.test(first.big.1 ~ Food, data=subset(master.release, trial==13)) #YES 
plot(first.big.1 ~ Food, data=subset(master.release, trial==13)) # BIG difference 
aggregate(first.big ~ Food, data=subset(master.release, trial==13), mean) # low food went 4 days earlier on ave.

# Date of max release (after onset) difer by trial? 
shapiro.test(subset(master.release, trial==13)$maxday.1^(1/3))
summary(aov(maxday.1^(1/3) ~ Temperature*Food, data=subset(master.release, trial==13))) # YES
TukeyHSD(aov(maxday.1^(1/3) ~ Temperature*Food, data=subset(master.release, trial==13)))
plot(maxday.1~temp_food_trial, data=mutate(subset(master.release, trial==13), temp_food_trial=as.factor(paste(Temperature, Food, trial, sep="_"))), main="Peak release ~ Temp:Food:Trial")

aggregate(maxday.1 ~ Temperature*Food, data=subset(master.release, trial==13), mean)
aggregate(maxday.1 ~ Temperature*Food, data=subset(master.release, trial==13), sd)

# Mean daily larvae released differ by treatment & trial?
shapiro.test(subset(master.release, trial==13)$mean.larvae)
test <- summary(aov(mean.larvae ~ Temperature*Food, data=subset(master.release, trial==13)))

test[[1]]["F value"] #to extract F value
test[[1]]["Pr(>F)"]

# Total release norm. differ by treatment? 
shapiro.test(subset(master.release, trial==13)$total.percap)
summary(aov(total.percap ~ Temperature*Food, data=subset(master.release, trial==13))) 
aggregate(total.percap ~ Temperature*Food, data=subset(master.release, trial==13), mean)
aggregate(total.percap ~ Temperature*Food, data=subset(master.release, trial==13), sd)

# Total release (not norm.) differ by treatment & trial? 
shapiro.test(subset(master.release, trial==13)$overall_Total^(1/2))
summary(aov(overall_Total^(1/2) ~ Temperature*Food, data=subset(master.release, trial==13)))

master.release$Food <- factor(master.release$Food, levels = rev(levels(master.release$Food)))

pdf("results/total-larvae-released-points.pdf", width=5, height = 5)
ggplot(subset(master.release, trial==13), aes(x=Food:Temperature, y=total.percap)) + 
  geom_jitter(width=0.2, size=4, shape=21, aes(color=Food:Temperature), stroke=1.5) + 
  theme_bw() + 
  labs(title="Total larvae released", y=("Total released (per adult)")) + 
  scale_color_manual(values=c('#92c5de','#f4a582','#0571b0','#ca0020'), name=element_blank(), labels = c("7°C+Low-food","10°C+Low-food", "7°C+High-food", "10°C+High-food")) +
  theme(text = element_text(size = 12, color="gray25")) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  scale_y_continuous(labels=scales::comma,lim=c(30000,170000)) + 
  scale_shape_manual(values=c(8,16), name="Days in treatment", labels=c("7-week", "12-week")) 
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

# 
pdf(file = "results/broodstock-survival.pdf", width = 6, height = 2.5)
brood.mort.13
dev.off()
