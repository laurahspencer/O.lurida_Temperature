
survival <- read.csv("data/Survival.csv", header=T, na.strings = "NA", stringsAsFactors = F, colClasses=
                       c(rep("character", times=2), "numeric", rep("factor", times=4), rep("numeric", times=4), "character", rep("numeric", times=5), rep("character", times=6)))
survival <- survival[c(-18:-23)]
survival$Date.stocked <- as.Date(survival$Date.stocked, format = "%m/%d/%y")
survival$Date.initial.count <- as.Date(survival$Date.initial.count, format = "%m/%d/%y")
survival$Date.imaged <- as.Date(survival$Date.imaged, format = "%m/%d/%y")
survival$TEMP <- survival$TRT
survival$TEMP <- gsub("\\<A\\>|\\<C\\>", "Cold", survival$TEMP)
survival$TEMP <- gsub("\\<B\\>|\\<D\\>", "Warm", survival$TEMP)
survival$FOOD <- survival$TRT
survival$FOOD <- gsub("\\<A\\>|\\<B\\>", "Low", survival$FOOD)
survival$FOOD <- gsub("\\<C\\>|\\<D\\>", "High", survival$FOOD)
survival$FOOD <- as.factor(survival$FOOD)
survival$TEMP <- as.factor(survival$TEMP)
survival$Dead.50.days <- (3*800)-survival$Live.50.days
survival$Dead.35.days <- 800-survival$Live.35.days
survival$TREAT <- as.factor(paste(survival$TEMP, "-", survival$FOOD))

100*mean(na.omit(survival$Live.50.days)/(800*3))
summary(100*(na.omit(survival$Live.50.days)/(800*3)))

100*mean(na.omit(survival$Live.35.days)/800)
100*median(na.omit(survival$Live.35.days)/800)
summary(100*(na.omit(survival$Live.35.days)/800))

plot(Live.50.days/(3*800) ~ Date.stocked, data=survival, col=TREAT, pch=8)

jpeg(file="results/boxplot-survival.jpeg", width = 900, height = 600)
plot(Live.50.days/(3*800) ~ TREAT, data=survival, col=c("skyblue3", "seagreen3", "indianred2",  "orange1"), main="Mean % survival across 12 groups per treatment", xlab="Treatment", ylab="Mean % survival", cex.lab=1.5, cex.main=1.5, par(mar=c(5,5,4.1,2.1)))
dev.off()

max(survival$Date.stocked)-min(survival$Date.stocked)

#levels(survival$TREAT) #color order ="skyblue3", "seagreen3", "indianred2",  "orange1"

#levels(survival$TRT.REP) <- c("A1\nCold-Low", "A2\nCold-Low",     # A1 and A2
                            #  "B1\nWarm-Low", "B2\nWarm-Low",     # B1 and B2
                            #  "C1\nCold-High", "C2\nCold-High",   # C1 and C2
                            #  "D1\nWarm-High", "D2\nWarm-High")   # D1 and D2

#levels(survival$TRT.REP) # color order: "seagreen3",  "orange1", skyblue3", "indianred2"
jpeg(file="results/boxplot-survival-rep.jpeg", width = 900, height = 600)
plot(Live.35.days/(800) ~ TRT.REP, data=survival, col=c("seagreen3", "seagreen3", "orange1","orange1", "skyblue3", "skyblue3", "indianred2",  "indianred2"), main="Mean % survival between treatment rep\n6 groups per rep", xlab="Treatment", ylab="Mean % survival", cex.lab=1.6, cex.main=1.5, par(mar=c(5,5,4.1,2.1)))
dev.off()

pdf(file="results/jitter-survival.pdf", width =7, height = 4.75)
ggplot(survival, aes(x=TREAT, y=100*Live.35.days/(800))) + 
  geom_violin(aes(color=TREAT)) + 
  geom_jitter(width=0.35, size=2.5, aes(color=TREAT)) + 
  theme_bw(base_size = 14) + 
  labs(title="% survival, by parental treatment", y=("% survival"), x="Treatment replicate") + 
  scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food")) + theme(text = element_text(size = 16)) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())
dev.off()


100*survival$Live.35.days/800

Survival.family.35 <- do.call(data.frame, aggregate(Live.35.days ~ Family+TEMP+FOOD+TRT.REP+TREAT, data = survival, FUN = function(x) c(mean = mean(x)/(800)*100, sd = sd(x)/(800)*100)))
names(Survival.family.35) <- c("Family", "TEMP", "FOOD", "TRT.REP", "TREAT", "Mean.Live.35", "SD.Live.35")

mean(Survival.family.35$SD.Live.35/Survival.family.35$Mean.Live.35) #CV within families
plot(x=Survival.family.35$Mean.Live.35, y=Survival.family.35$SD.Live.35)
summary(lm(Survival.family.35$SD.Live.35 ~ Survival.family.35$Mean.Live.35))

write.csv(file="data/Live.35.family.csv", x=Survival.family.35)

write.csv(file="data/Live.50.family.csv", x=(subset(survival, Live.50.days >0)))
 
Survival.family.50 <- subset(survival, Live.50.days >0)[c("Family", "TRT", "TRT.REP", "Live.50.days", "TEMP", "FOOD")]
Survival.family.50$Live.50.percent <- Survival.family.50$Live.50.days/(800*3)
write.csv(file="data/Live.50.family.csv", x=Survival.family.50)

# Model survival count data from day 35 (all silos separate)
glm.Date <- glm(cbind(Live.35.days, Dead.35.days) ~ Date.stocked, data=survival, quasibinomial)
Anova(glm.Date) #sign. date stocked factor 
plot(x=survival$Date.stocked, y=survival$Live.35.days/800)
summary(lm(survival$Live.35.days/800 ~ survival$Date.stocked)) #r-sq is very small, but p-value sign. 

summary(glm.TEMP <- glm(cbind(Live.35.days, Dead.35.days) ~ TEMP, data=survival, quasibinomial))
Anova(glm.TEMP)
summary(glm.FOOD <- glm(cbind(Live.35.days, Dead.35.days) ~ FOOD, data=survival, quasibinomial))
Anova(glm.FOOD)
summary(glm.TREAT <- glm(cbind(Live.35.days, Dead.35.days) ~ FOOD*TEMP, data=survival, quasibinomial))
Anova(glm.TREAT)
summary(glm.TREAT.Rep <- glm(cbind(Live.35.days, Dead.35.days) ~ TRT.REP, data=survival, quasibinomial))
Anova(glm.TREAT.Rep)

# Run mixed effects model, to include broodstock treatment rep as a random variable 

summary(glmer <- glmer(cbind(Live.35.days, Dead.35.days) ~ FOOD*TEMP + (1|TRT.REP), data=survival, family=binomial))
Anova(glmer)

# Pairwise comparison 

summary(glht(glm.TREAT.Rep, mcp(TRT.REP="Tukey")))
?glht

# Model survival count data from day 50 (silos combined by family) 
glm.Date.50 <- glm(cbind(Live.50.days, Dead.50.days) ~ Date.stocked, data=survival, quasibinomial)
summary(glm.Date.50)
glm.TEMP.50 <- glm(cbind(Live.50.days, Dead.50.days) ~ TEMP, data=survival, quasibinomial)
summary(glm.TEMP.50)
glm.FOOD.50 <- glm(cbind(Live.50.days, Dead.50.days) ~ FOOD, data=survival, quasibinomial)
summary(glm.FOOD.50)
glm.TREAT.50 <- glm(cbind(Live.50.days, Dead.50.days) ~ FOOD*TEMP, data=survival, quasibinomial)
summary(glm.TREAT.50)

# Plot - color by temp, symbol by food 
jpeg(file="results/survival-time-treat.jpeg", width = 700, height = 500)
ggplot(survival, aes(x=Date.stocked, y=100*(Live.35.days/800))) + geom_point(size=5, aes(color=TEMP, shape=FOOD)) + labs(title="% Survival by treatment and date larvae was released", y=("Percent Survival"), x=("Date Released")) + ylim(0,60) + scale_color_manual(values=c("royalblue2", "tomato1"))+ scale_shape_manual(values=c(16, 8)) + theme(text = element_text(size = 20))
dev.off()

plot_ly(data = survival, x = ~Date.stocked, y = ~100*(Live.50.days/(800*3)), type="scatter", mode="markers", marker=list(size=14), symbol=~FOOD, color=~TEMP, colors = c("royalblue2", "tomato1"), hovertext=~TRT.REP) %>%  #generate plotly plot
  layout(title="2018 O. lurida Pre-Conditioning Experiment\nPercent Survival",
         yaxis = list(title = '% Survival', size=26), xaxis=list(title="Larval release date", size=30),
         legend = list(x=.05, y=.95))

# plot each group separately 
ggplot(subset(survival, TRT == "A"), aes(x=Date.stocked, y=100*(Live.35.days/800))) + geom_point(size=3, aes(color=TRT.REP)) + labs(title="Treatment A: Low Food, Low Temp\n% survival by broodstock replicate and date larvae was released", y=("Percent Survival"), x=("Date Released")) + ylim(0,100)
ggplot(subset(survival, TRT == "B"), aes(x=Date.stocked, y=100*(Live.35.days/800))) + geom_point(size=3, aes(color=TRT.REP)) + labs(title="Treatment B: Low Food, High Temp\n% survival by broodstock replicate and date larvae was released", y=("Percent Survival"), x=("Date Released"))  + ylim(0,100)
ggplot(subset(survival, TRT == "C"), aes(x=Date.stocked, y=100*(Live.35.days/800))) + geom_point(size=3, aes(color=TRT.REP)) + labs(title="Treatment C: High Food, Low Temp\n% survival by broodstock replicate and date larvae was released", y=("Percent Survival"), x=("Date Released"))  + ylim(0,100)
ggplot(subset(survival, TRT == "D"), aes(x=Date.stocked, y=100*(Live.35.days/800))) + geom_point(size=3, aes(color=TRT.REP)) + labs(title="Treatment D: High Food, High Temp\n% survival by broodstock replicate and date larvae was released", y=("Percent Survival"), x=("Date Released")) + ylim(0,100) 

# Does survival correlate with %live/dead in newly released larvae? 
survival.collect <- merge(x=survival, y=Collection[c(5,14,15,16,20)], by.x="Family", by.y="Group", all.x =TRUE, all.y=FALSE)
plot_ly(data=survival.collect, x=~Perc.live, y=~Live.50.days, type="scatter", mode="text", text=~Family)
# B2-3 only outlier - Low percent live at release, and Low survival (but not the Lowest)

# Does survival correlate with # larvae collected, aka if competition was there more mortality?
plot(x=survival.collect$Live.Larvae, y=100*(survival.collect$Live.35.days/800), pch=16, cex=1.4, col="gray45", main="% larval survival ~ # larvae collected", xlab="# live larvae collected", ylab="% survival to post-set")

png(filename = "results/larval-surv-collected.png", width = 500, height = 325)
ggplot(survival.collect, aes(x=Live.Larvae, y=100*Live.35.days/(800))) + geom_point(size=2.5, aes(color=TEMP)) + theme_bw(base_size = 14) + labs(title="% larval survival ~ # larvae collected", y=("% survival"), x="# live larvae collected") + scale_color_manual(values=c('#0571b0','#ca0020'), name=element_blank(), labels = c("Cold", "Warm")) 
dev.off()

png(file="results/larvael-surv-time-col.png", width = 600, height = 325)
ggplot(subset(survival, TRT=="D"), aes(x=Date.stocked, y=100*(Live.35.days/800))) + theme_bw(base_size = 14) + geom_point(size=2.5, aes(color=TREAT)) + labs(title="% Survival by treatment and date larvae was released", y=("Percent Survival"), x=("Date Released")) + ylim(0,60) #+ scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food"))
dev.off()



# Run GLM 
summary(glm.collect <- glm(cbind(Live.35.days, Dead.35.days) ~ Live.Larvae, data=survival.collect, quasibinomial))
Anova(glm.collect) # yes a factor 

summary(glm.collect.dead <- glm(cbind(Live.35.days, Dead.35.days) ~ Dead.Larvae, data=survival.collect, quasibinomial))
Anova(glm.collect.dead) #not a factor 


# One big model with all factors / interactions tested 

summary(glm.collect.trt <- glm(cbind(Live.35.days, Dead.35.days) ~ Live.Larvae*Date.stocked*TEMP*FOOD, data=survival.collect, quasibinomial))
Anova(glm.collect.trt) 
plot(glm.collect.trt)

# Create a master dataframe, merged with larval size, to test effects on survival 
survival.collect$Sample.number <- as.factor(survival.collect$Sample.number) #first convert sample # to factor 
master <- merge(x=average.all[,c("Sample.number", "Length", "Width")], y=survival.collect, by="Sample.number", all.x=T, all.y=T)  #merge 

# Test each factor separately 
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ Date.stocked, data=master, quasibinomial)) # Date stocked sign.
summary(glm(cbind(Live.35.days, Dead.35.days) ~ Date.stocked, data=master, quasibinomial)) #positive eatimate
plot(100*(Live.35.days/800) ~ Date.stocked, data=master) #survival increases with time

Anova(glm(cbind(Live.35.days, Dead.35.days) ~ Live.Larvae, data=master, quasibinomial)) # live larvae collected sign.
summary(glm(cbind(Live.35.days, Dead.35.days) ~ Live.Larvae, data=master, quasibinomial)) #positive estimate
plot(100*(Live.35.days/800) ~ Live.Larvae, data=master) #survival increases with time

Anova(glm(cbind(Live.35.days, Dead.35.days) ~ Length, data=master, quasibinomial)) #not
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ Width, data=master, quasibinomial)) #not
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ I(Length*Width), data=master, quasibinomial)) #not 
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ TEMP, data=master, quasibinomial)) #not
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ FOOD, data=master, quasibinomial)) #not

# Test ALL factors in one model, no interaction
Anova(glm.all <- glm(cbind(Live.35.days, Dead.35.days) ~ Live.Larvae + Length + Width + Date.stocked + TEMP + FOOD, data=master, quasibinomial))

# Date stocked against treatments 
Anova(glm.all <- glm(cbind(Live.35.days, Dead.35.days) ~ Date.stocked*TEMP + Date.stocked*FOOD, data=master, quasibinomial))

Anova(glm(cbind(Live.35.days, Dead.35.days) ~ Live.Larvae*Date.stocked*TREAT.x, data=master, quasibinomial))
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ Live.Larvae + Date.stocked + TREAT.x:Live.Larvae + TREAT.x:Date.stocked, data=master, quasibinomial))

summary(glm(cbind(Live.35.days, Dead.35.days) ~ Live.Larvae:TREAT.x, data=master, quasibinomial))
summary(glm(cbind(Live.35.days, Dead.35.days) ~ Date.stocked:TREAT.x, data=master, quasibinomial))

# Do pairwise comparison to see if survival differed among treatments 
# COLD - High
# COLD - Low
# Warm - High
# Warm - Low

master$TREAT.x <- relevel(master$TREAT.x, ref = "Warm - High")
summary(glm(cbind(Live.35.days, Dead.35.days) ~ TREAT.x, data=master, quasibinomial))

# Does date stocked influence survival in both temperatures? (from plots, looks like just cold)
summary(glm(cbind(Live.35.days, Dead.35.days) ~ Date.stocked*FOOD, data=subset(master, TEMP=="Cold"), quasibinomial))
summary(glm(cbind(Live.35.days, Dead.35.days) ~ Date.stocked*FOOD, data=subset(master, TEMP=="Warm"), quasibinomial))

summary(glm.all <- glm(cbind(Live.35.days, Dead.35.days) ~ as.numeric(Date.stocked) + Live.Larvae, data=master, quasibinomial))

# Test these again, with temps separately 
summary(glm.all <- glm(cbind(Live.35.days, Dead.35.days) ~ as.numeric(Date.stocked) + Live.Larvae, data=subset(master, TEMP=="Cold"), quasibinomial))

summary(glm.all <- glm(cbind(Live.35.days, Dead.35.days) ~ as.numeric(Date.stocked) + Live.Larvae, data=subset(master, TEMP=="Warm"), quasibinomial)) 

# size ~ survival plots 
pdf(file = here::here("results", "larval-shell-length-survival.pdf"), width = 8.0, height = 5.25)
master %>% dplyr::select(Length, Live.35.days, TREAT.x) %>% 
  filter(!is.na(TREAT.x), !is.na(Length)) %>% 
  mutate(TREAT.x = factor(TREAT.x, levels=c("Cold - High", "Cold - Low", "Warm - High", "Warm - Low"))) %>%
ggplot(aes(x=Length, y=100*Live.35.days/(800))) + 
  geom_point(size=2.5, aes(color=TREAT.x)) + 
  labs(title="Larval survival ~ shell width upon release", 
       y=("% Survival"), 
       x="Shell Width (µm)") + 
  theme_bw(base_size = 14) + 
  scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'),
      name=element_blank(), 
      labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food"))
dev.off()

png(filename = "results/larval-shell-width-survival.png", width = 600, height = 325)
ggplot(master, aes(x=Width, y=100*Live.35.days/(800))) + geom_point(size=2.5, aes(color=TREAT.x)) + labs(title="Larval survival ~ shell height upon release", y=("% survival"), x="Shell height (µm)", size=14) + theme_bw(base_size = 14) + scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food"))
dev.off()

# NEED TO FIX 
# Plot mean survival against mean length 
plot(x=master$length, y=100*(master$Live.35.days/800), ylab="% Survival to Postset", xlab="Shell length upon release (um)", pch=16, cex=1.4, col="gray45", main="Larval survival ~ Shell length upon release")

# Plot mean survival against mean width 
plot(x=master$width, y=100*(master$Live.35.days/800), ylab="% Survival to Postset", xlab="Shell width upon release (um)", pch=16, cex=1.4, col="gray45", main="Larval survival ~ Shell width upon release")

str(survival)
mean(aggregate(Live.35.days ~ Family, data=survival, mean)[,2]) #average family sd = 26.2, mean = 85.9 
mean(aggregate(Live.35.days ~ TRT, data=survival, sd)[,2]) #average treatment sd = 83.0, mean = 84.9

sd(survival$Live.35.days)/mean(survival$Live.35.days)
