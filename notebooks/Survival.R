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

# General summary stats 
summary(100*(na.omit(survival$Live.35.days)/800)) #survival to 35 days 

summary(100*(na.omit(survival$Live.50.days)/(800*3))) #survival to 50 days 

# How long was larval collection period? 
max(survival$Date.stocked)-min(survival$Date.stocked)

# Plot survival by parental treatment 

# first reverse order of food factor 
survival$FOOD <- factor(survival$FOOD, levels = rev(levels(survival$FOOD)))

pdf(file="results/jitter-survival.pdf", width =4, height = 7)
ggplot(survival, aes(x=FOOD:TEMP, y=100*Live.35.days/(800))) + 
  geom_boxplot(aes(color=FOOD:TEMP), outlier.shape = NA) + 
  geom_jitter(width=0.35, size=2, aes(color=FOOD:TEMP)) + 
  theme_bw(base_size = 12) + 
  labs(title="% survival, by parental treatment", y=("% survival"), x="Treatment replicate") + 
  scale_color_manual(values=c('#92c5de', '#f4a582', '#0571b0','#ca0020'), name=element_blank(),
                     labels = c("Cold+\nLow Food", "Warm+\nLow Food", "Cold+\nHigh Food", "Warm+\nHigh Food")) + theme(text = element_text(size = 12), legend.position = "none") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())
dev.off()

# Create new dataframe that averages survival for each larval group (both percentage, and count)
Survival.family.35 <- do.call(data.frame, aggregate(Live.35.days ~ Family+TEMP+FOOD+TRT.REP+TREAT, data = survival, FUN = function(x) c(mean = mean(x)/(800)*100, sd = sd(x)/(800)*100, mean.live = mean(x), mean.dead = (800*3)-mean(x))))

names(Survival.family.35) <- c("Family", "TEMP", "FOOD", "TRT.REP", "TREAT", "Mean.Live.35", "SD.Live.35", "mean.live", "mean.dead")

# Plot survival again, this time  one point showing average survival of each group 
ggplot(Survival.family.35, aes(x=FOOD:TEMP, y=Mean.Live.35)) + 
  geom_boxplot(aes(color=FOOD:TEMP), outlier.shape = NA) + 
  geom_jitter(width=0.35, size=2, aes(color=FOOD:TEMP)) + 
  theme_bw(base_size = 12) + 
  labs(title="Mean % survival, by parental treatment", y=("% survival"), x="Treatment") + 
  scale_color_manual(values=c('#92c5de','#f4a582','#0571b0','#ca0020'), name=element_blank(), 
                      labels = c("Cold+\nLow Food", "Warm+\nLow Food", "Cold+\nHigh Food", "Warm+\nHigh Food")) + theme(text = element_text(size = 12), legend.position = "none") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

# Save survival dataframe
write.csv(file="data/Live.35.family.csv", x=Survival.family.35)

# Examine factor survival count data from day 35 (all silos separate)
glm.Date <- glm(cbind(Live.35.days, Dead.35.days) ~ Date.stocked, data=survival, quasibinomial)
Anova(glm.Date) #sign. date stocked factor 
plot(x=survival$Date.stocked, y=survival$Live.35.days/800)
summary(lm(survival$Live.35.days/800 ~ survival$Date.stocked)) #r-sq is very small, but p-value sign. 

# Merge survival data with collection data 
survival.collect <- merge(x=survival, y=Collection[c(5,14,15,16,20)], by.x="Family", by.y="Group", all.x =TRUE, all.y=FALSE)

# NOW create a master dataframe of larval survival, larval size, and collection data to test for various effects on survival 
survival.collect$Sample.number <- as.factor(survival.collect$Sample.number) #first convert sample # to factor 

master <- merge(x=average.all[,c("Sample.number", "Length", "Width")], 
                y=survival.collect, by="Sample.number", all.x=T, all.y=T) %>% 
  mutate(Sample.number=as.numeric(Sample.number)) %>%
  mutate(ethanol=ifelse(Sample.number < 50, 'YES', 'NO')) %>%
mutate(length.cor = ifelse(Sample.number<50, Length-5.620869, Length)) 

# Run GLMS to test for effects of other factors on survival  

# Test ALL possible factors in one model
Anova(glm.all <- glm(cbind(Live.35.days, Dead.35.days) ~ FOOD + TEMP + Live.Larvae + length.cor + Date.stocked, data=master, quasibinomial))

# Temp and Food not sign.
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ FOOD*TEMP, data=master, quasibinomial)) 

# Date stocked sign.
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ Date.stocked, data=master, quasibinomial)) 
summary(glm(cbind(Live.35.days, Dead.35.days) ~ Date.stocked, data=master, quasibinomial)) #positive eatimate
ggplot(master, aes(x=Date.stocked, y=100*Live.35.days/(800))) + 
  geom_point(size=2.5) + theme_bw(base_size = 14) + 
  labs(title="% larval survival ~ Date stocked", 
       y=("% survival"), 
       x="# live larvae collected") + 
  geom_smooth(method = "auto", color="gray50")

# live larvae collected sign.
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ Live.Larvae, data=master, quasibinomial)) 
summary(glm(cbind(Live.35.days, Dead.35.days) ~ Live.Larvae, data=master, quasibinomial)) #positive estimate
ggplot(master, aes(x=Live.Larvae, y=100*Live.35.days/(800))) + 
  geom_point(size=2.5) + theme_bw(base_size = 14) + 
  labs(title="% larval survival ~ # larvae collected", 
       y=("% survival"), 
       x="# live larvae collected") + 
  geom_smooth(method = "auto", color="gray50")

# Larval length

#test using random effects model to include use of ethanol as preservation - not sign.
Anova(lme(cbind(Live.35.days, Dead.35.days) ~ Length, random=~1|ethanol, data=na.omit(master)))

#test survival ~ length using the corrected length measurement (corr. for ethanol use, see Larval Size script)
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ length.cor, data=master, quasibinomial)) # tends to be sign.
summary(glm(cbind(Live.35.days, Dead.35.days) ~ length.cor, data=master, quasibinomial)) #positive estimate
plot(100*(Live.35.days/800) ~ length.cor, data=master) #survival increases with length

# Plot % survival ~ length, by use of ethanol to see what's going on 
ggplot(master, aes(x=Length, y=100*Live.35.days/(800))) + 
  geom_point(size=2.5) + theme_bw(base_size = 14) + 
  labs(title="% larval survival ~ # larvae collected", 
       y=("% survival"), 
       x="Larval shell width") + facet_wrap(~ethanol, scales = "free") +
  geom_smooth(method = "lm", color="gray50")

# survival ~ size + parental temperature 
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ TEMP/length.cor, data=master, quasibinomial)) 
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ TEMP*length.cor, data=master, quasibinomial)) 
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ FOOD*length.cor, data=master, quasibinomial)) 
summary(glm(cbind(Live.35.days, Dead.35.days) ~ TEMP*length.cor, data=master, quasibinomial)) 

Anova(glm(cbind(Live.35.days, Dead.35.days) ~ FOOD:TEMP*length.cor-FOOD:TEMP, data=subset(master, ethanol=="NO"), quasibinomial)) 

ggplot(na.omit(subset(master, ethanol=="NO")), aes(x=length.cor, y=100*Live.35.days/(800))) + 
  geom_point(size=2.5) + theme_bw(base_size = 14) + 
  labs(title="% larval survival ~ # larvae collected", 
       y=("% survival"), 
       x="Larval shell width") + facet_wrap(~TEMP, scales = "fixed") +
  geom_smooth(method = "lm", color="gray50")

ggplot(na.omit(master), aes(x=Date.stocked, y=100*Live.35.days/(800))) + 
  geom_point(size=2.5) + theme_bw(base_size = 14) + 
  labs(title="", 
       y=("SIZE"), 
       x="DATE") + facet_wrap(~FOOD:TEMP, scales = "fixed") +
  geom_smooth(method = "lm", color="gray50")

ggplot(na.omit(master), aes(x=length.cor, y=100*Live.35.days/(800))) + 
  geom_smooth(method = "lm", color="gray70", fill="gray85") +
  geom_point(size=2.5, shape=1, stroke=1.5, 
             aes(colour=as.numeric(Date.stocked))) + 
  theme_clean(base_size = 14) + 
  labs(title="", 
       y=("Larval survival (%)"), 
       x="Larval shell width (um)") + 
  facet_wrap(~TEMP, scales = "fixed") +
  scale_colour_gradient(low="yellow", high="red") + 
  theme(legend.position = "none")

anova(lm(length.cor ~ I(100*Live.35.days/(800))*TEMP, data=master))
summary(lm(length.cor ~ I(100*Live.35.days/(800)), data=subset(master, TEMP=="Cold")))
summary(lm(length.cor ~ I(100*Live.35.days/(800)), data=subset(master, TEMP=="Warm")))

# Inspect survival ~ size by treatment, and showing ethanol use 
# Decision: there are a couple outliers in the Cold:High and Cold:low groups that influence the trend- when removed the trend goes away. Also, use of ethanol definitely influences the size (even after corrected), which confounds the survival relationship. Therefore, do not pursue. 

ggplot(na.omit(subset(master, Live.35.days<500)), aes(x=length.cor, y=Live.35.days)) + 
  geom_smooth(method = "lm", color="gray70", fill="gray85") +
  geom_point(size=2.5, stroke=1.5, 
             aes(colour=FOOD, shape=ethanol)) + 
  theme_clean(base_size = 16) + 
  labs(title="", 
       y=("Larval survival (%)"), 
       x="Larval shell width (um)") + 
  facet_wrap(~TEMP:FOOD, scales = "fixed") +
  scale_color_manual(values=c("red","blue","green","brown"))
#  theme(legend.position = "none")

#plot FOR PAPER 

# Relevel food factors such that Low comes before High
master$FOOD <- factor(master$FOOD, levels = rev(levels(master$FOOD)))

# Test correlation between survival & size and find R2 of lm  
# NOTE: only use larvae that were NOT preserved with ethanol
test <- subset(master, ethanol=="NO") %>% dplyr::select(length.cor, Live.35.days, TREAT.x) %>% 
  filter(!is.na(TREAT.x), !is.na(length.cor)) 
cor.test(x=test$length.cor, y=test$Live.35.days)
hist(test$Live.35.days^.25)
shapiro.test(test$Live.35.days^.25)
anova(lm(Live.35.days^.25 ~ length.cor, data=test)) # not sign. 
summary(lm(Live.35.days^.25 ~ length.cor, data=test))

# Plot larval survival ~ length (corrected), with fitted lm  
png(file="results/larval-surv-size.png", width = 600, height = 325)
subset(master, ethanol=="NO") %>% dplyr::select(Length, Live.35.days, TREAT.x) %>% 
  filter(!is.na(TREAT.x), !is.na(Length)) %>% 
  ggplot(aes(x=Length, y=100*Live.35.days/(800))) + 
  geom_point(size=2, color="gray40") + 
  labs(title="Larval survival ~ shell width upon release", 
       y=("% Survival"), 
       x="Shell Length (µm)") + 
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  geom_smooth(method='lm', color="gray50") +
  annotate("text", x=175, y=40, color="gray30", 
            label="Adjusted R-squared=0.012\np-value=0.162")
dev.off()

# plot larval survival by date collecte/releasd 
png(file="results/larvael-surv-time-col.png", width = 600, height = 325)
ggplot(survival, aes(x=Date.stocked, y=100*(Live.35.days/800))) + theme_bw(base_size = 14) + geom_point(size=2.5, aes(color=TREAT)) + labs(title="% Survival by treatment and date larvae released", y=("Percent Survival"), x=("Date Released")) + ylim(0,60) + scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food"))
dev.off()



# NEED TO FIX 
# Plot mean survival against mean length 
plot(x=master$length.cor, y=100*(master$Live.35.days/800), ylab="% Survival to Postset", xlab="Shell length upon release (um)", pch=16, cex=1.4, col="gray45", main="Larval survival ~ Shell length upon release")

mean(aggregate(Live.35.days ~ Family, data=survival, mean)[,2]) #average family sd = 26.2, mean = 85.9 
mean(aggregate(Live.35.days ~ TRT, data=survival, sd)[,2]) #average treatment sd = 83.0, mean = 84.9

sd(survival$Live.35.days)/mean(survival$Live.35.days)

# Extra 

# Does % survival correlate with survival variance? 
mean(Survival.family.35$SD.Live.35/Survival.family.35$Mean.Live.35) #Average CV within families
plot(x=Survival.family.35$Mean.Live.35, y=Survival.family.35$SD.Live.35)
summary(lm(Survival.family.35$SD.Live.35 ~ Survival.family.35$Mean.Live.35)) #yes, variance within family increases with survival


# Survival plot for PCSGA presentation (has different colors)
survival %>% 
  mutate(FOOD = factor(FOOD, levels=c("Low", "High"))) %>% 
  ggplot(aes(x=TEMP:FOOD, y=100*Live.35.days/(800))) + 
  geom_boxplot(aes(color=TEMP:FOOD), outlier.shape = NA) + 
  geom_jitter(width=0.35, size=1.5, aes(color=TEMP:FOOD)) + 
  theme_bw(base_size = 12) + 
  labs(title="% survival to post-set\nby parental treatment", y=("% survival"), x="Treatment replicate") + 
  scale_color_manual(values=c('darkslateblue','#0571b0','chocolate2','#ca0020'), name=element_blank()) +  
  theme(text = element_text(size = 12),
        axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") +
  ylim(c(-5, 56))

write.csv(file="data/Live.50.family.csv", x=(subset(survival, Live.50.days >0)))
Survival.family.50 <- subset(survival, Live.50.days >0)[c("Family", "TRT", "TRT.REP", "Live.50.days", "TEMP", "FOOD")]
Survival.family.50$Live.50.percent <- Survival.family.50$Live.50.days/(800*3)
write.csv(file="data/Live.50.family.csv", x=Survival.family.50)

# Model survival count data from day 50 (silos combined by family) 
glm.Date.50 <- glm(cbind(Live.50.days, Dead.50.days) ~ Date.stocked, data=survival, quasibinomial)
summary(glm.Date.50)
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


ggplot(subset(survival, TRT == "A"), aes(x=Date.stocked, y=100*(Live.35.days/800))) + geom_point(size=3, aes(color=TRT.REP)) + labs(title="Treatment A: Low Food, Low Temp\n% survival by broodstock replicate and date larvae was released", y=("Percent Survival"), x=("Date Released")) + ylim(0,100) + theme_minimal()
ggplot(subset(survival, TRT == "B"), aes(x=Date.stocked, y=100*(Live.35.days/800))) + geom_point(size=3, aes(color=TRT.REP)) + labs(title="Treatment B: Low Food, High Temp\n% survival by broodstock replicate and date larvae was released", y=("Percent Survival"), x=("Date Released"))  + ylim(0,100) + theme_minimal()
ggplot(subset(survival, TRT == "C"), aes(x=Date.stocked, y=100*(Live.35.days/800))) + geom_point(size=3, aes(color=TRT.REP)) + labs(title="Treatment C: High Food, Low Temp\n% survival by broodstock replicate and date larvae was released", y=("Percent Survival"), x=("Date Released"))  + ylim(0,100) + theme_minimal()
ggplot(subset(survival, TRT == "D"), aes(x=Date.stocked, y=100*(Live.35.days/800))) + geom_point(size=3, aes(color=TRT.REP)) + labs(title="Treatment D: High Food, High Temp\n% survival by broodstock replicate and date larvae was released", y=("Percent Survival"), x=("Date Released")) + ylim(0,100) + theme_minimal()

p5 <- ggplotly(ggplot(survival, aes(x=Date.stocked, y=100*(Live.35.days/800))) + geom_point(size=3, aes(color=TRT.REP)) + labs(title="Survival @ 5 weeks ~ Date Stocked", y=("Percent Survival"), x=("Date Released")) + ylim(0,100) + theme_minimal() + facet_wrap(~TRT))
htmlwidgets::saveWidget(as_widget(p5), file = "larval-survival-time-5weeks-trt-reps.html")

p7 <- ggplotly(ggplot(survival, aes(x=Date.stocked, y=100*(Live.50.days/(800*3)))) + geom_point(size=3, aes(color=TRT.REP)) + labs(title="Survival @ 7 weeks ~ Date Stocked", y=("Percent Survival"), x=("Date Released")) + ylim(0,100) + theme_minimal() + facet_wrap(~TRT))
htmlwidgets::saveWidget(as_widget(p7), file = "larval-survival-time-7weeks-trt-reps.html")

plot_ly(data=survival.collect, x=~Perc.live, y=~Live.50.days, type="scatter", mode="text", text=~Family)
# B2-3 only outlier - Low percent live at release, and Low survival (but not the Lowest)

# Does survival correlate with # larvae collected, aka if competition was there more mortality?
plot(x=survival.collect$Live.Larvae, y=100*(survival.collect$Live.35.days/800), pch=16, cex=1.4, col="gray45", main="% larval survival ~ # larvae collected", xlab="# live larvae collected", ylab="% survival to post-set")

png(filename = "results/larval-surv-collected.png", width = 500, height = 325)
ggplot(survival.collect, aes(x=Live.Larvae, y=100*Live.35.days/(800))) + geom_point(size=2.5, aes(color=TEMP)) + theme_bw(base_size = 14) + labs(title="% larval survival ~ # larvae collected", y=("% survival"), x="# live larvae collected") + scale_color_manual(values=c('#0571b0','#ca0020'), name=element_blank(), labels = c("Cold", "Warm")) 
dev.off()


# Possible plot for paper, not sure 
pdf(file = here::here("results", "larval-shell-length-survival.pdf"), width = 7.0, height = 4.5)
master %>% dplyr::select(length.cor, Live.35.days, TREAT.x, FOOD, TEMP) %>% 
  filter(!is.na(TREAT.x), !is.na(length.cor)) %>% 
  # mutate(TREAT.x = factor(TREAT.x, levels=c("Cold - High", "Cold - Low", "Warm - High", "Warm - Low"))) %>%
  ggplot(aes(x=length.cor, y=100*Live.35.days/(800))) + 
  geom_point(size=2.5, aes(color=FOOD:TEMP)) + 
  labs(title="Larval survival ~ shell width upon release", 
       y=("% Survival"), 
       x="Shell Width (µm)") + 
  theme_bw(base_size = 12) + 
  scale_color_manual(values=c('#92c5de','#f4a582','#0571b0','#ca0020'),
                     name=element_blank(),
                     labels = c("7°C+Low-food", "10°C+Low-food", "7°C+High-food", "10°C+High-food")) +
  theme(legend.position = c(0.15, 0.81), axis.title = element_text(color="gray25"), plot.title = element_text(color="gray25"))
dev.off()