# ----- import size data and asses by treatment  

new.length <- read.csv("data/New_larvae-measurements.csv", header=T, stringsAsFactors = F)[c(1:7)]
new.length$Sample <- as.factor(new.length$Sample)
new.length.ann <- merge(x=new.length, y=Collection[,c("Group", "Sample.number")], by.x="Sample", by.y="Sample.number", all.x=T, all.y=F)
new.length.ann <- merge(x=new.length.ann[,c("Sample", "mm", "Length.Width", "Group")],y=survival[!duplicated(survival$Family),][c("Date.stocked", "Family", "TRT", "TRT.REP", "TEMP", "FOOD", "TREAT")], by.x="Group", by.y="Family")

# Plot length by treatment 
ggplot(subset(new.length.ann, Length.Width=="length"), aes(x=TREAT, y=1000*mm)) + geom_boxplot(aes(fill=TREAT)) + labs(title="shell length ~ treatment", y=("shell length"), x="treatment") + scale_fill_manual(values=c("skyblue3", "seagreen3", "indianred2",  "orange1")) + theme(text = element_text(size = 16)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

# Plot width by treatment 
ggplot(subset(new.length.ann, Length.Width=="width"), aes(x=TREAT, y=1000*mm)) + geom_boxplot(aes(fill=TREAT)) + labs(title="shell width ~ treatment", y=("shell width"), x="treatment") + scale_fill_manual(values=c("skyblue3", "seagreen3", "indianred2",  "orange1")) + theme(text = element_text(size = 16)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

summary(aov(mm^2 ~ TEMP*FOOD, data=subset(new.length.ann, Length.Width=="width")))
summary(aov(mm^2 ~ TEMP, data=subset(new.length.ann, Length.Width=="width"))) # TEMPERATURE AFFECTS WIDTH (tested on all data)

# length is non-normal, can't transform. Use non-parametric 
shapiro.test(log(subset(new.length.ann, Length.Width=="length")$mm)^2)
summary(aov(mm^2 ~ TEMP*FOOD, data=subset(new.length.ann, Length.Width=="length")))
summary(aov(mm^2 ~ TEMP, data=subset(new.length.ann, Length.Width=="length"))) # TEMPERATURE AFFECTS WIDTH 
kruskal.test(subset(new.length.ann, Length.Width=="length" & TEMP=="COLD")$mm, g = subset(new.length.ann, Length.Width=="length" & TEMP=="COLD")$FOOD) #no effect of food within cold group 
kruskal.test(subset(new.length.ann, Length.Width=="length" & TEMP=="WARM")$mm, g = subset(new.length.ann, Length.Width=="length" & TEMP=="WARM")$FOOD) #Food had an effect within warm group 
kruskal.test(subset(new.length.ann, Length.Width=="length" & FOOD=="HIGH")$mm, g = subset(new.length.ann, Length.Width=="length" & FOOD=="HIGH")$TEMP) #temp had an effect within high food group 
kruskal.test(subset(new.length.ann, Length.Width=="length" & FOOD=="LOW")$mm, g = subset(new.length.ann, Length.Width=="length" & FOOD=="LOW")$TEMP) #temp had an effect within high food group 




# Did parental treatment influence larval size upon release? 

# length 

# test anova assumptions - remove duplicate entries for each sample 
shapiro.test(master[!duplicated(master$Sample), ]$length)
bartlett.test(x=master[!duplicated(master$Sample), ]$length, g=master[!duplicated(master$Sample), ]$TRT)

anova(test <- lm(length ~ FOOD*TEMP, data=master[!duplicated(master$Sample), ]))
anova(test <- lm(length ~ TEMP, data=master[!duplicated(master$Sample), ])) 
TukeyHSD(aov(length ~ TREAT.x, data=master[!duplicated(master$Sample), ]))   

# width 

#test anova assumptions 
shapiro.test(master[!duplicated(master$Sample), ]$width)
bartlett.test(x=master[!duplicated(master$Sample), ]$width, g=master[!duplicated(master$Sample), ]$TRT)

anova(lm(width ~ FOOD*TEMP, data=master[!duplicated(master$Sample), ])) # Temperature influenced shell width  
TukeyHSD(aov(width ~ FOOD*TEMP, data=master[!duplicated(master$Sample), ]))  

#---------- Mean larval size by group, test whether differences are retained in mean data  

length.mean <- aggregate(mm ~ Sample + Length.Width + TEMP + FOOD + TREAT, new.length.ann, mean)
shapiro.test(subset(length.mean, Length.Width == "length")$mm)
shapiro.test(subset(length.mean, Length.Width == "width")$mm)
bartlett.test(subset(length.mean, Length.Width == "length")$mm, g=subset(length.mean, Length.Width == "length")$TREAT)
bartlett.test(subset(length.mean, Length.Width == "width")$mm, g=subset(length.mean, Length.Width == "width")$TREAT)

summary(aov(mm ~ FOOD*TEMP, data=subset(length.mean, Length.Width == "length")))  # no diff  
TukeyHSD(aov(mm ~ FOOD*TEMP, data=subset(length.mean, Length.Width == "length")))  # no diff  
summary(aov(mm ~ FOOD*TEMP, data=subset(length.mean, Length.Width == "width")))   # no diff 

# size range, all 
summary(subset(length.mean, Length.Width=="length")$mm)
summary(subset(length.mean, Length.Width=="width")$mm)

# mean size, by treatment 
aggregate(mm ~ FOOD + TEMP, mean, data=subset(length.mean, Length.Width=="length"))
aggregate(mm ~ FOOD + TEMP, mean, data=subset(length.mean, Length.Width=="width"))


# Plot mean width by treatment 
png(file="results/larval-shell-width.png", width=300, height=300)
ggplot(subset(length.mean, Length.Width=="width"), aes(x=TREAT, y=1000*mm)) + geom_boxplot(lwd=0.8, aes(col=TREAT)) + theme_bw(base_size = 11) + labs(title="Larval shell width\nby parental treatment", y=("shell width (µm)"), x="treatment") + scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food")) + theme(text = element_text(size = 12)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) + geom_jitter(width=0.2, aes(col=TREAT))
dev.off()

# Plot mean length by treatment 
png(file="results/larval-shell-length.png", width=300, height=300)
ggplot(subset(length.mean, Length.Width=="length"), aes(x=TREAT, y=1000*mm)) + geom_boxplot(lwd=0.8, aes(col=TREAT)) + theme_bw(base_size = 11) + labs(title="Larval shell length\nby parental treatment", y=("shell length (µm)"), x="treatment") + scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food")) + theme(text = element_text(size = 12)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) + geom_jitter(width=0.2, aes(col=TREAT))
dev.off()
