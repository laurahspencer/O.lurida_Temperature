library(dplyr)
library(plotly)
library(ggplot2)

survival <- read.csv("data/Survival.csv", header=T, na.strings = "NA", stringsAsFactors = F, colClasses=
                       c(rep("character", times=2), "numeric", rep("factor", times=4), rep("numeric", times=4), "character", rep("numeric", times=5), rep("character", times=6)))
survival <- survival[c(-18:-23)]
survival$Date.stocked <- as.Date(survival$Date.stocked, format = "%m/%d/%y")
survival$Date.initial.count <- as.Date(survival$Date.initial.count, format = "%m/%d/%y")
survival$Date.imaged <- as.Date(survival$Date.imaged, format = "%m/%d/%y")
survival$TEMP <- survival$TRT
survival$TEMP <- gsub("\\<A\\>|\\<C\\>", "COLD", survival$TEMP)
survival$TEMP <- gsub("\\<B\\>|\\<D\\>", "WARM", survival$TEMP)
survival$FOOD <- survival$TRT
survival$FOOD <- gsub("\\<A\\>|\\<B\\>", "LOW", survival$FOOD)
survival$FOOD <- gsub("\\<C\\>|\\<D\\>", "HIGH", survival$FOOD)
survival$FOOD <- as.factor(survival$FOOD)
survival$TEMP <- as.factor(survival$TEMP)
survival$Dead.50.days <- (3*800)-survival$Live.50.days
survival$Dead.35.days <- 800-survival$Live.35.days
survival$TREAT <- as.factor(paste(survival$TEMP, "-", survival$FOOD))

100*mean(na.omit(survival$Live.50.days)/(800*3))
summary(100*(na.omit(survival$Live.50.days)/(800*3)))

plot(Live.50.days/(3*800) ~ Date.stocked, data=survival, col=TREAT, pch=8)

jpeg(file="results/boxplot-survival.jpeg", width = 900, height = 600)
plot(Live.50.days/(3*800) ~ TREAT, data=survival, col=c("skyblue3", "seagreen3", "indianred2",  "orange1"), main="Mean % survival across 12 groups per treatment", xlab="Treatment", ylab="Mean % survival", cex.lab=1.5, cex.main=1.5, par(mar=c(5,5,4.1,2.1)))
dev.off()

#levels(survival$TREAT) #color order ="skyblue3", "seagreen3", "indianred2",  "orange1"

levels(survival$TRT.REP) <- c("A1\nCOLD-LOW", "A2\nCOLD-LOW",     # A1 and A2
                              "B1\nWARM-LOW", "B2\nWARM-LOW",     # B1 and B2
                              "C1\nCOLD-HIGH", "C2\nCOLD-HIGH",   # C1 and C2
                              "D1\nWARM-HIGH", "D2\nWARM-HIGH")   # D1 and D2
#levels(survival$TRT.REP) # color order: "seagreen3",  "orange1", skyblue3", "indianred2"
jpeg(file="results/boxplot-survival-rep.jpeg", width = 900, height = 600)
plot(Live.50.days/(3*800) ~ TRT.REP, data=survival, col=c("seagreen3", "seagreen3", "orange1","orange1", "skyblue3", "skyblue3", "indianred2",  "indianred2"), main="Mean % survival between treatment rep\n6 groups per rep", xlab="Treatment", ylab="Mean % survival", cex.lab=1.6, cex.main=1.5, par(mar=c(5,5,4.1,2.1)))
dev.off()

jpeg(file="results/jitter-survival.jpeg", width = 700, height = 500)
ggplot(survival, aes(x=TREAT, y=100*Live.50.days/(3*800))) + geom_jitter(width=0.35, size=6, aes(color=TREAT)) + labs(title="Mean % survival, by treatment\n12 groups per treatment", y=("Mean % survival"), x="Treatment replicate") + scale_color_manual(values=c("skyblue3", "seagreen3", "indianred2",  "orange1")) + theme(text = element_text(size = 22)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
dev.off()


Survival.family.35 <- do.call(data.frame, aggregate(Live.35.days ~ Family+TEMP+FOOD+TRT.REP+TREAT, data = survival, FUN = function(x) c(mean = mean(x)/(800*3)*100, sd = sd(x)/(800*3)*100)))
names(Survival.family.35) <- c("Family", "TEMP", "FOOD", "TRT.REP", "TREAT", "Mean.Live.35", "SD.Live.35")
write.csv(file="data/Live.35.family.csv", x=Survival.family.35)

write.csv(file="data/Live.50.family.csv", x=(subset(survival, Live.50.days >0)))
 
Survival.family.50 <- subset(survival, Live.50.days >0)[c("Family", "TRT", "TRT.REP", "Live.50.days", "TEMP", "FOOD")]
Survival.family.50$Live.50.percent <- Survival.family.50$Live.50.days/(800*3)
write.csv(file="data/Live.50.family.csv", x=Survival.family.50)

library(reshape2)
dcast(Survival.family.50, month + day ~ variable)

# Model survival count data from day 35 (all silos separate)
glm.Date <- glm(cbind(Live.35.days, Dead.35.days) ~ Date.stocked, data=survival, quasibinomial)
summary(glm.Date)
glm.TEMP <- glm(cbind(Live.35.days, Dead.35.days) ~ TEMP, data=survival, quasibinomial)
summary(glm.TEMP)
glm.FOOD <- glm(cbind(Live.35.days, Dead.35.days) ~ FOOD, data=survival, quasibinomial)
summary(glm.FOOD)
glm.TREAT <- glm(cbind(Live.35.days, Dead.35.days) ~ FOOD*TEMP, data=survival, quasibinomial)
summary(glm.TREAT)

# Model survival count data from day 50 (silos combined by family) 
glm.Date.50 <- glm(cbind(Live.50.days, Dead.50.days) ~ Date.stocked, data=survival, quasibinomial)
summary(glm.Date.50)
glm.TEMP.50 <- glm(cbind(Live.50.days, Dead.50.days) ~ TEMP, data=survival, quasibinomial)
summary(glm.TEMP.50)
glm.FOOD.50 <- glm(cbind(Live.50.days, Dead.50.days) ~ FOOD, data=survival, quasibinomial)
summary(glm.FOOD.50)
glm.TREAT.50 <- glm(cbind(Live.50.days, Dead.50.days) ~ FOOD*TEMP, data=survival, quasibinomial)
summary(glm.TREAT.50)

library(ggplot2)
# Plot - color by temp, symbol by food 
jpeg(file="results/survival-time-treat.jpeg", width = 700, height = 500)
ggplot(survival, aes(x=Date.stocked, y=100*(Live.35.days/800))) + geom_point(size=5, aes(color=TEMP, shape=FOOD)) + labs(title="% Survival by treatment and date larvae was released", y=("Percent Survival"), x=("Date Released")) + ylim(0,60) + scale_color_manual(values=c("royalblue2", "tomato1"))+ scale_shape_manual(values=c(16, 8)) + theme(text = element_text(size = 20))
dev.off()

#Plot - color by treatment 
jpeg(file="results/survival-time-treat-col.jpeg", width = 700, height = 500)
ggplot(survival, aes(x=Date.stocked, y=100*(Live.35.days/800))) + geom_point(size=5, aes(color=TREAT)) + labs(title="% Survival by treatment and date larvae was released", y=("Percent Survival"), x=("Date Released")) + ylim(0,60) + scale_color_manual(values=c("skyblue3", "seagreen3", "indianred2",  "orange1"))+ scale_shape_manual(values=c(16, 8)) + theme(text = element_text(size = 20))
dev.off()

library(plotly)
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
# B2-3 only outlier - low percent live at release, and low survival (but not the lowest)

