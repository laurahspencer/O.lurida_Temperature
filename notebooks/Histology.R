# After redoing separate male/female stages, re-do analysis. 
histo <- read.csv("data/histology.csv", header=T, stringsAsFactors = T, na.strings = c("NA", " NA ", "TBD"))

#Convert a few columns to factors & reorder

histo$Week <- as.factor(histo$Week)
histo$Bag.. <- as.factor(histo$Bag..)
histo$FEMSTAGE.COL <- as.factor(histo$FEMSTAGE.COL)
histo$MALSTAGE.COL <- as.factor(histo$MALSTAGE.COL)
histo$SEX <- factor(histo$SEX, levels=c("I", "M", "HPM", "H", "HPF", "F"))
histo$TREAT <- factor(histo$TREAT, levels=c("PRE", "A" , "B", "C", "D", "Wild"))
histo$Date <- factor(histo$Date, levels=c("11/30/17", "12/20/17", "1/23/18", "1/4/18", "2/27/18", "2/9/18", "3/13/18", "3/23/18"))

summary(histo$SEX)
(56+194)/(8+25+45+57+56+194)

# % female 
table(subset(histo, Week==0)$SEX.COL)
25/(25+7+2+6) #62.5%, week 0 

table(subset(histo, Week==0)$FEMSTAGE.COL)
(13+9)/(4+12+13+9+2)

table(histo$SEX.COL) 
250/(250+57+8+70) # 64.9%, all weeks all treatments (including wild)

table(subset(histo, TREAT!="Wild")$SEX.COL) 
203/(203+45+8+59) # 64.4%, all weeks all treatments (not including wild)

table(subset(histo, TREAT=="Wild")$SEX.COL)
47/(47+12+11) #67.1%, Wild only 

# Prepare contingency tables 

# SEX
print(CT.sex.all <- table(histo$TREAT, histo$SEX))
print(CT.sex.Week <- table(histo$TREAT, histo$SEX, histo$Week))
print(CT.sex.FL <- table(subset(histo, Week=="0" | Week=="13")$TREAT, subset(histo, Week=="0" | Week=="13")$SEX))

# MALE STAGE
print(CT.MALSTAGE.all <- table(histo$TREAT, histo$MALSTAGE.COL))
print(CT.MALSTAGE.Week <- table(histo$TREAT, histo$MALSTAGE.COL, histo$Week))
print(CT.MALSTAGE.FL <- table(subset(histo, Week=="0" | Week=="13")$TREAT, subset(histo, Week=="0" | Week=="13")$MALSTAGE.COL))

# FEMALE STAGE
print(CT.fem.stage.all <- table(histo$TREAT, histo$FEMSTAGE.COL))
print(CT.fem.stage.Week <- table(histo$TREAT, histo$FEMSTAGE.COL, histo$Week))
print(CT.fem.stage.FL <- table(subset(histo, Week=="0" | Week=="13")$TREAT, subset(histo, Week=="0" | Week=="13")$FEMSTAGE.COL))

# CT for WILD only 
CT.W.sex <- table(subset(histo, TREAT=="Wild" | TREAT=="PRE")$Week, subset(histo, TREAT=="Wild" | TREAT=="PRE")$SEX)
CT.W.MALSTAGE <- table(subset(histo, TREAT=="Wild" | TREAT=="PRE")$Week, subset(histo, TREAT=="Wild" | TREAT=="PRE")$MALSTAGE.COL)
CT.W.fem.stage <- table(subset(histo, TREAT=="Wild" | TREAT=="PRE")$Week, subset(histo, TREAT=="Wild" | TREAT=="PRE")$FEMSTAGE.COL)

# CT for Treatment A only 
CT.A.sex <- table(subset(histo, TREAT=="A" | TREAT=="PRE")$Week, subset(histo, TREAT=="A" | TREAT=="PRE")$SEX)
CT.A.MALSTAGE <- table(subset(histo, TREAT=="A" | TREAT=="PRE")$Week, subset(histo, TREAT=="A" | TREAT=="PRE")$MALSTAGE.COL)
CT.A.fem.stage <- table(subset(histo, TREAT=="A" | TREAT=="PRE")$Week, subset(histo, TREAT=="A" | TREAT=="PRE")$FEMSTAGE.COL)

# CT for Treatment B only 
CT.B.sex <- table(subset(histo, TREAT=="B" | TREAT=="PRE")$Week, subset(histo, TREAT=="B" | TREAT=="PRE")$SEX)
CT.B.MALSTAGE <- table(subset(histo, TREAT=="B" | TREAT=="PRE")$Week, subset(histo, TREAT=="B" | TREAT=="PRE")$MALSTAGE.COL)
CT.B.fem.stage <- table(subset(histo, TREAT=="B" | TREAT=="PRE")$Week, subset(histo, TREAT=="B" | TREAT=="PRE")$FEMSTAGE.COL)

# CT for Treatment C only 
CT.C.sex <- table(subset(histo, TREAT=="C" | TREAT=="PRE")$Week, subset(histo, TREAT=="C" | TREAT=="PRE")$SEX)
CT.C.MALSTAGE <- table(subset(histo, TREAT=="C" | TREAT=="PRE")$Week, subset(histo, TREAT=="C" | TREAT=="PRE")$MALSTAGE.COL)
CT.C.fem.stage <- table(subset(histo, TREAT=="C" | TREAT=="PRE")$Week, subset(histo, TREAT=="C" | TREAT=="PRE")$FEMSTAGE.COL)

# CT for Treatment D only 
CT.D.sex <- table(subset(histo, TREAT=="D" | TREAT=="PRE")$Week, subset(histo, TREAT=="D" | TREAT=="PRE")$SEX)
CT.D.MALSTAGE <- table(subset(histo, TREAT=="D" | TREAT=="PRE")$Week, subset(histo, TREAT=="D" | TREAT=="PRE")$MALSTAGE.COL)
CT.D.fem.stage <- table(subset(histo, TREAT=="D" | TREAT=="PRE")$Week, subset(histo, TREAT=="D" | TREAT=="PRE")$FEMSTAGE.COL)

# Significance testing between treatments 

# All treatment weeks prior to reproductive conditioning combined (3-13):

print(sex.all <- table(subset(histo, Week!="0" & Week!="15" & Week!="16")$TREAT, subset(histo, Week!="0" & Week!="15" & Week!="16")$SEX))
chisq.test(sex.all[-1,], simulate.p.value = T, B = 10000)  #no sign. diff

# Collapsed sex designations (HPF = F, etc.)
print(sex.all.col <- table(subset(histo, Week!="0" & Week!="15" & Week!="16")$TREAT, subset(histo, Week!="0" & Week!="15" & Week!="16")$SEX.COL))
chisq.test(sex.all.col[-1,-1], simulate.p.value = T, B = 10000)  #no sign. diff

# Male gonad stage 
print(malstage.all <- table(subset(histo, Week!="0" & Week!="15" & Week!="16")$TREAT, subset(histo, Week!="0" & Week!="15" & Week!="16")$MALSTAGE.COL))
chisq.test(malstage.all[-1,], simulate.p.value = T, B = 10000)  #different, p=0.0118. How?

# Female gonad stage 
print(femalstage.all <- table(subset(histo, Week!="0" & Week!="15" & Week!="16")$TREAT, subset(histo, Week!="0" & Week!="15" & Week!="16")$FEMSTAGE.COL))
chisq.test(femalstage.all[-1,], simulate.p.value = T, B = 10000)  #no sign. diff


# Week 3 

print(sex.wk3 <- table(subset(histo, Week=="3")$TREAT, subset(histo, Week=="3")$SEX))
chisq.test(sex.wk3[-1,-1], simulate.p.value = T, B = 10000)  #no sign. diff

# Collapsed sex designations (HPF = F, etc.)
print(sex.wk3.col <- table(subset(histo, Week=="3")$TREAT, subset(histo, Week=="3")$SEX.COL))
fisher.test(sex.wk3.col[-1,-1], simulate.p.value = T, B = 10000)  #no sign. diff

# Male gonad stage 
print(malstage.wk3 <- table(subset(histo, Week=="3")$TREAT, subset(histo, Week=="3")$MALSTAGE.COL))
chisq.test(malstage.wk3[-1,], simulate.p.value = T, B = 10000)  #no sign. diff

# Female gonad stage 
print(femalstage.wk3 <- table(subset(histo, Week=="3")$TREAT, subset(histo, Week=="3")$FEMSTAGE.COL))
chisq.test(femalstage.wk3[-1,], simulate.p.value = T, B = 10000)  #no sign. diff


# Week 5 

print(sex.wk5 <- table(subset(histo, Week=="5")$TREAT, subset(histo, Week=="5")$SEX))
chisq.test(sex.wk5[-1,-1], simulate.p.value = T, B = 10000)  #no sign. diff

# Collapsed sex designations (HPF = F, etc.)
print(sex.wk5.col <- table(subset(histo, Week=="5")$TREAT, subset(histo, Week=="5")$SEX.COL))
fisher.test(sex.wk5.col[-1,-1], simulate.p.value = T, B = 10000)  #no sign. diff

# Male gonad stage 
print(malstage.wk5 <- table(subset(histo, Week=="5")$TREAT, subset(histo, Week=="5")$MALSTAGE.COL))
chisq.test(malstage.wk5[-1,], simulate.p.value = T, B = 10000)  #no sign. diff

# Female gonad stage 
print(femalstage.wk5 <- table(subset(histo, Week=="5")$TREAT, subset(histo, Week=="5")$FEMSTAGE.COL))
fisher.test(femalstage.wk5[-1,], simulate.p.value = T, B = 10000)  #no sign. diff


# Week 8 

print(sex.wk8 <- table(subset(histo, Week=="8")$TREAT, subset(histo, Week=="8")$SEX))
chisq.test(sex.wk8[-1,-1], simulate.p.value = T, B = 10000)  #no sign. diff

# Collapsed sex designations (HPF = F, etc.)
print(sex.wk8.col <- table(subset(histo, Week=="8")$TREAT, subset(histo, Week=="8")$SEX.COL))
fisher.test(sex.wk8.col[-1,-1], simulate.p.value = T, B = 10000)  #no sign. diff

# Male gonad stage 
print(malstage.wk8 <- table(subset(histo, Week=="8")$TREAT, subset(histo, Week=="8")$MALSTAGE.COL))
chisq.test(malstage.wk8[-1,], simulate.p.value = T, B = 10000)  #no sign. diff

# Female gonad stage 
print(femalstage.wk8 <- table(subset(histo, Week=="8")$TREAT, subset(histo, Week=="8")$FEMSTAGE.COL))
chisq.test(femalstage.wk8[-1,], simulate.p.value = T, B = 10000)  #no sign. diff


# Week 10 

print(sex.wk10 <- table(subset(histo, Week=="10")$TREAT, subset(histo, Week=="10")$SEX))
chisq.test(sex.wk10[-1,], simulate.p.value = T, B = 10000)  #no sign. diff

# Collapsed sex designations (HPF = F, etc.)
print(sex.wk10.col <- table(subset(histo, Week=="10")$TREAT, subset(histo, Week=="10")$SEX.COL))
chisq.test(sex.wk10.col[-1,-1], simulate.p.value = T, B = 10000)  #no sign. diff

# Male gonad stage 
print(malstage.wk10 <- table(subset(histo, Week=="10")$TREAT, subset(histo, Week=="10")$MALSTAGE.COL))
chisq.test(malstage.wk10[-1,], simulate.p.value = T, B = 10000)  #no sign. diff

# Female gonad stage 
print(femalstage.wk10 <- table(subset(histo, Week=="10")$TREAT, subset(histo, Week=="10")$FEMSTAGE.COL))
fisher.test(femalstage.wk10[-1,], simulate.p.value = T, B = 10000)  #no sign. diff



# Week 13 (end of treatment)

print(sex.wk13 <- table(subset(histo, Week=="13")$TREAT, subset(histo, Week=="13")$SEX))
chisq.test(sex.wk13[-1,], simulate.p.value = T, B = 10000)  #can't use chi-sq
fisher.test(sex.wk13[-1,], simulate.p.value = T, B = 10000)  #no sign. diff

# Collapsed sex designations (HPF = F, etc.)
print(sex.wk13.col <- table(subset(histo, Week=="13")$TREAT, subset(histo, Week=="13")$SEX.COL))
chisq.test(sex.wk13.col[-1,-1], simulate.p.value = T, B = 10000)  #no sign. diff

# Male gonad stage 
print(malstage.wk13 <- table(subset(histo, Week=="13")$TREAT, subset(histo, Week=="13")$MALSTAGE.COL))
chisq.test(malstage.wk13[-1,], simulate.p.value = T, B = 10000)  # diff! 
require(rcompanion)
pairwiseNominalIndependence(malstage.wk13[2:6,],fisher = TRUE,gtest  = FALSE, chisq  = FALSE, digits = 3)


# Female gonad stage 
print(femalstage.wk13 <- table(subset(histo, Week=="13")$TREAT, subset(histo, Week=="13")$FEMSTAGE.COL))
chisq.test(femalstage.wk13[-1,], simulate.p.value = T, B = 10000)  #no sign. diff


# Week 15 (first reproductive conditioning sampling)

print(sex.wk15 <- table(subset(histo, Week=="15")$TREAT, subset(histo, Week=="15")$SEX))
chisq.test(sex.wk15[2:5,], simulate.p.value = T, B = 10000)  #cant use

# Collapsed sex designations (HPF = F, etc.)
print(sex.wk15.col <- table(subset(histo, Week=="15")$TREAT, subset(histo, Week=="15")$SEX.COL))
chisq.test(sex.wk15.col[2:5,-1], simulate.p.value = T, B = 10000)  #no sign. diff

# Male gonad stage 
print(malstage.wk15 <- table(subset(histo, Week=="15")$TREAT, subset(histo, Week=="15")$MALSTAGE.COL))
chisq.test(malstage.wk15[2:5,], simulate.p.value = T, B = 10000)  #Yes diff. 

# Female gonad stage 
print(femalstage.wk15 <- table(subset(histo, Week=="15")$TREAT, subset(histo, Week=="15")$FEMSTAGE.COL))
chisq.test(femalstage.wk15[2:5,-5], simulate.p.value = T, B = 10000)  #no sign. diff


# Week 16 (last sampling, 3rd week of reproductive conditioning

print(sex.wk16 <- table(subset(histo, Week=="16")$TREAT, subset(histo, Week=="16")$SEX))
fisher.test(sex.wk16[2:5,], simulate.p.value = T, B = 10000)  #yes diff 
chisq.test(sex.wk16[2:5,], simulate.p.value = T, B = 10000)  #no sign. diff

# Collapsed sex designations (HPF = F, etc.)
print(sex.wk16.col <- table(subset(histo, Week=="16")$TREAT, subset(histo, Week=="16")$SEX.COL))
chisq.test(sex.wk16.col[2:5,-1], simulate.p.value = T, B = 10000)  #no sign. diff

# Male gonad stage 
print(malstage.wk16 <- table(subset(histo, Week=="16")$TREAT, subset(histo, Week=="16")$MALSTAGE.COL))
chisq.test(malstage.wk16[2:5,], simulate.p.value = T, B = 10000)  #Yes diff. 

# Female gonad stage 
print(femalstage.wk16 <- table(subset(histo, Week=="16")$TREAT, subset(histo, Week=="16")$FEMSTAGE.COL))
chisq.test(femalstage.wk16[2:5,], simulate.p.value = T, B = 10000)  #no sign. diff


# Weeks 15 & 16 (2 conditioning weeks) 

print(sex.wkcond <- table(subset(histo, Week=="15" | Week=="16")$TREAT, subset(histo, Week=="15" | Week=="16")$SEX))
fisher.test(sex.wkcond[2:5,], simulate.p.value = T, B = 10000)  #yes diff 
chisq.test(sex.wkcond[2:5,], simulate.p.value = T, B = 10000)  #no sign. diff

# Collapsed sex designations (HPF = F, etc.)
print(sex.wkcond.col <- table(subset(histo, Week=="15" | Week=="16")$TREAT, subset(histo, Week=="15" | Week=="16")$SEX.COL))
chisq.test(sex.wkcond.col[2:5,-1], simulate.p.value = T, B = 10000)  #no sign. diff

# Male gonad stage 
print(malstage.wkcond <- table(subset(histo, Week=="15" | Week=="16")$TREAT, subset(histo, Week=="15" | Week=="16")$MALSTAGE.COL))
chisq.test(malstage.wkcond[2:5,], simulate.p.value = T, B = 10000)  #Yes diff
pairwiseNominalIndependence(malstage.wkcond[2:5,],fisher = TRUE,gtest  = FALSE, chisq  = FALSE, digits = 3)


# Female gonad stage 
print(femalstage.wkcond <- table(subset(histo, Week=="15" | Week=="16")$TREAT, subset(histo, Week=="15" | Week=="16")$FEMSTAGE.COL))
chisq.test(femalstage.wkcond[2:5,], simulate.p.value = T, B = 10000)  #no sign. diff


# Compare Weeks 0 to 13

print(sex.wkcond <- table(histo$Week, histo$SEX))
fisher.test(sex.wkcond[c(1,6),], simulate.p.value = T, B = 10000)  #yes diff 
chisq.test(sex.wkcond[c(1,6),], simulate.p.value = T, B = 10000)  #no sign. diff

# Collapsed sex designations (HPF = F, etc.)
print(sex.wkcond.col <- table(histo$Week, histo$SEX.COL))
chisq.test(sex.wkcond.col[c(1,6),-1], simulate.p.value = T, B = 10000)  #no sign. diff

# Male gonad stage 
print(malstage.wkcond <- table(histo$Week, histo$MALSTAGE.COL))
chisq.test(malstage.wkcond[c(1,6),-1], simulate.p.value = T, B = 10000)  #Yes diff. 

# Female gonad stage 
print(femalstage.wkcond <- table(histo$Week, histo$FEMSTAGE.COL))
chisq.test(femalstage.wkcond[c(1,6),-1], simulate.p.value = T, B = 10000)  #no sign. diff

# Female gonad stage, weeks 13 vs. 16
print(femalstage.wkcond <- table(histo$Week, histo$FEMSTAGE.COL))
chisq.test(femalstage.wkcond[c(6,8),-1], simulate.p.value = T, B = 10000)  #no sign. diff


# Assess ripe oocyte size 

oocyte <- read.csv("data/2019-04-02_Oocyte-length.csv", header=T, stringsAsFactors = T)
oocyte.meanlength <- aggregate(Length ~ Sample+TEMP+FOOD+Week+TREAT+TREAT.NAME, oocyte, mean)

oocyte.size <- oocyte %>%
  group_by(Sample, Week, TEMP, FOOD, TREAT.NAME) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))

oocyte.size$cv <- oocyte.size$sd_length/oocyte.size$mean_length
summary(oocyte.size$cv) #average variation within an individual 

summary(oocyte$Length)
sd(oocyte$Length)
aggregate(Length ~ TEMP+Week, oocyte, mean)
aggregate(Length ~ TEMP, oocyte, sd)
aggregate(Sample ~ TEMP+FOOD, oocyte, FUN = length)
length(unique(subset(oocyte, TREAT=="D")$Sample)) # number of samples for each treatment

# Compare oocyte length by treatment, using mean size per individual 
summary(aov(Length ~ TEMP*FOOD, data=subset(oocyte.meanlength, TEMP!="Wild")))
TukeyHSD(aov(Length ~ TEMP*FOOD, data=subset(oocyte.meanlength, TEMP!="Wild")))
summary(aov(Length ~ TEMP, data=subset(oocyte.meanlength, TEMP!="Wild"))) 
summary(aov(Length ~ FOOD, data=subset(oocyte.meanlength, TEMP!="Wild"))) 

summary(aov(Length ~ TEMP, data=subset(oocyte.meanlength, TEMP!="Wild" & Week!=13))) 
summary(aov(Length ~ TEMP, data=subset(oocyte.meanlength, TEMP!="Wild" & Week!=15))) 
summary(aov(Length ~ TEMP, data=subset(oocyte.meanlength, TEMP!="Wild" & Week!=16))) 

# test whether i can just use 12 oocytes 
test <- oocyte %>% group_by(Sample, Week, TEMP, FOOD, TREAT.NAME) %>% sample_n(size = 12)
oocyte.meanlength.test <- aggregate(Length ~ Sample+TEMP+FOOD+Week+TREAT+TREAT.NAME, test, mean)

# Compare oocyte length by treatment, using mean size per individual 
summary(aov(Length ~ TEMP*FOOD, data=subset(oocyte.meanlength.test, TEMP!="Wild")))

#Plot mean oocyte length for each oyster sample, by treatment 
pdf(file="results/Stage3-oocyte-size.pdf", width=5, height=6)
ggplot(data=subset(oocyte.meanlength, TREAT!="Wild"), aes(x=as.factor(Week), y=Length, fill=TEMP)) + geom_boxplot() + theme_bw(base_size = 13) + ggtitle(label = "Stage 3 oocyte size\nweeks 13-16") + ylab("Maximum oocyte length (um)") + scale_fill_manual(values=c("slategray", "slategray1", "lightcoral",  "rosybrown1"), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food")) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
dev.off()

pdf(file="results/Stage3-oocyte-size.pdf", width=5.5, height=6)
ggplot(data=subset(oocyte.meanlength, TREAT!="Wild"), aes(x=TREAT.NAME, y=Length, col=TREAT.NAME)) + geom_boxplot(lwd=0.5) + theme_bw(base_size = 11) + ggtitle(label = "Ripe oocyte size") + ylab("Maximum oocyte length (um)") + scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("7°C+high-food", "7°C+low-food", "10°C+high-food", "10°C+low-food")) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) + geom_jitter(width = .2,aes(shape=as.factor(Week)),size=2) + scale_shape_manual(values=c(15,16,17), name="Date Oocyte Sampled", labels=c("February 27", "March 13", "March 23"))
dev.off()
#c(1,2,0)

ggplot(aes(x=TEMP:FOOD, y=Length)) + 
  geom_boxplot(lwd=0.5, aes(col=TEMP:FOOD), outlier.shape = NA)

# For PCSGA presentation 
subset(oocyte.meanlength, TREAT!="Wild") %>%
mutate(FOOD = factor(FOOD, levels=c("LOW", "HIGH"))) %>% 
ggplot(aes(x=TEMP:FOOD, y=Length)) + 
  geom_boxplot(lwd=0.5, aes(col=TEMP:FOOD), outlier.shape = NA) + 
  theme_bw(base_size = 11) + 
  ggtitle(label = "Ripe oocyte size") + 
  ylab("Maximum oocyte length (um)") + 
  scale_color_manual(values=c('darkslateblue','#0571b0','chocolate2','#ca0020'), name=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  geom_jitter(width = .2,size=2, aes(col=TEMP:FOOD)) + 
  ylim(c(60,110)) + 
  theme(legend.position = "none") #+ scale_shape_manual(values=c(15,16,17), name="Date Oocyte Sampled", labels=c("February 27", "March 13", "March 23"))

# with wild too (but only 1 wild sample)
ggplot(data=oocyte, aes(x=TREAT.NAME, y=Length, fill=TREAT.NAME)) +
  geom_boxplot() + ggtitle(label = "Stage 3 oocyte size, weeks 13-16") + ylab("Maximum oocyte length (um)") + scale_fill_manual(values=c("skyblue3", "seagreen3", "indianred2",  "orange1", "gray30")) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

# Plot % of each sex over time 
test <- table(subset(histo, TREAT!="Wild")$Week, subset(histo, TREAT!="Wild")$SEX.COL)[,-1]
test <- as.data.frame(prop.table(test, margin = 1))
ggplot(data=test, aes(x=Var1, y=Freq, group=Var2, col=Var2)) + geom_line() + theme_bw(base_size = 12) + geom_point(size=3) + ylab("Frequency")
  
# Plot % of each sex over time 
test <- table(subset(histo, TREAT!="Wild")$Week, subset(histo, TREAT!="Wild")$SEX.COL, subset(histo, TREAT!="Wild")$TREAT)[-1,-1,2:5]
test <- as.data.frame(prop.table(test, margin = 1))
ggplot(data=subset(test, Var3=="A"), aes(x=Var1, y=Freq, group=Var2, col=Var2)) + geom_line() + theme_bw(base_size = 12) + geom_point(size=3) + ylab("Frequency")
ggplot(data=subset(test, Var3=="C"), aes(x=Var1, y=Freq, group=Var2, col=Var2)) + geom_line() + theme_bw(base_size = 12) + geom_point(size=3) + ylab("Frequency")
ggplot(data=subset(test, Var3=="B"), aes(x=Var1, y=Freq, group=Var2, col=Var2)) + geom_line() + theme_bw(base_size = 12) + geom_point(size=3) + ylab("Frequency")
ggplot(data=subset(test, Var3=="D"), aes(x=Var1, y=Freq, group=Var2, col=Var2)) + geom_line() + theme_bw(base_size = 12) + geom_point(size=3) + ylab("Frequency")


# Plot male stage over time 
test <- table(subset(histo, TREAT!="Wild")$Week, subset(histo, TREAT!="Wild")$MALSTAGE.COL)
test <- as.data.frame(prop.table(test, margin = 1))
ggplot(data=subset(test, Var2!=0), aes(x=Var1, y=Freq, group=Var2, col=Var2)) + geom_line() + theme_bw(base_size = 12) + geom_point(size=3) + ylab("Frequency") + scale_color_manual(values=c("#6baed6","#3182bd","#08519c","#bdd7e7"))

# Plot female stage over time 
test <- table(subset(histo, TREAT!="Wild")$Week, subset(histo, TREAT!="Wild")$FEMSTAGE.COL)
test <- as.data.frame(prop.table(test, margin = 1))
ggplot(data=subset(test, Var2!=0), aes(x=Var1, y=Freq, group=Var2, col=Var2)) + geom_line() + theme_bw(base_size = 12) + geom_point(size=3) + ylab("Frequency") + scale_color_manual(values=c("#df65b0","#dd1c77","#980043","#d7b5d8")) 


### =================================================================== 
### ------------------------------------------ Line plots for paper!!! -------------

# Plot # female, male over time, by treatment (exclude H)

treats <- as.character(unique(subset(histo, TREAT!="PRE")$TREAT))

sexcol.treats <- list()
for (i in 1:length(treats)) {
  sexcol.treats[[i]] <- as.data.frame(table(subset(histo, TREAT==treats[i] | TREAT=="PRE" )$Week, subset(histo, TREAT==treats[i] | TREAT=="PRE" )$SEX.COL), margin=1)
}

# Important for plotting - set the wild samples from weeks 15 & 16 to NA
sexcol.treats[[1]]$Freq[which(sexcol.treats[[1]]$Var1 == c(15,16))] <- NA

p <- list()
# plot wild samples, w/o vertical lines 
p[[1]] <- ggplot(data=subset(sexcol.treats[[1]], Var1!=0  & (Var2=="F" | Var2=="M" | Var2=="H")), aes(x=Var1, y=Freq, group=Var2, col=Var2)) + 
  geom_line() + theme_bw(base_size = 12) + geom_point(size=2.5) + 
  theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title = element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylim(c(0,10)) + 
  scale_color_manual(values=c("gray10", "gray45", "gray75")) #+ 
#scale_color_manual(values=c("#b3cde3", "#8c96c6", "#88419d"))

for (i in 2:length(treats)){
  p[[i]] <- ggplot(data=subset(sexcol.treats[[i]], Var1!=0 & (Var2=="F" | Var2=="M" | Var2=="H")), aes(x=Var1, y=Freq, group=Var2, col=Var2)) + 
    geom_line() + theme_bw(base_size = 12) + geom_point(size=2.5) + 
    theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title = element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylim(c(0,10)) + 
    geom_vline(xintercept = 5.1, linetype="dashed", color = "gray50", size=.5) + 
    scale_color_manual(values=c("gray10", "gray45", "gray75")) #+ 
  #scale_color_manual(values=c("#b3cde3", "#8c96c6", "#88419d"))
}
pdf(file = "results/sex-line-plots.pdf", width = 4.5, height = 8)
do.call(grid.arrange, c(c(p[1], p[2], p[3], p[4], p[5]), list(ncol=1, nrow=5)))
dev.off()

pdf(file = "results/sex-line-plots-legend.pdf", width = 4.5, height = 2)
ggplot(data=subset(sexcol.treats[[i]], Var1!="0" & (Var2=="F" | Var2=="M" | Var2=="H")), aes(x=Var1, y=Freq, group=Var2, col=Var2)) + geom_line() + theme_bw(base_size = 12) + geom_point(size=2.5) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "bottom",
        legend.title=element_text(size=13), legend.text = element_text(size=13)) + 
  guides(color = guide_legend(override.aes = list(size=6))) + scale_x_discrete(labels= c("Dec 20", "Jan 4", "Jan 23", "Feb 9", "Feb 27", "Mar 13", "Mar 23")) + 
  scale_color_manual(values=c("gray10", "gray45", "gray75"), name=element_blank(), labels = c("Female / HPF", "Male / HPM", "Hermaph.")) #+
#scale_color_manual(values=c("#b3cde3", "#8c96c6", "#88419d"), name=element_blank(), labels = c("Female / HPF", "Male / HPM", "Hermaph."))
dev.off()

# Plot male stage over time 

# For both male and female gamete stage, a column that combines stages 2 and 3, and replaces stage 0 with NA (none present)
histo <- histo %>% 
  mutate(MALSTAGE.COL2 = droplevels(replace(MALSTAGE.COL, MALSTAGE.COL=="2", "3")), 
         FEMSTAGE.COL2 = droplevels(replace(FEMSTAGE.COL, FEMSTAGE.COL=="2", "3"))) #%>%
# mutate(MALSTAGE.COL2 = droplevels(replace(MALSTAGE.COL2, MALSTAGE.COL=="0", "NA")), 
#        FEMSTAGE.COL2 = droplevels(replace(FEMSTAGE.COL2, FEMSTAGE.COL=="0", "NA")))

# Calculate # of samples in each category, and plot 

malestage.treats <- list()
for (i in 1:length(treats)) {
  malestage.treats[[i]] <- as.data.frame(table(subset(histo, c(TREAT==treats[i] | TREAT=="PRE"))$Week, subset(histo, c(TREAT==treats[i] | TREAT=="PRE"))$MALSTAGE.COL2), margin=1)
  malestage.treats[[i]]["Var2"] <- droplevels(malestage.treats[[i]]["Var2"])
}

# Important for plotting - set the wild samples from weeks 15 & 16 to NA
malestage.treats[[1]]$Freq[which(malestage.treats[[1]]$Var1 == c(15,16))] <- NA

p <- list()
p[[1]] <- ggplot(data=subset(malestage.treats[[1]], Var1!="0" & Var2!="0"), aes(x=Var1, y=Freq, group=Var2, col=Var2)) + 
  scale_y_continuous(position = "right", limits = c(0,10), breaks=c(0,5,10)) + 
  geom_line() + theme_bw(base_size = 12) + geom_point(size=2.5) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none", plot.title = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  scale_color_manual(values=c("gray45", "gray10", "gray75")) 
#scale_color_manual(values=c("#bdc9e1", "#67a9cf", "#02818a"))

for (i in 2:length(treats)){
  p[[i]] <- ggplot(data=subset(malestage.treats[[i]], Var1!="0" & Var2!="0"), aes(x=Var1, y=Freq, group=Var2, col=Var2)) + 
    scale_y_continuous(position = "right", limits = c(0,10), breaks=c(0,5,10)) + 
    geom_line() + theme_bw(base_size = 12) + geom_point(size=2.5) + 
    geom_vline(xintercept = 5.1, linetype="dashed", color = "gray50", size=.5) + 
    theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none", plot.title = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
    scale_color_manual(values=c("gray45", "gray10", "gray75"))
  #scale_color_manual(values=c("#bdc9e1", "#67a9cf", "#02818a"))
}
pdf(file = "results/malestage-line-plots.pdf", width = 4.8, height = 8)
do.call(grid.arrange, c(c(p[1], p[2], p[3], p[4], p[5]), list(ncol=1, nrow=5)))
dev.off()

# Generate plot just for legend 
pdf(file = "results/malestage-line-plots-legend.pdf", width = 4.5, height = 2)
ggplot(data=subset(malestage.treats[[i]], Var1!="0" & Var2!="0"), aes(x=Var1, y=Freq, group=Var2, col=Var2)) + geom_line() + theme_bw(base_size = 12) + geom_point(size=2.5) + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "bottom", axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.title=element_text(size=14), legend.text = element_text(size=14))  + guides(color = guide_legend(override.aes = list(size=6))) + scale_x_discrete(labels= c("Dec 20", "Jan 4", "Jan 23", "Feb 9", "Feb 27", "Mar 13", "Mar 23")) + 
  scale_color_manual(values=c("gray45", "gray10", "gray75"),name=element_blank(), labels = c("Early", "Advanced/Ripe", "Spawned")) #+
#scale_color_manual(values=c("#bdc9e1", "#67a9cf", "#02818a"),name=element_blank(), labels = c("Early", "Advanced/Ripe", "Spawned"))
dev.off()

# plot female stage over time 

femalestage.treats <- list()
for (i in 1:length(treats)) {
  femalestage.treats[[i]] <- as.data.frame(table(subset(histo, TREAT==treats[i] | TREAT=="PRE" )$Week, subset(histo, TREAT==treats[i] | TREAT=="PRE" )$FEMSTAGE.COL2), margin=1)
}

# Important for plotting - set the wild samples from weeks 15 & 16 to NA
femalestage.treats[[1]]$Freq[which(femalestage.treats[[1]]$Var1 == c(15,16))] <- NA

p <- list()

p[[1]] <- ggplot(data=subset(femalestage.treats[[1]], Var1!="0" & Var2!="0"), aes(x=Var1, y=Freq, group=Var2, col=Var2)) + geom_line() + theme_bw(base_size = 12) + geom_point(size=2) + ylim(c(0,10)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none", plot.title = element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  scale_color_manual(values=c("gray45", "gray10", "gray75")) 
#scale_color_manual(values=c("#fdbe85", "#fd8d3c", "#d94701")) 

for (i in 2:length(treats)){
  p[[i]] <- ggplot(data=subset(femalestage.treats[[i]], Var1!="0" & Var2!="0"), aes(x=Var1, y=Freq, group=Var2, col=Var2)) + geom_line() + theme_bw(base_size = 12) + geom_point(size=2) + ylim(c(0,10)) + 
    geom_vline(xintercept = 5.1, linetype="dashed", color = "gray50", size=.5) + 
    theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none", plot.title = element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
    scale_color_manual(values=c("gray45", "gray10", "gray75")) 
  #scale_color_manual(values=c("#fdbe85", "#fd8d3c", "#d94701")) 
}
pdf(file = "results/femalestage-line-plots.pdf", width = 4.5, height = 8)
do.call(grid.arrange, c(c(p[1], p[2], p[3], p[4], p[5]), list(ncol=1, nrow=5)))
dev.off()

pdf(file = "results/femalestage-line-plots-legend.pdf", width = 4.5, height = 2)
ggplot(data=subset(femalestage.treats[[i]], Var1!="0" & Var2!="0"), aes(x=Var1, y=Freq, group=Var2, col=Var2)) + geom_line() + theme_bw(base_size = 12) + geom_point(size=2) + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "bottom", axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.title=element_text(size=14), legend.text = element_text(size=14)) + guides(color = guide_legend(override.aes = list(size=6))) + scale_x_discrete(labels= c("Dec 20", "Jan 4", "Jan 23", "Feb 9", "Feb 27", "Mar 13", "Mar 23")) + 
  scale_color_manual(values=c("gray45", "gray10", "gray75"),name=element_blank(), labels = c("Early", "Advanced/Ripe", "Spawned")) #+
#scale_color_manual(values=c("#fdbe85", "#fd8d3c", "#d94701"),name=element_blank(), labels = c("Early", "Advanced/Ripe", "Spawned"))
dev.off()



#### extra plots for fun 
 
library(ggplot2)
library(dplyr)

ggplot(data=subset(histo, MALE.STAGE !="NA"), aes(y=FEMALE.STAGE, x=MALE.STAGE, col=SEX.COL)) + theme_bw(base_size = 12) + geom_jitter(size=2) + xlab("sperm stage") + ylab("egg stage")  

ggplot(data=subset(histo, MALE.STAGE !="NA"), aes(x=FEMALE.STAGE, y=Est..Tissue.Weight..g., col=SEX.COL)) + theme_bw(base_size = 12) + geom_jitter(size=2) 

ggplot(data=subset(histo, MALE.STAGE !="NA"), aes(x=MALE.STAGE, y=Est..Tissue.Weight..g., col=SEX.COL)) + theme_bw(base_size = 12) + geom_jitter(size=2) 


# 
# # Histograms of female stage 
# p <- list()
# for (i in 1:length(treats)){
#   p[[i]] <- ggplot(data=subset(histo, TREAT==treats[i] & FEMSTAGE.COL2!="NA"), aes(FEMSTAGE.COL2, fill=FEMSTAGE.COL2)) + geom_histogram( stat="count", binwidth = 1) + theme_bw(base_size = 12) + 
#     theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none", plot.title = element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + scale_fill_manual(values=c("gray90", "gray45", "gray10", "gray75"), drop=FALSE) + scale_x_discrete(drop=FALSE)
# }
# pdf(file = "results/femalestage-hist.pdf", width = 2, height = 8)
# do.call(grid.arrange, c(c(p[1], p[4], p[2], p[5], p[3]), list(ncol=1, nrow=5)))
# dev.off()


# Plot sex and male/female stage over time for each treatment


# WILD - MUD BAY BEACH SAMPLINGS 

jpeg(file = "results/gonad-barplots-wild.jpeg", width = 850, height = 550)

par(mfrow = c(1, 3), mar=c(1, 2, 3, 0), oma=c(4,4,4,0), col="gray30")

print(barplot(t(prop.table(CT.W.sex[1:6,], 1)), main="Dominant gonad sex", xlab="Week Sampled", ylab="% Sampled", las=1, col=c("gray75", "#08519c", "#3182bd", "purple3","mediumorchid3", "#df65b0"),  cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, args.legend = list(x = "topright", bty = "n", inset=c(-1, 0), title="Gonad Sex", cex=1.5)))
chisq.test(CT.W.sex[1:6,], simulate.p.value = T, B = 10000) 

colnames(CT.W.MALSTAGE) <- c("None present (0)", "Early (1)", "Advanced (2)", "Ripe (3)", "Spawned (4)")
print(barplot(t(prop.table(CT.W.MALSTAGE[1:6,], 1)), main="Male", xlab=NA, ylab=NA, las=1, col=c("#eff3ff","#6baed6","#3182bd","#08519c","#bdd7e7"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-.5, 0), title="Gonad Stage", cex=1.5)))
chisq.test(CT.W.MALSTAGE[1:6,], simulate.p.value = T, B = 10000) 

colnames(CT.D.fem.stage) <- c("None present (0)", "Early (1)", "Advanced (2)", "Ripe (3)", "Spawned (4)")
print(barplot(t(prop.table(CT.D.fem.stage[1:6,], 1)), main="Female", xlab=NA, ylab=NA, las=1, col=c("#f1eef6","#df65b0","#dd1c77","#980043","#d7b5d8"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-.5, 0), title="Gonad Stage", cex=1.5)))
chisq.test(CT.W.fem.stage[1:6,], simulate.p.value = T, B = 10000) 

mtext("Mud Bay", side=3, line = .5, outer = T, cex = 2, col = "gray30")
mtext("% Sampled", side=2, line = 1, outer = T, cex = 1.5, col = "gray30")
mtext("Week Sampled", side=1, line = 2.5, outer = T, cex = 1.5, col = "gray30")

dev.off()

# A = 7C low nutrition

jpeg(file = "results/gonad-barplots-7C-low.jpeg", width = 1015, height = 550)

par(mfrow = c(1, 3), mar=c(1, 2, 3, 0), oma=c(4,4,4,0), col="gray30")

print(barplot(t(prop.table(CT.A.sex, 1)), main="Dominant gonad sex", xlab="Week Sampled", ylab="% Sampled", las=1, col=c("gray75", "#08519c", "#3182bd", "purple3","mediumorchid3", "#df65b0"),  cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, args.legend = list(x = "topright", bty = "n", inset=c(-1, 0), title="Gonad Sex", cex=1.5)))
chisq.test(CT.A.sex, simulate.p.value = T, B = 10000) 

colnames(CT.A.MALSTAGE) <- c("None present (0)", "Early (1)", "Advanced (2)", "Ripe (3)", "Spawned (4)")
print(barplot(t(prop.table(CT.A.MALSTAGE, 1)), main="Male", xlab=NA, ylab=NA, las=1, col=c("#eff3ff","#6baed6","#3182bd","#08519c","#bdd7e7"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-.5, 0), title="Gonad Stage", cex=1.5)))
chisq.test(CT.A.MALSTAGE, simulate.p.value = T, B = 10000) 

colnames(CT.A.fem.stage) <- c("None present (0)", "Early (1)", "Advanced (2)", "Ripe (3)", "Spawned (4)")
print(barplot(t(prop.table(CT.A.fem.stage, 1)), main="Female", xlab=NA, ylab=NA, las=1, col=c("#f1eef6","#df65b0","#dd1c77","#980043","#d7b5d8"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-.5, 0), title="Gonad Stage", cex=1.5)))
chisq.test(CT.A.fem.stage, simulate.p.value = T, B = 10000) 

mtext("7°C - Low Nutrition", side=3, line = .5, outer = T, cex = 2, col = "gray30")
mtext("% Sampled", side=2, line = 1, outer = T, cex = 1.5, col = "gray30")
mtext("Week Sampled", side=1, line = 2.5, outer = T, cex = 1.5, col = "gray30")

dev.off()

# C = 7C high nutrition

jpeg(file = "results/gonad-barplots-7C-high.jpeg", width = 1015, height = 550)

par(mfrow = c(1, 3), mar=c(1, 2, 3, 0), oma=c(4,4,4,0), col="gray30")

print(barplot(t(prop.table(CT.C.sex, 1)), main="Dominant gonad sex", xlab="Week Sampled", ylab="% Sampled", las=1, col=c("gray75", "#08519c", "#3182bd", "purple3","mediumorchid3", "#df65b0"),  cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, args.legend = list(x = "topright", bty = "n", inset=c(-1, 0), title="Gonad Sex", cex=1.5)))
chisq.test(CT.C.sex, simulate.p.value = T, B = 10000) 

colnames(CT.C.MALSTAGE) <- c("None present (0)", "Early (1)", "Advanced (2)", "Ripe (3)", "Spawned (4)")
print(barplot(t(prop.table(CT.C.MALSTAGE, 1)), main="Male", xlab=NA, ylab=NA, las=1, col=c("#eff3ff","#6baed6","#3182bd","#08519c","#bdd7e7"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-.5, 0), title="Gonad Stage", cex=1.5)))
chisq.test(CT.C.MALSTAGE, simulate.p.value = T, B = 10000) 

colnames(CT.C.fem.stage) <- c("None present (0)", "Early (1)", "Advanced (2)", "Ripe (3)", "Spawned (4)")
print(barplot(t(prop.table(CT.C.fem.stage, 1)), main="Female", xlab=NA, ylab=NA, las=1, col=c("#f1eef6","#df65b0","#dd1c77","#980043","#d7b5d8"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-.5, 0), title="Gonad Stage", cex=1.5)))
chisq.test(CT.C.fem.stage, simulate.p.value = T, B = 10000) 

mtext("7°C - High Nutrition", side=3, line = .5, outer = T, cex = 2, col = "gray30")
mtext("% Sampled", side=2, line = 1, outer = T, cex = 1.5, col = "gray30")
mtext("Week Sampled", side=1, line = 2.5, outer = T, cex = 1.5, col = "gray30")

dev.off()


# B = 10C low nutrition

jpeg(file = "results/gonad-barplots-10C-low.jpeg", width = 1015, height = 550)

par(mfrow = c(1, 3), mar=c(1, 2, 3, 0), oma=c(4,4,4,0), col="gray30")

print(barplot(t(prop.table(CT.B.sex, 1)), main="Dominant gonad sex", xlab="Week Sampled", ylab="% Sampled", las=1, col=c("gray75", "#08519c", "#3182bd", "purple3","mediumorchid3", "#df65b0"),  cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, args.legend = list(x = "topright", bty = "n", inset=c(-1, 0), title="Gonad Sex", cex=1.5)))
chisq.test(CT.B.sex, simulate.p.value = T, B = 10000) 

colnames(CT.B.MALSTAGE) <- c("None present (0)", "Early (1)", "Advanced (2)", "Ripe (3)", "Spawned (4)")
print(barplot(t(prop.table(CT.B.MALSTAGE, 1)), main="Male", xlab=NA, ylab=NA, las=1, col=c("#eff3ff","#6baed6","#3182bd","#08519c","#bdd7e7"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-.5, 0), title="Gonad Stage", cex=1.5)))
chisq.test(CT.B.MALSTAGE, simulate.p.value = T, B = 10000) 

colnames(CT.B.fem.stage) <- c("None present (0)", "Early (1)", "Advanced (2)", "Ripe (3)", "Spawned (4)")
print(barplot(t(prop.table(CT.B.fem.stage, 1)), main="Female", xlab=NA, ylab=NA, las=1, col=c("#f1eef6","#df65b0","#dd1c77","#980043","#d7b5d8"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-.5, 0), title="Gonad Stage", cex=1.5)))
chisq.test(CT.B.fem.stage, simulate.p.value = T, B = 10000) 

mtext("10°C - Low Nutrition", side=3, line = .5, outer = T, cex = 2, col = "gray30")
mtext("% Sampled", side=2, line = 1, outer = T, cex = 1.5, col = "gray30")
mtext("Week Sampled", side=1, line = 2.5, outer = T, cex = 1.5, col = "gray30")

dev.off()

# D = 10C high nutrition

jpeg(file = "results/gonad-barplots-10C-high.jpeg", width = 1015, height = 550)

par(mfrow = c(1, 3), mar=c(1, 2, 3, 0), oma=c(4,4,4,0), col="gray30")

print(barplot(t(prop.table(CT.D.sex, 1)), main="Dominant gonad sex", xlab="Week Sampled", ylab="% Sampled", las=1, col=c("gray75", "#08519c", "#3182bd", "purple3","mediumorchid3", "#df65b0"),  cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, args.legend = list(x = "topright", bty = "n", inset=c(-1, 0), title="Gonad Sex", cex=1.5)))
chisq.test(CT.D.sex, simulate.p.value = T, B = 10000)  #sign.

colnames(CT.D.MALSTAGE) <- c("None present (0)", "Early (1)", "Advanced (2)", "Ripe (3)", "Spawned (4)")
print(barplot(t(prop.table(CT.D.MALSTAGE, 1)), main="Male", xlab=NA, ylab=NA, las=1, col=c("#eff3ff","#6baed6","#3182bd","#08519c","#bdd7e7"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-.5, 0), title="Gonad Stage", cex=1.5)))
chisq.test(CT.D.MALSTAGE, simulate.p.value = T, B = 10000) 

colnames(CT.D.fem.stage) <- c("None present (0)", "Early (1)", "Advanced (2)", "Ripe (3)", "Spawned (4)")
print(barplot(t(prop.table(CT.D.fem.stage, 1)), main="Female", xlab=NA, ylab=NA, las=1, col=c("#f1eef6","#df65b0","#dd1c77","#980043","#d7b5d8"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-.5, 0), title="Gonad Stage", cex=1.5)))
chisq.test(CT.D.fem.stage, simulate.p.value = T, B = 10000) 

mtext("10°C - High Nutrition", side=3, line = .5, outer = T, cex = 2, col = "gray30")
mtext("% Sampled", side=2, line = 1, outer = T, cex = 1.5, col = "gray30")
mtext("Week Sampled", side=1, line = 2.5, outer = T, cex = 1.5, col = "gray30")

dev.off()


# PLOTS FOR LEGENDS 

pdf(file = "results/gonad-barplots-sex-legend.pdf", width = 8.5, height = 6.5)
par(mar=c(1, 2, 3, 20), col="gray30")
colnames(CT.D.sex) <- c("Undifferentiated", "Male", "Male dominant", "Hermaphroditic", "Female dominant", "Female")
print(barplot(t(prop.table(CT.D.sex, 1)), main="Dominant gonad sex", xlab="Week Sampled", ylab="% Sampled", las=1, col=c("gray75", "#08519c", "#3182bd", "purple3","mediumorchid3", "#df65b0"),  cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.3, cex.names = 1.3, legend.text = T, args.legend = list(x = "topright", bty = "n", inset=c(-1, 0), title="Gonad Sex", cex=1.5)))
dev.off()

pdf(file = "results/gonad-barplots-male-stage-legend.pdf", width = 8.5, height = 6.5)
par(mar=c(1, 2, 3, 20), col="gray30")
print(barplot(t(prop.table(CT.D.MALSTAGE, 1)), main="Male", xlab=NA, ylab=NA, las=1, col=c("#eff3ff","#6baed6","#3182bd","#08519c","#bdd7e7"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.3, cex.names = 1.3, legend.text = T, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-1, 0), title="Male Gonad Stage", cex=1.5)))
dev.off()

pdf(file = "results/gonad-barplots-female-stage-legend.pdf", width = 8.5, height = 6.5)
par(mar=c(1, 2, 3, 20), col="gray30")
print(barplot(t(prop.table(CT.D.fem.stage, 1)), main="Female", xlab=NA, ylab=NA, las=1, col=c("#f1eef6","#df65b0","#dd1c77","#980043","#d7b5d8"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.3, cex.names = 1.3, legend.text = T, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-1, 0), title="Female Gonad Stage", cex=1.5)))
dev.off()
