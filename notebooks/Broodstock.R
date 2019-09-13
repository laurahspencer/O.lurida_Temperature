# Adult size and mortality 

str(histo) # use "histo" dataframe, imported in the Gonad-histology.R script 

mean(subset(histo, TREAT!="Wild")$Length..cm., na.rm=TRUE)*10
sd(subset(histo, TREAT!="Wild")$Length..cm., na.rm=TRUE)*10

size.adult <- histo %>%
  filter(!is.na(TREAT)) %>%
  group_by(TREAT, Week) %>%
  summarize(mean_length = mean(Length..cm., na.rm = TRUE), sd_length = sd(Length..cm., na.rm = TRUE),
            mean_weight = mean(Est..Tissue.Weight..g., na.rm = TRUE), sd_weight = sd(Est..Tissue.Weight..g., na.rm = TRUE))
size.adult[(nrow(size.adult)+1):(nrow(size.adult)+5),] <- size.adult[1,]
size.adult[(nrow(size.adult)-4):(nrow(size.adult)),"TREAT"] <- as.factor(c("A", "B", "C", "D", "Wild"))

# subset(size.adult, Week!=0 & TREAT!="Wild")
# plot mean weight over time 
pdf("results/broodstock-weight.pdf", width = 7, height = 3.5)
ggplot(data=subset(size.adult, TREAT!="PRE" & mean_weight<3.5), aes(x=Week, y=mean_weight, group=TREAT, col=TREAT)) + geom_line() + theme_bw(base_size = 12) + geom_point(size=3) + ylab("Mean weight (g)") + ggtitle(label = "Mean wet tissue weight in broodstock") + theme(axis.title.x = element_blank()) +
  scale_color_manual(values=c("#92c5de",
                        "#ca0020","#0571b0","#f4a582", "gray30"),
                        name="Treatment",
                        breaks=c("C", "A", "B", "D", "Wild"),
                        labels=c("7°C+high-food", "7°C+low-food", "10°C+high-food", "10°C+low-food", "Wild")) + scale_x_discrete(labels= c("Nov 30", "Dec 20", "Jan 4", "Jan 23", "Feb 9", "Feb 27", "Mar 3", "Mar 23")) + 
  geom_vline(xintercept = 4.1, linetype="solid", color = "gray50", size=.5) + 
  geom_vline(xintercept = 6.1, linetype="dashed", color = "gray50", size=.5)
#+geom_errorbar(aes(ymin=mean_weight-sd_weight, ymax=mean_weight+sd_weight), width=.1)
dev.off()

# Does weight change over time, diff by treat? 
anova(lm(Est..Tissue.Weight..g. ~ as.numeric(Week)+factor(TREAT), data=subset(histo, TREAT!="Wild")))
summary(lm(Est..Tissue.Weight..g. ~ as.numeric(Week)+factor(TREAT), data=subset(histo, TREAT!="Wild"))) # overall, tissue weight same by treat.
anova(lm(Est..Tissue.Weight..g. ~ as.numeric(Week), data=subset(histo, TREAT!="Wild" & TREAT!="PRE"))) 
summary(lm(Est..Tissue.Weight..g. ~ as.numeric(Week), data=subset(histo, TREAT!="Wild" & TREAT!="PRE"))) 
# weight = 3.28 - 0.0929x (where x=week)

# did wild change?
anova(lm(Est..Tissue.Weight..g. ~ as.numeric(Week), data=subset(histo, TREAT=="Wild"))) 
summary(lm(Est..Tissue.Weight..g. ~ as.numeric(Week), data=subset(histo, TREAT=="Wild"))) 

# plot mean length over time 

# All
ggplot(data=size.adult, aes(x=Week, y=mean_length, group=TREAT, col=TREAT)) +
  geom_errorbar(aes(ymin=mean_length-sd_length, ymax=mean_length+sd_length), width=.1) +
  geom_line()+
  geom_point()

ggplot(data=subset(histo, TREAT=="A"), aes(x=Week, y=Est..Tissue.Weight..g.)) +
  geom_boxplot(fill="seagreen3") + ggtitle(label = "Wet tissue weight (g), Cold-Low") +
  ylim(1, 7)  

ggplot(data=subset(histo, TREAT=="B"), aes(x=Week, y=Est..Tissue.Weight..g.)) +
  geom_boxplot(fill="orange1") + ggtitle(label = "Wet tissue weight (g), Warm-Low") +
  ylim(1, 7)  

ggplot(data=subset(histo, TREAT=="C"), aes(x=Week, y=Est..Tissue.Weight..g.)) +
  geom_boxplot(fill="skyblue3") + ggtitle(label = "Wet tissue weight (g), Cold-High") +
  ylim(1, 7)  

ggplot(data=subset(histo, TREAT=="D"), aes(x=Week, y=Est..Tissue.Weight..g.)) +
  geom_boxplot(fill="indianred2") + ggtitle(label = "Wet tissue weight (g), Warm-High") +
  ylim(1, 7)

ggplot(data=subset(histo, TREAT=="Wild"), aes(x=Week, y=Est..Tissue.Weight..g.)) +
  geom_boxplot(fill="gray50") + ggtitle(label = "Wet tissue weight (g), Mud Bay") +
  ylim(1, 7)  


plot_ly(data = size.adult, x = ~Week, y = ~mean_weight, type="scatter", mode="markers", marker=list(size=14), color=~TREAT, colors = c("skyblue3", "seagreen3", "indianred2",  "orange1"), hovertext=~TREAT) %>%  #generate plotly plot
  layout(title="Shell length",
         yaxis = list(title = 'shell length', size=26), xaxis=list(title="Week", size=30),
         legend = list(x=.05, y=.95))



# Does length change over time, diff by treat? 
anova(lm(Length..cm. ~ as.numeric(Week)+factor(TREAT), data=subset(histo, TREAT!="Wild" & TREAT!="PRE")))
summary(lm(Length..cm. ~ as.numeric(Week)+factor(TREAT), data=subset(histo, TREAT!="Wild" & TREAT!="PRE"))) # overall, tissue weight same by treat.
summary(lm(Length..cm. ~ as.numeric(Week), data=subset(histo, TREAT!="Wild" & TREAT!="PRE"))) 
anova(lm(Length..cm. ~ as.numeric(Week), data=subset(histo, TREAT!="Wild" & TREAT!="PRE"))) 

# wild?
summary(lm(Length..cm. ~ as.numeric(Week), data=subset(histo, TREAT=="Wild"))) 


# plot individual weeks 

# Week 3
ggplot(data=subset(histo, Week==3), aes(x=TREAT, y=Length..cm., col=TREAT)) + geom_boxplot() + ggtitle(label="Week 3 Length")

ggplot(data=subset(histo, Week==3), aes(x=TREAT, y=Est..Tissue.Weight..g., col=TREAT)) + geom_boxplot() + ggtitle(label="Week 3 Tissue Weight")

# Week 5
ggplot(data=subset(histo, Week==5), aes(x=TREAT, y=Length..cm., col=TREAT)) + geom_boxplot() + ggtitle(label="Week 5 Length")

ggplot(data=subset(histo, Week==5), aes(x=TREAT, y=Est..Tissue.Weight..g., col=TREAT)) + geom_boxplot() + ggtitle(label="Week 5 Tissue Weight")

# Week 8
ggplot(data=subset(histo, Week==8), aes(x=TREAT, y=Length..cm., col=TREAT)) + geom_boxplot() + ggtitle(label="Week 8 Length")

ggplot(data=subset(histo, Week==8), aes(x=TREAT, y=Est..Tissue.Weight..g., col=TREAT)) + geom_boxplot() + ggtitle(label="Week 8 Tissue Weight")

# Week 10
ggplot(data=subset(histo, Week==10), aes(x=TREAT, y=Length..cm., col=TREAT)) + geom_boxplot() + ggtitle(label="Week 10 Length")

ggplot(data=subset(histo, Week==10), aes(x=TREAT, y=Est..Tissue.Weight..g., col=TREAT)) + geom_boxplot() + ggtitle(label="Week 10 Tissue Weight")

# Week 13
ggplot(data=subset(histo, Week==13), aes(x=TREAT, y=Length..cm., col=TREAT)) + geom_boxplot() + ggtitle(label="Week 13 Length")

ggplot(data=subset(histo, Week==13), aes(x=TREAT, y=Est..Tissue.Weight..g., col=TREAT)) + geom_boxplot() + ggtitle(label="Week 13 Tissue Weight")

histo$TREAT <- relevel(histo$TREAT, ref = "C")
hist(subset(histo, Week==13)$Est..Tissue.Weight..g.)
shapiro.test(subset(histo, Week==13)$Est..Tissue.Weight..g.)
anova(lm(Est..Tissue.Weight..g. ~ TREAT, data=subset(histo, Week==13)))
summary(lm(Est..Tissue.Weight..g. ~ TREAT, data=subset(histo, Week==13)))


# Week 15
ggplot(data=subset(histo, Week==15), aes(x=TREAT, y=Length..cm., col=TREAT)) + geom_boxplot() + ggtitle(label="Week 15 Length")

ggplot(data=subset(histo, Week==15), aes(x=TREAT, y=Est..Tissue.Weight..g., col=TREAT)) + geom_boxplot() + ggtitle(label="Week 15 Tissue Weight")


# Broodstock mortality 
brood.mortality <- read.csv("data/broodstock-mortality.csv", header=T, stringsAsFactors = F)
brood.mortality$Date <- as.Date(brood.mortality$Date, format = "%m/%d/%y")
brood.mortality <- melt(data = brood.mortality, id.vars = "Date", value.name = "Alive", variable.name = "TREAT.rep")
brood.mortality$TREAT <- brood.mortality$TREAT.rep
brood.mortality$TREAT <- gsub("1", "", brood.mortality$TREAT)
brood.mortality$TREAT <- as.factor(gsub("2", "", brood.mortality$TREAT))


png(file="results/broodstock-survival.png", width=700, height=300)
ggplot(data=brood.mortality, aes(x=Date, y=100*Alive, group=TREAT.rep, col=TREAT)) + theme_bw(base_size = 13) + ggtitle("Broodstock survival over time") + geom_line()+ geom_point() + xlab("Date") + ylab("% Alive") + scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food")) + scale_y_continuous(breaks=c(50, 60, 70, 80, 90, 100))
dev.off()


# Figure out how  to do survival analysis with my % survival data. Maybe just do a binomial glm? 

survdiff(data = brood.mortality, formula = Surv(Alive) ~ TREAT)
