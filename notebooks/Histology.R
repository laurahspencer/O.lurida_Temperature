Histology <- read.csv("data/histology.csv", header=T, stringsAsFactors = T, na.strings = c("NA", "TBD", ""))
str(Histology)
Histology$FEMALE.STAGE <- as.factor(Histology$FEMALE.STAGE)
Histology$MALE.STAGE <- as.factor(Histology$MALE.STAGE)

### Stats
# Resource: http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r
library("gplots")
df <- droplevels(df)

########## Compare sex and dominant stage between temperatures at February sampling 
levels(Histology$SEX)

CT.Date.SEX <- table(Histology$Date, Histology$SEX)
CT.Treat.SEX <- table(Histology$TREAT, Histology$SEX)
chisq.test(CT.Treat.SEX, simulate.p.value = T)

CT.Date.STAGE.F <- table(Histology$Date, Histology$FEMALE.STAGE)
CT.Treat.STAGE.F <- table(Histology$TREAT, Histology$FEMALE.STAGE)
chisq.test(CT.Treat.STAGE.F, simulate.p.value = T)

CT.Date.STAGE.M <- table(Histology$Date, Histology$MALE.STAGE)
CT.Treat.STAGE.M <- table(Histology$TREAT, Histology$MALE.STAGE)
chisq.test(CT.Treat.STAGE.M, simulate.p.value = T)

CT.Date.Atresia <- table(Histology$Date, Histology$SIGNIFICANT.ATRESIA)
CT.Treat.Atresia <- table(Histology$SIGNIFICANT.ATRESIA, Histology$TREAT)
chisq.test(CT.Treat.Atresia)

# Contingency table, showing stages in 10C and 6C treatment groups 
balloonplot(t(CT.Treat.Atresia), label = T, show.margins = FALSE)

summary(Histology$SEX)
131/(131+33+34+21+3+17) # % female only
17/(131+33+34+21+3+17) # % male only
(131+34)/(131+33+34+21+3+17) # % female + HPF
(17+21)/(131+33+34+21+3+17) # % male + HPM

