## Get functions and data
source("~/Documents/Past subjects Psychology/PSYC3010/R_files/functions3010.R")
data <- read_sav("/home/ben/Documents/interdisciplinary project/James Holmes_Expt 2_DATA.sav")
data <- as.data.frame(JH_data)

data$id <- 1:nrow(data)
data$Condition <- as.numeric(data$Condition)
data$Diary <- data$Days_One_to_Seven_Image_Based_Intrusions_in_Intrusion_Diary
data$IPT <- data$Intrusion_Provocation_Task_Intrusions

#Find tetris in T only or T+R
#2 and three into true false, then true false into 1, 2
data$Tetris <- factor(ifelse(data$Condition %in% c(2, 3), 2, 1), levels=c(1, 2),
                      labels = c("No Tetris", "Tetris"))

#Same for R or T+R
data$Reactivation <- factor(ifelse(data$Condition %in% c(2, 4), 2, 1), levels=c(1, 2),
                      labels = c("No Reactivation", "Reactivation"))



## Look at scores, graphs etc
aggregate(Diary ~ Tetris + Reactivation, data, function(x)
  c(mean=mean(x), sd=sd(x), n=length(x)))

aggregate(Diary ~ Tetris, data, function(x)
  c(mean=mean(x), sd=sd(x), n=length(x)))

aggregate(Diary ~ Reactivation, data, function(x)
  c(mean=mean(x), sd=sd(x), n=length(x)))

mean(data$Diary)
sd(data$Diary)

aov_barplot(Diary ~ Tetris + Reactivation, data)

#good one
aov_barplot(Diary ~ Tetris + Reactivation, data)

## Conduct ANOVA
Diary_anova <- aov_ez(
  id = "id",
  dv = "Diary",
  between = c("Tetris","Reactivation"),
  data = data,
  anova_table = list(es = "pes")
)

Diary_anova

##Do again for IPT
## Look at scores, graphs etc
aggregate(IPT ~ Tetris + Reactivation, data, function(x)
  c(mean=mean(x), sd=sd(x), n=length(x)))

aggregate(IPT ~ Tetris, data, function(x)
  c(mean=mean(x), sd=sd(x), n=length(x)))

aggregate(IPT ~ Reactivation, data, function(x)
  c(mean=mean(x), sd=sd(x), n=length(x)))

mean(data$IPT)
sd(data$IPT)

aov_barplot(IPT ~ Tetris + Reactivation, data)

#good one
aov_barplot(IPT ~ Tetris + Reactivation, data)

## Conduct ANOVA
IPT_anova <- aov_ez(
  id = "id",
  dv = "IPT",
  between = c("Tetris","Reactivation"),
  data = data,
  anova_table = list(es = "pes")
)

IPT_anova

