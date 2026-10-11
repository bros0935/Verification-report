
## Get functions and data
source("~/Documents/Past subjects Psychology/PSYC3010/R_files/functions3010.R")
JH_data <- read_sav("/home/ben/Documents/interdisciplinary project/James Holmes_Expt 2_DATA.sav")
JH_data <- as.data.frame(JH_data)


##Screen and recode data
# Condition: 1 = No-Task Control, 2 = Reactivation+Tetris, 3 = Tetris Only, 4 = Reactivation Only
JH_data$id <- 1:nrow(JH_data)
JH_data$Condition <- as.numeric(JH_data$Condition)
JH_data$Diary <- JH_data$Days_One_to_Seven_Image_Based_Intrusions_in_Intrusion_Diary
JH_data$IPT <- JH_data$Intrusion_Provocation_Task_Intrusions

descript_hist(JH_data[, c("Diary", "IPT")])

JH_data$Tetris <- factor(ifelse(JH_data$Condition %in% c(2, 3), 2, 1), levels=c(1, 2),
                         labels = c("No Tetris", "Tetris"))
JH_data$Reactivation <- factor(ifelse(JH_data$Condition %in% c(2, 4), 2, 1), levels=c(1, 2),
                               labels = c("No Reactivation", "Reactivation"))

table(JH_data$Tetris, JH_data$Reactivation)


#####    DV 1: Diary intrusions (Days 1-7)    ####

## Look at scores, graphs etc
aggregate(Diary ~ Tetris + Reactivation, JH_data, function(x)
  c(mean=mean(x), sd=sd(x), n=length(x)))

aggregate(Diary ~ Tetris, JH_data, function(x)
  c(mean=mean(x), sd=sd(x), n=length(x)))

aggregate(Diary ~ Reactivation, JH_data, function(x)
  c(mean=mean(x), sd=sd(x), n=length(x)))

mean(JH_data$Diary)
sd(JH_data$Diary)

aov_barplot(Diary ~ Tetris + Reactivation, JH_data)

## Conduct ANOVA
Diary_anova <- aov_ez(
  id = "id",
  dv = "Diary",
  between = c("Tetris","Reactivation"),
  data = JH_data,
  anova_table = list(es = "pes")
)
Diary_anova


#####    DV 2: IPT intrusions (Day 7)    ####

## Look at scores, graphs etc
aggregate(IPT ~ Tetris + Reactivation, JH_data, function(x)
  c(mean=mean(x), sd=sd(x), n=length(x)))

aggregate(IPT ~ Tetris, JH_data, function(x)
  c(mean=mean(x), sd=sd(x), n=length(x)))

aggregate(IPT ~ Reactivation, JH_data, function(x)
  c(mean=mean(x), sd=sd(x), n=length(x)))

mean(JH_data$IPT)
sd(JH_data$IPT)

aov_barplot(IPT ~ Tetris + Reactivation, JH_data)

## Conduct ANOVA
IPT_anova <- aov_ez(
  id = "id",
  dv = "IPT",
  between = c("Tetris","Reactivation"),
  data = JH_data,
  anova_table = list(es = "pes")
)
IPT_anova
