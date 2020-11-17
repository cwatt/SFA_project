# 10/29/19

rm(list=ls())

setwd("C:/Users/Cassi/Google Drive/08 Research/07 SFA experiments 2019/Phase1/water_content")

library("ggplot2")

data <- read.csv("SFAphase1_waterloss.csv")

head(data)
data$rep <- as.factor(data$rep)
class(data$rep)

data$ucosm <-as.factor(data$ucosm)

ucosm1 <- data[data$ucosm %in% 1,]

graph <- ggplot(data,aes(x=date, y=contents)) +
  geom_point() +
  facet_wrap(~ucosm) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
graph

graph2 <- ggplot(data, aes(x=date, y=contents)) +
  geom_boxplot() +
  facet_wrap(~rep) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
graph2
