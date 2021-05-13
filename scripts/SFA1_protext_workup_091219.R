# Cassi Wattenburger
# 9/12/19
# SFA1 soil wash sample protein extraction workup

rm(list=ls())

setwd("/home/cassi/SFAphase1/protein_extractions")

library("tidyverse")
library("broom")

# Import raw data

raw1 <- read.csv("SFA1_prot_1.1-20.2_091219_edit.csv")
raw2 <- read.csv("SFA1_prot_21.1-40.2_091219_edit.csv")
raw3 <- read.csv("SFA1_prot_41.1-cntrls_091219_edit.csv")

raw.list <- list(raw1, raw2, raw3)

# Substract blank (Standard I) average from other standards and unknowns

blank.list <- list()
  
for (i in 1:length(raw.list)) {
  blank.list[[i]] <- filter(raw.list[[i]], ID=="I")
}

for (i in 1:length(blank.list)) {
  blank.list[[i]] <- summarize(blank.list[[i]], avg=mean(OD_Values))
}

for (i in 1:length(raw.list)) {
  raw.list[[i]] <- mutate(raw.list[[i]], corrected=OD_Values-blank.list[[i]]$avg)
}

# Average technical replicates for protein conc. measurement

raw.summary <- list()

for (i in 1:length(raw.list)) {
  raw.summary[[i]] <- raw.list[[i]] %>% group_by(Sample, ID, Type) %>% summarize(avg=mean(corrected), sd=sd(corrected))
}

###################
# Standards work up

std.list <- list()

for (i in 1:length(raw.summary)) {
  std.list[[i]] <- subset(raw.summary[[i]], Type=="Standard")
}

for (i in 1:length(std.list)) { 
  std.list[[i]]["Actual"] <- c(1500, 1000, 750, 500, 250, 125, 25, 0) # actual conc. values of standards
}

# Std linear regressions

model.list <- list() 

for (i in 1:length(std.list)) {
  model.list[[i]] <- lm(avg ~ Actual, data=std.list[[i]])
}

r2.list <- list() # store R2 of models

for (i in 1:length(model.list)) {
  r2.list[[i]] <- summary(model.list[[i]])$r.squared 
}

for (i in 1:length(model.list)) {
  model.list[[i]] <- as.data.frame(tidy(model.list[[i]])) # store intercept and slope
}

# Graph stds

title.list <- list("Standards for samples 1-20", "Standards for samples 21-40", "Standards for samples 41-51 and controls")

stdgraph.list <- list()

for (i in 1:length(std.list)) {
  stdgraph.list[[i]] <- ggplot(std.list[[i]], aes(x=Actual, y=avg)) +
    geom_point() +
    geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd)) +
    geom_smooth(linetype=2) +
    labs(x="Measurement", y="Actual conc. (ug/mL)", title=title.list[[i]]) +
    geom_text(aes(label=ID, hjust=-0.5, vjust=1)) +
    annotate("text", x=50, y=2.1, label="R2=") +
    annotate("text", x=200, y=2, label=r2.list[[i]])
}

##################
# Unknowns work up

unk.list <- list()

for (i in 1:length(raw.summary)) {
  unk.list[[i]] <- filter(raw.summary[[i]], Type=="Sample" | Type=="Positive" | Type=="Negative")
}

# calculate conc. protein based on standards

protconc.list <- list()

for (i in 1:length(unk.list)) {
  protconc.list[[i]] <- mutate(unk.list[[i]], Conc = (avg-model.list[[i]]$estimate[1])/model.list[[i]]$estimate[2]) # DOUBLE CHECK THIS
}

# average extraction duplicates

for (i in 1:length(protconc.list)) {
  protconc.list[[i]] <- protconc.list[[i]] %>% group_by(ID, Type) %>% summarize(concavg=mean(Conc), sd=sd(Conc))
}

# combine data frames

prot.data <- do.call('rbind', protconc.list)
prot.data$ID <- as.numeric(prot.data$ID)

# Average positive and negative controls
negavg <- mean(filter(prot.data, Type=="Negative")$concavg)
posavg <- mean(filter(prot.data, Type=="Positive")$concavg)

prot.data2 <- filter(prot.data, Type=="Sample")
prot.data2$ID <- as.numeric(prot.data2$ID)


# graph

prot.graph <- ggplot(filter(prot.data2, Type=="Sample"), aes(x=ID, y=concavg)) +
  geom_point() +
  geom_errorbar(aes(ymin=concavg-sd, ymax=concavg+sd), width=0.2) +
  labs(title="Sample protein concentrations", y="Conc. (ug/mL)", x="Sample ID") +
  geom_hline(yintercept=negavg, linetype=2) +
  geom_hline(yintercept=posavg) +
  scale_x_continuous(breaks = seq(1, 51, by = 1)) +
  theme_bw()
prot.graph  


###############
# Normalization
# Samples will be diluted so that final protein concentrations are approximately equal

# remove sample 51 (protein extraction blank)
prot.data3 <- prot.data2[!prot.data2$ID==51,]
  
# lowest sample concentration
min <- tapply(prot.data3$concavg, prot.data3$Type, min)
min <- as.numeric(min)

min.sample <- filter(prot.data3, concavg==min)
min.sample

