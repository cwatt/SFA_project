# SFA Phase 1 - GCMS data
# Cassi Wattenburger
# Started 10/30/19

library("ggplot2")
library("stringr")
library("plyr")
library("dplyr")
library("broom")

rm(list=ls())

setwd("C:/Users/Cassi/Google Drive/08 Research/07 SFA experiments 2019/Phase1/GCMS/data")

###############################################################
# Import data

# File paths
dir <- list.dirs(path=getwd(), full.names=TRUE, recursive=FALSE)

files <- list()

for (i in dir) {
  files[[i]] <- list.files(path=i, pattern="*.txt", full.names=TRUE, recursive=FALSE)
}
# NOTE: foo files were removed from directories (do not contain data)

# Exclude data from 11/01/19, sampling was compromised by plugged needles
files <- files[names(files) != "C:/Users/Cassi/Google Drive/08 Research/07 SFA experiments 2019/Phase1/GCMS/data/110119"] 

# Extract area data from each file
name.list <- list()
area.list <- list() 

for (i in names(files)) {
  for (x in files[[i]]) {
    name <- str_extract(x, "data/.+/.+.txt") # not catching those dates for some reason?
    name <- gsub("data/.+/", "", name)
    name <- gsub(".txt", "", name)
    name.list[[i]] <- append(name.list[[i]], name)
    
    t <- read.table(x, header=FALSE, fill=TRUE, sep="\t", skip=8, nrows=3, comment.char="", stringsAsFactors=FALSE) # error, no lines
    t <- t[1,10] # TIC area 
    t <- as.data.frame(t)
    area.list[[i]] <- append(area.list[[i]], t)
  }
}

# Convert into dataframes
raw.data <- list()

for (i in names(files)) {
  raw.data[[i]] <- cbind(area.list[[i]])  
  rownames(raw.data[[i]]) <- name.list[[i]]
  #raw.data[[i]] <- cbind(raw.data[[i]])
  colnames(raw.data[[i]]) <- "area"
}

# Date info
date <- list()

for (i in names(files)) {
  date[[i]] <- str_extract(i, "data/.+")
  date[[i]] <- gsub("data/", "", date[[i]])
}

############################################################################
# Standards

# Extract standards from each time point
stds <- list()

for (i in names(raw.data)) {
  match <- c("std.*")
  stds[[i]] <- raw.data[[i]][rownames(raw.data[[i]]) %in% grep(paste(match, collapse="|"), rownames(raw.data[[i]]), value=TRUE),]
  stds[[i]] <- cbind(stds[[i]])
  colnames(stds[[i]]) <- "area"
  stds[[i]] <- mutate(data.frame(stds[[i]]), std = rownames(data.frame(stds[[i]]))) # add std names as column
  stds[[i]] <- mutate(data.frame(stds[[i]]), date = date[[i]]) # add date info
}

# Smoosh all stds into one dataframe
stds.df <- data.frame()

for (i in names(stds)) {
  stds.df <- rbind(stds.df, stds[[i]])
}

# Add ppm info to stds dataframe
stds.df$ppm <- ifelse(stds.df$std == "std0", 0, ifelse(stds.df$std == "std1", 1000, ifelse(stds.df$std == "std2", 2000,
                            ifelse(stds.df$std == "std3", 4000, ifelse(stds.df$std == "std4" | stds.df$std == "std4.1" | stds.df$std == "std4.2" |
                            stds.df$std == "std4.3" | stds.df$std == "std4.4" | stds.df$std == "std4.5" | stds.df$std == "std4.6", 10000,
                            ifelse(stds.df$std == "std5", 20000, ifelse(stds.df$std == "std6", 40000, ifelse(stds.df$std == "std7", 100000, 
                            ifelse(stds.df$std == "std8", 300000, ifelse(stds.df$std == "std9", 500000, NA))))))))))

# Remove stds 8-9 up to 11/09/19 and all std 10 (oversaturated GCMS)
stds.df <- stds.df[!stds.df$std %in% "std10",]
stds.df <- stds.df[!((stds.df$date < 111119 & stds.df$std=="std8") | (stds.df$date < 111119 & stds.df$std=="std9")),]

# Adjust for 90% split ratio stds 8 and 9 (too concentrated for GCMS to handle), other stds run at 50% split ratio
# multiplier = 1.8, from 0.5x=0.9
stds.df$area <- as.numeric(stds.df$area)
stds.df <- mutate(stds.df, area.conv = ifelse(std=="std8", area*1.8, ifelse(std=="std9", area*1.8, area*1)))

# Isolate replicates of std 4
stds4.df <- stds.df[stds.df$std %in% grep("std4+", stds.df$std, value=TRUE),]

# Remove std 4 reps from stds data for regression
stds.df <- stds.df[!stds.df$std %in% grep("std4.+", stds.df$std, value=TRUE),]


# Graph standards

# by date
stds.date.graph <- ggplot(stds.df, aes(x=ppm, y=area.conv)) +
  geom_point() +
  facet_wrap(~date) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label=std),hjust=0, vjust=0)
stds.date.graph

# 11/19/19 and 11/26/19 are outliers

# all together
stds.all.graph <- ggplot(stds.df, aes(x=ppm, y=area.conv)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label=std),hjust=0, vjust=0)
stds.all.graph

# Identify outliers

stds.out.graphs <- list()

for (i in unique(stds.df$std)) {
  stds.out.graphs[[i]] <- ggplot(stds.df[stds.df$std==i,], aes(x=ppm, y=area.conv)) +
    geom_boxplot() +
    facet_grid(.~std, scales="free") +
    geom_text(aes(label=date), hjust=1, vjust=0, size=3)
}

stds.nona.df <- na.omit(stds.df)
stds.summary.df <- ddply(stds.nona.df, .(std), summarize, mean=mean(area.conv), sd=sd(area))

stds.summary.df <- mutate(stds.summary.df, cutoff=sd*2) # defining outliers outside 2 standard deviations of mean
  
stds.summary.df <- mutate(stds.summary.df, cutoff_min=mean-cutoff)
stds.summary.df <- mutate(stds.summary.df, cutoff_max=mean+cutoff)

stds.rm.df <- data.frame()

for (i in stds.summary.df$std) {
  x <- filter(stds.nona.df[stds.nona.df$std==i,], !(area.conv < stds.summary.df$cutoff_min[stds.summary.df$std==i] | area.conv > stds.summary.df$cutoff_max[stds.summary.df$std==i])) # remove outliers
  stds.rm.df <- rbind(x, stds.rm.df)
}
# removed 16 outliers

# Check for systematic trends in standards
day.anova <- summary(aov(area.conv ~ date, stds.rm.df))
# P = 0.6, affect of day on std data not significant

# graph without outliers
stds.rm.graph <- ggplot(stds.rm.df, aes(x=ppm, y=area.conv)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label=std),hjust=0, vjust=0)
stds.rm.graph

# Linear regression
stds.conv <- lm(area.conv ~ ppm, data=stds.rm.df) #TROUBLESHOOOT intercept of this looks funky?????????
summary(stds.conv)
# R-squared = 0.985

plot(area.conv ~ ppm, data=stds.rm.df, main="Calibration curve") + abline(stds.conv)

stds.conv.df <- as.data.frame(tidy(stds.conv))

# Look at std 4 replicate data
stds4.graph <- ggplot(stds4.df, aes(x=date, y=area)) +
  geom_point() +
  geom_text(aes(label=std),hjust=0, vjust=0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
stds4.graph
# pretty good except for the one outlier, don't trust to use as drift check though, will have to rely on randomized design instead

# Ambient data
amb <- list()
match <- c("ambient")

for (i in names(raw.data)) {
  amb[[i]] <- raw.data[[i]][rownames(raw.data[[i]]) %in% grep(paste(match, collapse="|"), rownames(raw.data[[i]]), value=TRUE),]
  amb[[i]] <- cbind(amb[[i]])
  colnames(amb[[i]]) <- "area"
  amb[[i]] <- mutate(data.frame(amb[[i]]), date = date[[i]]) # add date info
}

amb.df <- data.frame()

for (i in names(amb)) {
  amb.df <- rbind(amb.df, amb[[i]])
}

amb.df$area <- as.numeric(amb.df$area)

# Remove outliers (due to cross-contamination between samples before I realized it was a problem)

amb.avg <- mean(amb.df$area)
amb.sd <- sd(amb.df$area)
amb_min <- amb.avg-(2*amb.sd)
amb_max <- amb.avg+(2*amb.sd)

amb.rm.df <- filter()

amb.rm.df <- filter(amb.df, !(amb.df$area < amb_min | amb.df$area > amb_max))
# removed one outlier

# Calculate ambient air value
amb.value <- mean(amb.rm.df$area)

# Graph ambient data
amb.graph <- ggplot(amb.rm.df, aes(x=date, y=area)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
amb.graph


###############################################
# Samples

# Extract samples from each time point
samples <- list()

for (i in names(raw.data)) {
  match <- c("std.*", "ambient")
  samples[[i]] <- raw.data[[i]][!rownames(raw.data[[i]]) %in% grep(paste(match, collapse="|"), rownames(raw.data[[i]]), value=TRUE),]
  samples[[i]] <- cbind(samples[[i]])
  colnames(samples[[i]]) <- "area"
  samples[[i]] <- mutate(data.frame(samples[[i]]), id = rownames(data.frame(samples[[i]]))) # add sample names as column
  samples[[i]] <- mutate(data.frame(samples[[i]]), date = date[[i]]) # add date info
}

# Smoosh all samples into one dataframe
samples.df <- data.frame()

for (i in names(samples)) {
  samples.df <- rbind(samples.df, samples[[i]])
}

# NOTE: samples 5 and 51 on 11/11/19 were accidentally switched, labelling was corrected on .txt files prior to import

# Add replicate, sample number, day metadata

# split into replicates
rep1 <- c(1:51)
rep2 <- c(52:102)
rep3 <- c(103:153)

rep1.df <- samples.df[samples.df$id %in% rep1,]
rep2.df <- samples.df[samples.df$id %in% rep2,]
rep3.df <- samples.df[samples.df$id %in% rep3,]

# replicates
rep1.df$id <- as.numeric(rep1.df$id)
rep1.df <- rep1.df[order(rep1.df$id),]
rep1.df <- rep1.df[order(rep1.df$date),]

rep2.df$id <- as.numeric(rep2.df$id)
rep2.df <- rep2.df[order(rep2.df$id),]
rep2.df <- rep2.df[order(rep2.df$date),]

rep3.df$id <- as.numeric(rep3.df$id)
rep3.df <- rep3.df[order(rep3.df$id),]
rep3.df <- rep3.df[order(rep3.df$date),]

rep1.df$rep <- paste(rep(1, nrow(rep1.df)))
rep2.df$rep <- paste(rep(2, nrow(rep2.df)))
rep3.df$rep <- paste(rep(3, nrow(rep3.df)))

# sample number
################################### be careful of labelling here
samplenum <- c(1:51)
rep1.df$sample <- paste(rep(samplenum, length(unique(rep1.df$date))))
rep2.df$sample <- paste(rep(samplenum, length(unique(rep2.df$date))))
rep3.df$sample <- paste(rep(samplenum, length(unique(rep3.df$date))))

# day
rep1.df <- mutate(rep1.df, day = ifelse(rep1.df$date==102619, 0, ifelse(rep1.df$date==102719, 1, ifelse(rep1.df$date==102919, 3, 
                                ifelse(rep1.df$date==103119, 5, ifelse(rep1.df$date==110219, 7, ifelse(rep1.df$date==110419, 9, 
                                ifelse(rep1.df$date==110719, 12, ifelse(rep1.df$date==111119, 17, ifelse(rep1.df$date==111419, 19, 
                                ifelse(rep1.df$date==111819, 23, ifelse(rep1.df$date==112119, 26, ifelse(rep1.df$date==112519, 30, NA)))))))))))))

rep2.df <- mutate(rep2.df, day = ifelse(rep2.df$date==102719, 0, ifelse(rep2.df$date==102819, 1, ifelse(rep2.df$date==103019, 3, 
                                ifelse(rep2.df$date==110119, 5, ifelse(rep2.df$date==110319, 7, ifelse(rep2.df$date==110519, 9,
                                ifelse(rep2.df$date==110819, 12, ifelse(rep2.df$date==111219, 17, ifelse(rep2.df$date==111519, 19, 
                                ifelse(rep2.df$date==111919, 23, ifelse(rep2.df$date==112219, 26, ifelse(rep2.df$date==112619, 30, NA)))))))))))))

rep3.df <- mutate(rep3.df, day = ifelse(rep3.df$date==102819, 0, ifelse(rep3.df$date==102919, 1, ifelse(rep3.df$date==103119, 3, 
                                ifelse(rep3.df$date==110219, 5, ifelse(rep3.df$date==110419, 7, ifelse(rep3.df$date==110619, 9, 
                                ifelse(rep3.df$date==110919, 12, ifelse(rep3.df$date==111319, 17, ifelse(rep3.df$date==111619, 19,
                                ifelse(rep3.df$date==112019, 23, ifelse(rep3.df$date==112319, 26, ifelse(rep3.df$date==112719, 30, NA)))))))))))))

# Rejoin data
samples.df <- rbind(rep1.df, rep2.df, rep3.df)

samples.df$sample <- as.factor(samples.df$sample)
samples.df$rep <- as.factor(samples.df$rep)
samples.df$date <- as.factor(samples.df$date)
samples.df$area <- as.numeric(samples.df$area)

# Distinguish negative controls
samples.df <- mutate(samples.df, type = ifelse(samples.df$sample==51, "negative", "sample"))

# Apply conversion to ppm
samples.df <- mutate(samples.df, ppm = (area - stds.conv.df[1,2])/stds.conv.df[2,2])

# Apply conversion to avg ambient
amb.value.ppm <- (amb.value-stds.conv.df[1,2])/stds.conv.df[2,2]

# Convert to mol CO2/mL
samples.df <- mutate(samples.df, molCO2=ppm/1000) # 1. convert ppm (mg/L) to g/L
samples.df$molCO2 <- samples.df$molCO2/44.009 # 2. convert to mol CO2/L
samples.df$molCO2 <- samples.df$molCO2*0.001 # 3. convert to mol/mL
samples.df$molCO2 <- samples.df$molCO2*XXXX ######################## 4. total mol in space of ucosm













# Average by replicate and day
avg.df <- ddply(samples.df, .(sample, day, type), summarize, mean=mean(ppm), sd=sd(ppm))

# Graph samples
samples.graph <- ggplot(avg.df, aes(x=sample, y=mean, color=type)) +
  geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd)) +
  facet_grid(~day)
samples.graph

samples.graph2 <- ggplot(samples.df, aes(x=sample, y=area, color=type)) +
  geom_point() +
  geom_text(aes(label=sample)) +
  #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd)) +
  facet_grid(rep~day)
samples.graph2










# Bottle volumes (mL)
bottle.empty <- 37.3
bottle.sand <- 34.17

# Function to convert ppm into mol C ################## want someone else to double check this calc (does new volume matter or just pressure????)
calc_mol_C = function(ppm,ion,volume.L) {
  # moles (n) = PV / RT
  
  temp.K=294.261
  pressure.atm=0.94 # (account for shared volume between serum bottle and sample vial)
  R=0.08206
  
  ppm = as.numeric(ppm)
  ion = as.numeric(ion)
  mol.volume = (pressure.atm * volume.L) / (R * temp.K)  # mol gas in container
  mol.CO2 = mol.volume * (ppm / 1000000)  
  mol.C <- mol.CO2  # fraction of mol as CO2 ~ mol C
  
  return(mol.C)
}


