######
rm(list=ls())
# final plotting from AgentSeal, DEPONS and Tracking
library(ggplot2)
library(scales)
library(gridExtra)
library(rgdal)
library(broom)
library(viridis)
library(sf)
library(dismo)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(tidyverse)
par(mfrow=c(2,6), mar=c(4,4.5,4,1))
options(scipen=999)

setwd("C:/Users/mec21/OneDrive - University of St Andrews/SMRU/DEB_SMRUC/writing/DEB manuscript/Final draft for submission/Submitted to GitHub/Movement/Input")
# hs - agent seal (from agent seal analysis.r)
# loading AgentSeal results for final plotting
load(file="aggMovNHSum.RData")

AShsH30 <- (aggMovNHSum$high30P)
AShsL30 <- (aggMovNHSum$low30P)

AShsH30m <- mean(aggMovNHSum$high30P)
AShsL30m <- mean(aggMovNHSum$low30P)
AShsH30md <- median(aggMovNHSum$high30P)
AShsL30md <- median(aggMovNHSum$low30P)


hsL = hist(AShsL30,plot=FALSE) # or hist(x,plot=FALSE) to avoid the plot of the histogram
hsL$density = hsL$counts/sum(hsL$counts)*100

hsH = hist(AShsH30,plot=FALSE, breaks=seq(0,1,0.1)) # or hist(x,plot=FALSE) to avoid the plot of the histogram
hsH$density = hsH$counts/sum(hsH$counts)*100

# loading porpoise data
load(file="aggProp.RData")
TRhpL30 <- aggProp[,2]
TRhpH30 <- aggProp[,5]

TRhpH30m <- mean(TRhpH30)
TRhpL30m <- mean(TRhpL30)
TRhpH30md <- median(TRhpH30)
TRhpL30md <- median(TRhpL30)

hpL = hist(TRhpL30,plot=FALSE, breaks=seq(0,1,0.1)) # or hist(x,plot=FALSE) to avoid the plot of the histogram
hpL$density = hpL$counts/sum(hpL$counts)*100

hpH = hist(TRhpH30,plot=FALSE, breaks=seq(0,1,0.1)) # or hist(x,plot=FALSE) to avoid the plot of the histogram
hpH$density = hpH$counts/sum(hpH$counts)*100

# hs -tracking
# loading tracking data
load(file="aggMovHS.RData")
TRhsH30 <- (aggMov$pH)
TRhsL30 <- (aggMov$pL)

TRhsH30m <- mean(aggMov$pH)
TRhsL30m <- mean(aggMov$pL)
TRhsH30md <- median(aggMov$pH)
TRhsL30md <- median(aggMov$pL)

hs2L = hist(TRhsL30,plot=FALSE, breaks=seq(0,1,0.1)) # or hist(x,plot=FALSE) to avoid the plot of the histogram
hs2L$density = hsL$counts/sum(hsL$counts)*100

hs2H = hist(TRhsH30,plot=FALSE, breaks=seq(0,1,0.1)) # or hist(x,plot=FALSE) to avoid the plot of the histogram
hs2H$density = hs2H$counts/sum(hs2H$counts)*100

# gs -tracking
# loading data
load(file="aggMovGS.RData")
TRgsH30 <- (aggMov$pH)
TRgsL30 <- (aggMov$pL)

TRgsH30m <- mean(aggMov$pH)
TRgsL30m <- mean(aggMov$pL)
TRgsH30md <- median(aggMov$pH)
TRgsL30md <- median(aggMov$pL)

gsL = hist(TRgsL30,plot=FALSE, breaks=seq(0,1,0.1)) # or hist(x,plot=FALSE) to avoid the plot of the histogram
gsL$density = gsL$counts/sum(gsL$counts)*100

gsH = hist(TRgsH30,plot=FALSE, breaks=seq(0,1,0.1)) # or hist(x,plot=FALSE) to avoid the plot of the histogram
gsH$density = gsH$counts/sum(gsH$counts)*100

# combining low and high on one graph
gs <- cbind.data.frame(c(TRgsH30,TRgsL30),c(rep("High",length(TRgsH30)),rep("Low",length(TRgsL30))) )
colnames(gs) <- c("pdist","area")
gs$species <- "GS"
gs$Source <- "Telemetry"

hs <- cbind.data.frame(c(TRhsH30,TRhsL30),c(rep("High",length(TRhsH30)),rep("Low",length(TRhsL30))) )
colnames(hs) <- c("pdist","area")
hs$species <- "HS"
hs$Source <- "Telemetry"

hsm <- cbind.data.frame(c(AShsH30,AShsL30),c(rep("High",length(AShsH30)),rep("Low",length(AShsL30))) )
colnames(hsm) <- c("pdist","area")
hsm$species <- "HS"
hsm$Source <- "Model"

hp <- cbind.data.frame(c(TRhpH30,TRhpL30),c(rep("High",length(TRhpH30)),rep("Low",length(TRhpL30))) )
colnames(hp) <- c("pdist","area")
hp$species <- "HP"
hp$Source <- "Model"

all <- rbind(gs,hs,hsm,hp)

######################################################
### making a new graph with bins fitting the Figure 3

all$prob_occ_bin <- NA
all$prob_occ_bin[all$pdist<0.02] <- 0
all$prob_occ_bin[all$pdist>=0.02 & all$pdist<0.07] <- 0.05
all$prob_occ_bin[all$pdist>=0.07 & all$pdist<0.15] <- 0.1
all$prob_occ_bin[all$pdist>=0.15 & all$pdist<0.3] <- 0.2
all$prob_occ_bin[all$pdist>=0.3 & all$pdist<0.5] <- 0.4
all$prob_occ_bin[all$pdist>=0.5 & all$pdist<0.7] <- 0.6
all$prob_occ_bin[all$pdist>=0.7 & all$pdist<0.9] <- 0.8
all$prob_occ_bin[all$pdist>=0.9] <- 1

binsss <- c(0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8,1)

all2 <- all

all2 <- all2 %>% group_by(area, species, Source, prob_occ_bin) %>%
  summarise(n = n()) %>%
  mutate(perc = n / sum(n)) 


# 
# ### plotting

# for plotting  have to make sure all values are present for bins and set their values to zero. 
# otherwise, the bars have different width

all2$area <- as.factor(all2$area)
all2$species <- as.factor(all2$species)
all2$Source <- as.factor(all2$Source)
all2$prob_occ_bin <- as.factor(all2$prob_occ_bin)

lvls <- expand.grid(lapply(as.data.frame(all2[, c('area', 'species','Source','prob_occ_bin')]), levels))

# I dont have all levels
tel <- lvls[lvls$species %in% c("GS","HS") & lvls$Source=="Telemetry",]
mod <- lvls[lvls$species %in% c("HP","HS") & lvls$Source=="Model",]
lvls <- rbind(tel,mod)
all2 <- merge(all2, lvls, all=TRUE)

pdist <- ggplot(all2) +
  geom_col(aes(x=as.factor(prob_occ_bin),y=perc*100, fill=area),position = "dodge", width=0.5)+
  ylab("Percentage") + xlab("Probability of exposure")+
  scale_fill_manual(values=c("red","black"))+
  facet_wrap(species~Source, ncol = 2)+
  #scale_x_continuous(limits = c(0,1), breaks = c(0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8,1))+
  labs(fill="Densities")+
  theme(legend.position="bottom")
pdist
