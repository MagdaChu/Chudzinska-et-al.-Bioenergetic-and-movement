rm(list = ls())
# input species
library(parallel)
library(doParallel)
library(foreach)

# set your working directory where you store in the input code
# NOTE THAT YOU HAVE TO REDEFINE YOUR DIRECTORY AROUND LINE 80 AGAIN
setwd("")

# define number of computer cores you want to use. The code is written to run parallel simulations

numCores <- detectCores()-2 # default is all computer cores - 2

# for full combination of p_dist and disturbance.effect
qdisturbance.effect.h <- c(0,1,2,4,6,8,12) # disturbance effect in hours
p_dist1 <- c(0.05,0.1,0.2,0.4,0.6,0.8,1)
comb <- expand.grid(disturbance.effect.h = disturbance.effect.h ,p_dist1 = p_dist1)
trr <- comb[comb$disturbance.effect.h==0,]

# removing all combinations with disturbance = 0h. Only one is needed
comb <- comb[c(-8,-15,-22,-29,-36,-43),]

# upload piling schedule. You can upload your own, just make sure you follow the format as in the example
piling.file <- "DEB_PS1.csv" 
pile <- read.csv(file = piling.file, header = TRUE)
All_Dates <- strptime(pile[,3],format="%d/%m/%Y")
All_Durations <- pile[,2]
Dates <- unique(All_Dates)

load("Porpoise_ParamsFromABC_superNonstrict.RData")

# determine Julian date for first day in piling file
first.day <- Dates$yday[1] + 1
age.affected <- c(17)

source('PorpoiseDEB_Params.R')
first.julian_days <- which(julian_days==first.day)
start_pile <- first.julian_days[age.affected] # first day of piling at the age when porpoise is affected
difStartP2B <- first.day + ((365-birthday) + 10)
start_adcalf_birthdeath_new4s <- start_pile - difStartP2B
###

sim_number_full <- 2002 # define your modelled population size

st <- Sys.time()
params <- list()
finalResults100 <- list()
finalResultsDisturb <- list()
# to output the actual time of death for adults and babies
finalResults100DeathListAd <- list()
finalResultsDisturbDeathListAd <- list()
finalResults100DeathListPup <- list()
finalResultsDisturbDeathListPup <- list()

for ( dd in 1:nrow(comb)){ #nrow(comb) 

    for (difparam in 1:100){ # in this loop, each disturbance scenario is run 100 times, each time drawing from different, plausible parameter combination

doParallel::registerDoParallel(cl <- makeCluster(spec = numCores, type = "PSOCK",setup_strategy = "sequential"))
  
  #set.seed(difparam)     
  gamma <- paramFilt$gamma[difparam]
  Tr <- paramFilt$Tr[difparam]
  mu_s <- paramFilt$mu_s[difparam]
  Sigma_M_full <- paramFilt$Sigma_M_full[difparam]
  Rmean <- paramFilt$Rmean[difparam]
  Tn <- paramFilt$Tn[difparam]
  
  
  #set.seed(Sys.time())
  
  clusterExport(cl, list("gamma","Tr","mu_s","Sigma_M_full","Rmean", "numCores","Tn","comb","dd","difparam","start_pile","age.affected","pile","sim_number_full","start_adcalf_birthdeath_new4s")) #, envir = env "eta", "xi_m", "skipping_point",
  
results <- parallel::clusterEvalQ(cl = cl, expr = { 
  setwd("")
  

  spec <- 'HP'

  
  # # specify no disturbance
  is_disturbance <- TRUE
  model.piling <- TRUE
  dose.response <- FALSE

  # all files and parameters related to piling

  All_Dates <- strptime(pile[,3],format="%d/%m/%Y")
  All_Durations <- pile[,2]

  # combine duration of all piling activities that start on the same day (should not be the case if we model monopiles = one pile per day)
  Dates <- unique(All_Dates)
  Durations <- 0
  for (i in 1:length(Dates)){Durations[i] <- sum(All_Durations[which(All_Dates==Dates[i])])}
  years.affected <- length(unique(Dates$year))

  # determine Julian date for first day in piling file
  first.day <- Dates$yday[1] + 1
  # calculate number of days since first piling event to use as index
  days.to.sample1 <- as.integer(round(as.numeric(difftime(Dates[2:length(Dates)],Dates[1]))))
  
  disturbance.effect <- round(comb$disturbance.effect.h[dd]/24,2) # disturbance effect as a fraction of 24h (rounding is only to satisfy my esthetic needs)
  disturbance.effect <- 1 - disturbance.effect 
  
  p_dist <- rep(comb$p_dist1[dd],length(Durations))

  #for dose response
  distss <- c(12,8,6,4,2,1)
  probss <- c(0.025, 0.05, 0.075, 0.13, 0.25, 0.47)
  
  # set one_age_class to TRUE if all simulated animals are to survive to the maximum age
  one_age_class <- TRUE
  #age.affected <- c(7)
  maximum_age <- (age.affected + 2)*365
 
  # number of females to be simulated
  sim_number_temp <- sim_number_full
  sim_number <- sim_number_full/numCores
  # should R vary from day to day?
  stochastic_R <- FALSE
  
  if (one_age_class) max_age <- age.affected + years.affected + 10
  if (max_age > 30) max_age <- 30

  
 source('PorpoiseDEB_Params_abc.R')
  

  
  #### from dublin array porpoise code
  
  first.julian_days <- which(julian_days==first.day)
  start_pile <- first.julian_days[age.affected] # first day of piling at the age when porpoise is affected
  days.to.sample <- start_pile
  days.to.sample[2:length(Dates)] <- start_pile + days.to.sample1 # actual piling duration
  start1 <- which(julian_days==1)
  
  # from now calculations for pup/calf death and birth, adult survival and conception are as follow
  
  ### spanning over 4 annual cycles ###
  # pup/calf birth and deaths and adult deaths should start 
  # at first birthing season before the piling starts 
  # and continue after the first birthing season after piling ends
  # to just before the next next birthing season starts
  
  difStartP2B <- first.day + ((365-birthday) + 10)
  start_adcalf_birthdeath_new4s <- start_pile - difStartP2B
  # end_calf_new <- start_pile + (3 * 365) # so calculations end one year after the end of last piling day
  difStopP2B <- (birthday - 10)- first.day
  end_adcalf_birthdeath_new4s <- start_pile + (years.affected * 365) + difStopP2B# so calculations end one year after the end of last piling day
  #(end_adcalf_birthdeath_new4s-start_adcalf_birthdeath_new4s)/365
  
  # conception should start from first conception before piling to 
  # first conception after piling
  
  difStartP2C <- (365-implant_start+10) + first.day
  start_calf_conc_new4s <- start_pile - difStartP2C
  # end_calf_new <- start_pile + (3 * 365) # so calculations end one year after the end of last piling day
  difStopP2C <- implant_start-10 - first.day
  end_calf_conc_new4s <- start_pile + (years.affected * 365) + difStopP2C
  #(end_calf_conc_new4s-start_calf_conc_new4s)/365
  
  source('PorpoiseDEB_Life_Cycle.R')
  
  #calf_deaths_during_piling3s <- calf_births_during_piling3s  <- conceptions_during_piling3s <- ad_deaths_starvation_during_piling3s <- rep(c(0),sim_number)
  calf_deaths_during_piling4s <- calf_births_during_piling4s  <- conceptions_during_piling4s <- ad_deaths_starvation_during_piling4s <- when_born <- when_died_calves <- when_died_adults <- when_conceived <- rep(c(0),sim_number)
  when_born <- when_died_calves <- when_died_adults <- when_conceived <- list()
  
  for (L3 in 1:sim_number) {
    
    # calculations over 4 annual cycles
    calf_deaths_during_piling4s[L3] <- length(which(calf_deaths[start_adcalf_birthdeath_new4s:end_adcalf_birthdeath_new4s,L3]>0))
    calf_births_during_piling4s[L3]<- length(which(calf_age[start_adcalf_birthdeath_new4s:end_adcalf_birthdeath_new4s,L3]==1))
    conceptions_during_piling4s[L3] <- length(which(conceive_record[start_calf_conc_new4s:end_calf_conc_new4s,L3]>0))
    #calf_survivors[L3] <- length(which(calf_age[start_calf_new:end_calf_new,i]==max_age_calf))
    if (starve_death_day[[L3]]<end_adcalf_birthdeath_new4s & starve_death_day[L3] >= start_adcalf_birthdeath_new4s)  ad_deaths_starvation_during_piling4s[L3] <- ad_deaths_starvation_during_piling4s[L3] + 1
    
    
    when_died_calves[[L3]] <- which(calf_deaths[start_adcalf_birthdeath_new4s:end_adcalf_birthdeath_new4s,L3]>0)
    when_born[[L3]] <- which(calf_age[start_adcalf_birthdeath_new4s:end_adcalf_birthdeath_new4s,L3]==1)
    when_conceived[[L3]] <- which(conceive_record[start_calf_conc_new4s:end_calf_conc_new4s,L3]>0)  
    
    
  }
  
  when_died_calves <- unlist(when_died_calves)
  if (length(when_died_calves)==0) when_died_calves <- 0
  when_born <- unlist(when_born)
  if (length(when_born)==0) when_born <- 0
  when_conceived <- unlist(when_conceived)
  if (length(when_conceived)==0) when_conceived <- 0
  
  when_died_adults <- starve_death_day
  if (length(when_died_adults)==0) when_died_adults <- 0 # this should not happen as they all die when reaching max age
  
  list(
    calf_deaths_during_piling4s,
    calf_births_during_piling4s,
    conceptions_during_piling4s,
    ad_deaths_starvation_during_piling4s,
    when_died_calves,
    when_born,
    when_conceived,
    when_died_adults
  )
  


})
parallel::stopCluster(cl) 

## merging results as now each element of results

calf_deaths4s <- c(results[[1]][[1]])
calf_births4s <- c(results[[1]][[2]])
conceptions4s <- c(results[[1]][[3]])
ad_deaths4s <- c(results[[1]][[4]])

when_died_calves <- c(results[[1]][[5]])
when_born <- c(results[[1]][[6]])
when_conceived <- c(results[[1]][[7]])
when_died_adults <- c(results[[1]][[8]])


for (nn in 2:numCores){
  calf_deaths4s <- c(calf_deaths4s,results[[nn]][[1]])
  calf_births4s <- c(calf_births4s,results[[nn]][[2]])
  conceptions4s <- c(conceptions4s,results[[nn]][[3]])
  ad_deaths4s <- c(ad_deaths4s,results[[nn]][[4]])
  
  when_died_calves <- c(when_died_calves,results[[nn]][[5]])
  when_born <- c(when_born,results[[nn]][[6]])
  when_conceived <- c(when_conceived,results[[nn]][[7]])
  when_died_adults <- c(when_died_adults,results[[nn]][[8]])
  
}

calf_survivals4s <- sum(calf_deaths4s)/sum(calf_births4s)
nadults <- length(when_died_adults)

birth_rates4s <- sum(calf_births4s)/nadults
adult_mortality4s <- sum(ad_deaths4s)/nadults
fertility4s<- sum(conceptions4s)/nadults
fertility4s2<- sum(conceptions4s)/sum(calf_births4s)


endYear1bd <- 365#mean_birthday- 10 - first.day
endYear2bd <- endYear1bd + 365
endYear3bd <- endYear1bd + (2*365)
endYear4bd <- endYear1bd + (3*365)
endYear5bd <- endYear1bd + (4*365)

# conception
endYear1c <- 365#implant_day - 10 + (365 - first.day)
endYear2c <- endYear1c+365
endYear3c <- endYear1c+(2*365)
endYear4c <- endYear1c+(3*365)
endYear5c <- endYear1c+(4*365)

died <- when_died_calves
died[died>0 & died<endYear1bd] <- 1
died[died>=endYear1bd & died<endYear2bd] <- 2
died[died>=endYear2bd & died<endYear3bd] <- 3
died[died>=endYear3bd & died<endYear4bd] <- 4
died[died>=endYear4bd & died<endYear5bd] <- 5
died <- factor(died,levels=c(1:5))
died <- as.data.frame(table(died))
colnames(died) <- c("Year","Freq")

born <- when_born
born[born>0 & born<endYear1bd] <- 1
born[born>=endYear1bd & born<endYear2bd] <- 2
born[born>=endYear2bd & born<endYear3bd] <- 3
born[born>=endYear3bd & born<endYear4bd] <- 4
born[born>=endYear4bd & born<endYear5bd] <- 5
born <- factor(born,levels=c(1:5))
born <- as.data.frame(table(born))
colnames(born) <- c("Year","Freq")

conc <- when_conceived
conc[conc>0 & conc<endYear1c] <- 1
conc[conc>=endYear1c & conc<endYear2c] <- 2
conc[conc>=endYear2c & conc<endYear3c] <- 3
conc[conc>=endYear3c & conc<endYear4c] <- 4
conc[conc>=endYear4c & conc<endYear5c] <- 5
conc <- factor(conc,levels=c(1:5))
conc <- as.data.frame(table(conc))
colnames(conc) <- c("Year","Freq")

# for adult deaths, it is real time in the model so start of the calculation is as start_adcalf_birthdeath_new4s

endYear1ad <- start_adcalf_birthdeath_new4s+365
endYear2ad <- endYear1ad + 365
endYear3ad <- endYear1ad + (2*365)
endYear4ad <- endYear1ad + (3*365)
endYear5ad <- endYear1ad + (4*365)

diedAd <- when_died_adults
diedAd <- diedAd[diedAd!=0]
diedAd[diedAd>0 & diedAd<endYear1ad] <- 1
diedAd[diedAd>=endYear1ad & diedAd<endYear2ad] <- 2
diedAd[diedAd>=endYear2ad & diedAd<endYear3ad] <- 3
diedAd[diedAd>=endYear3ad & diedAd<endYear4ad] <- 4
diedAd[diedAd>=endYear4ad & diedAd<endYear5ad] <- 5
diedAd <- factor(diedAd,levels=c(1:5))
diedAd <- as.data.frame(table(diedAd))
colnames(diedAd) <- c("Year","Freq")

when <- as.data.frame(matrix(c(died$Freq,born$Freq,conc$Freq,diedAd$Freq), ncol=20))
colnames(when) <- c("DiedY1","DiedY2","DiedY3","DiedY4","DiedY5","BornY1","BornY2","BornY3","BornY4","BornY5","ConcY1","ConcY2","ConcY3","ConcY4","ConcY5","DiedAdY1","DiedAdY2","DiedAdY3","DiedAdY4","DiedAdY5")
finalResults100[[difparam]] <- cbind.data.frame(calf_survivals4s,birth_rates4s,adult_mortality4s,fertility4s,fertility4s2, when,difparam,dd,comb[dd,]) #calf_survivalsAsJohn, calf_survivals3s,birth_rates3s,adult_mortality3s,fertility3s,
finalResults100DeathListAd[[difparam]] <- when_died_adults
finalResults100DeathListPup[[difparam]] <- when_died_calves
print(difparam)

}
  finalResultsDisturb[[dd]] <- do.call("rbind",finalResults100)
  finalResultsDisturbDeathListAd[[dd]] <- finalResults100DeathListAd
  finalResultsDisturbDeathListPup[[dd]] <- finalResults100DeathListPup
  print(dd)
} # end of big disturbance loop

stp <- Sys.time()
stp-st # should be around 38 hours

save(finalResultsDisturb, file="Porpoise_S9_NewPeriod4vitalratesCalculations.RData")
save(finalResultsDisturbDeathListAd, file="Porpoise_S9_NewPeriod4vitalratesCalculationsDeathListA.RData")
save(finalResultsDisturbDeathListPup, file="Porpoise_S9_NewPeriod4vitalratesCalculationsDeathListPup.RData")

