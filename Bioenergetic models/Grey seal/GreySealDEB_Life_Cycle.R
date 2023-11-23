

########## SETUP INITIAL VECTORS & MATRICES ################

max_days <- rep(maximum_days,sim_number)
# years.disturb      <- ceiling(maximum_days/length(Dive.loss.vector)) + 1
# the vector "birthday" will hold the day of the year on which each simulated female was born
birthday <- 0
# These matrices will hold daily values for all the variables for each simulated female
Sa_sim <- GR_sim <- rho_sim <- rho_preg <- Cm_sim <- Ir_sim <- CLact_sim <- surplus_sim <- starv_prop <- matrix(c(0),nrow = maximum_days,ncol = sim_number)
# F_sim is total reserves, F_preg is threshold level of reserves to become pregnant, 
F_sim <- F_preg <- matrix(c(0),nrow = maximum_days,ncol = sim_number)
# preg_state indicates days on which a female can become pregnant
# conceive_record indicates day on which each female actually gives birth
preg_state <- conceive_record <- matrix(c(0),nrow = maximum_days,ncol = sim_number)

# set up matrices for calf and foetus related variables
S_foetus_sim <- calf_demand <- GR_foetus_sim <- matrix(c(0), nrow = maximum_days, ncol = sim_number)
calf_day_record <- calf_deaths <- calf_age <- matrix(c(0), nrow = maximum_days, ncol = sim_number)
Ir_C <- Im_C <- Cm_C <- surplus_C <- rho_C_sim <- F_C <- Sa_C <- GR_C <- matrix(c(0), nrow = maximum_days, ncol = sim_number)

# summary statistics for each simulated female
repro_success <- life_expectancy <- first_calf <- time_as_adult <- rep(c(0),sim_number)
did_female_starve <- starve_death_day <- rep(c(0),sim_number)
ICI <- birth_rate <- calves_born <- calf_starve_deaths <- calf_survival <- calf_survival1 <- rep(c(0),sim_number)

# initialise matrix values to modify intake as a result of disturbance
#########################CHECK THIS###################################
intake_mod <- matrix(c(1),nrow=maximum_days,ncol=sim_number)
disturbance.days <- 0
disturbances <- 0

for (i_sim in 1:sim_number){
  # determine actual birthday of female i_sim
  birthday[i_sim] <- as.integer(rnorm(1,mean_birthday,sd_birthday))
  #day 1 is first day of simulation
  day1 <- birthday[i_sim]
  if (day1 > 365) day1 <- day1 - 365
  julian_days <- c(day1:365)
  remaining <- length(julian_days) + 1
  julian_days[remaining:365] <- c(1:(day1-1))
  # replicate julian_days
  max_years <- ceiling(maximum_days/365) + 1
  julian_days <- rep(julian_days,max_years)
  start_cycle <- which(julian_scale==julian_days[1])
  # create time series of resource density values, with a separate vector for calves so that variation is independent  
    R_cycle <- scalar[start_cycle[1]:(start_cycle[1]+364)]
    R_sim <- R_sim_calf <- Rmean*(rep(R_cycle,max_years))
   
    if (stochastic_R) {
        adult_modifiers <- rbeta(length(R_sim),a_beta,b_beta)/mu_beta
        calf_modifiers <- rbeta(length(R_sim),a_beta,b_beta)/mu_beta
        R_sim <- R_sim*adult_modifiers
        R_sim_calf <- R_sim*calf_modifiers
        }

        if (one_age_class) {life_expectancy[i_sim] <- maximum_age} else {
      # determine life expectancy from age_day1 by sampling from cumulative survival curve
      check <- runif(1,0,cum_Surv[age_day1])
      life_expectancy[i_sim] <- min(which(cum_Surv<=check))
          }
    max_days[i_sim] <- life_expectancy[i_sim] - age_day1
    # max_days[i_sim] <- life_expectancy[i_sim]
    
    # if (is_disturbance){
    #   start1 <- which(julian_days==j.days.to.sample[1])
    #   start.day <- start1[which.birthday]
    #   actual.days <- start.day + (j.days.to.sample - j.days.to.sample[1])
    #     if (model.piling) {
    #     # choose individual history
    #     bad.days <- sapply(p_dist, function(x) rbinom(1,1,x))
    #     intake_mod[days.to.sample[which(bad.days==1)],i_sim]<-disturbance.effect
    #       
    #     # indiv <- sample(c(1:nrow(p_dist)),1,replace=FALSE)
    #     # bad.days <- sapply(p_dist[indiv,], function(x) rbinom(1,1,x))
    #     # intake_mod[actual.days[which(bad.days==1)],i_sim]<-(1-disturbance.effect)
    #     } else {
    #       for (i_dist in 1:years.affected){
    #         bad.days <- sort(sample(days.to.sample,disturbance.number,replace=FALSE))
    #         all.bad.days <- which(julian_days==bad.days[1])
    #         j_bad <- all.bad.days[age.affected + i_dist - 1] - bad.days[1] + bad.days
    #         intake_mod[j_bad,i_sim]<-disturbance.effect}
    #       
    #         # bad.days <- sort(sample(actual.days,disturbance.number,replace=FALSE))
    #         # intake_mod[bad.days,i_sim]<-(1-disturbance.effect)
    #         }
    #     disturbances[i_sim] <- length(which(intake_mod[,i_sim]<1))
    #     ######## modify resource density to take account of disturbance (if any)
    #     R_sim[1:length(intake_mod[,i_sim])] <- R_sim[1:length(intake_mod[,i_sim])]*intake_mod[,i_sim]
    #     # assume calf is disturbed on same days as its mother (might want to modify this)
    #     R_sim_calf[1:length(intake_mod[,i_sim])] <- R_sim_calf[1:length(intake_mod[,i_sim])]*intake_mod[,i_sim]
    #     }
    if (is_disturbance){
      if (model.piling) {
        if (dose.response){
          bad.days <- sapply(p_dist, function(x) rbinom(1,1,x))
          samSize <- length(bad.days[bad.days==1])
          de <- sample(x = distss,
                       size = samSize,
                       replace = TRUE,
                       prob = probss)
          disturbance.effect <- 1 - round(de/24,2) # disturbance effect as a fraction of 24h (rounding is only to satisfy my esthetic needs)
          
          
          intake_mod[days.to.sample[which(bad.days==1)],i_sim]<-disturbance.effect
          
        }
        else
        {
        bad.days <- sapply(p_dist, function(x) rbinom(1,1,x))
        intake_mod[days.to.sample[which(bad.days==1)],i_sim]<-disturbance.effect
        }
        } else {
        for (i_dist in 1:years.affected){
          bad.days <- sort(sample(days.to.sample,disturbance.number,replace=FALSE))
          all.bad.days <- which(julian_days==bad.days[1])
          j_bad <- all.bad.days[age.affected + i_dist - 1] - bad.days[1] + bad.days
          intake_mod[j_bad,i_sim]<-disturbance.effect}
      }
      disturbances[i_sim] <- length(which(intake_mod[,i_sim]<1))
      ######## modify resource density to take account of disturbance (if any)
      
          R_sim[1:length(intake_mod[,i_sim])] <- R_sim[1:length(intake_mod[,i_sim])]*intake_mod[,i_sim]
          # assume calf is disturbed on same days as its mother (might want to modify this)
          R_sim_calf[1:length(intake_mod[,i_sim])] <- R_sim_calf[1:length(intake_mod[,i_sim])]*intake_mod[,i_sim]
      
      # R_sim <- R_sim*intake_mod[,i_sim]
      # # assume calf is disturbed on same days as its mother (might want to modify this)
      # R_sim_calf <- R_sim_calf*intake_mod[,i_sim]
    }
    
    # fill vector with age on each simulated day
    age_days <- rep(c(0),max_days[i_sim])
    age_days[1:(life_expectancy[i_sim]-age_day1+1)] <- c(age_day1:life_expectancy[i_sim])
    alive <- which(age_days>0)

  # transfer values for core mass and growth costs for female from lifetime calculations
  Sa_sim[1,i_sim] <- Sa[age_day1]
  GR_sim[1:length(alive),i_sim] <- GR[age_days[alive]]

  # calculate initial reserve level 
    F_sim[1,i_sim] <- rho_start*Sa_sim[1,i_sim]/(1-rho_start )
    
    # set default value for calf_birth_day
    calf_birth_day <- 0
    
    # loop over simulation days
    #####################################################            
  for (i_day in 1:max_days[i_sim])  {
 
    # all females > min_age are assumed to get pregnant
    if ((age_days[i_day]>min_age) & (julian_days[i_day]==implant_day)) {
    ################ CODE FOR FEMALES THAT ARE ABLE TO BECOME PREGNANT ON i_day ############ 
    # check if ovum implants successfully
    if (rbinom(1,1,fert_success)==1){
      # detemine life expectancy of foetus/calf
    tmp  <- f.foetusandcalf_life(i_sim) 
    conceive <- i_day
    foetal_life <- tmp$life_expect_f
    # calf_days is life expectancy of calf, it was set to 0 in f.foetusandcalf_life if foetus was predicted to die
    calf_days <- tmp$life_expect_c
    conceive_record[conceive,i_sim] <- foetal_life
    # "calf_birth_day" is day on which female gives birth OR day on which foetus dies
    calf_birth_day <- conceive + foetal_life
      # insert foetal daily growth rate values into appropriate elements of GR_foetus_sim
      # but first check if there is enough time for female to give birth before end of simulation
      if (calf_birth_day > max(alive)){
      GR_foetus_sim[conceive:(max(alive)-1),i_sim] <- GR_foetus[1:length(conceive:(max(alive)-1))]
      calf_days <- 0} 
        else {GR_foetus_sim[conceive:(calf_birth_day-1),i_sim] <- GR_foetus[1:foetal_life]}
        }
        # end of calculations for females that implant successfully
      }
      # end of initial calculations relating to pregnancy
      
    #### if it's decision_day, check that there is a foetus
    if ((julian_days[i_day] == decision_day_of_year) & (GR_foetus_sim[i_day,i_sim]>0)){
     # check to see if female mass allows pregnancy to continue
      is_birth <- f.prob_of_pupping((F_sim[i_day,i_sim] + Sa_sim[i_day,i_sim]))
      # check if pregnancy should be terminated because of poor condition
      # if so, zero all future foetal growth and foetal size
      if (is_birth==0)  {
        calf_days <- -1
        GR_foetus_sim[i_day:max_days[i_sim],i_sim] <- S_foetus_sim[i_day:max_days[i_sim],i_sim] <- 0} else {
          # check that foetus will survive until parturition
          if (calf_days > 0) {
            # check that female will survive until end of lactation
            if ((calf_days+calf_birth_day) > max_days[i_sim]) calf_days <- max_days[i_sim] - calf_birth_day
            # modify R_sim to reflect fact that female doesn't feed during lactation and for 4 days before (Pomeroy et al 1999:239)
            end_lact <- ifelse(calf_days<Tl, (calf_birth_day+calf_days),(calf_birth_day+Tl))
            R_sim[(calf_birth_day-4):end_lact] <- 0
            # insert calf ages, and calf growth rate  
            calf_age[calf_birth_day:(calf_birth_day+calf_days-1),i_sim] <- c(1:calf_days)
            GR_C[calf_birth_day:(calf_birth_day+calf_days-1),i_sim] <- GR[1:calf_days]}
            } 
      
        # calf_day_record = -1 if it's aborted because of poor condition
        calf_day_record[calf_birth_day, i_sim] <- calf_days
          }
          # end of decision day calculations
  
    ###### calculate energy intake of mother
    # check if mother is still alive
    if(age_days[i_day]>0){
      rho_sim[i_day,i_sim] <- (F_sim[i_day,i_sim])/(F_sim[i_day,i_sim] + Sa_sim[i_day,i_sim])
      if (rho_sim[i_day,i_sim] <= 0) rho_sim[i_day,i_sim] <- 0.000001
      Ir_sim[i_day,i_sim] <- assim_effic[age_days[i_day]]*R_sim[i_day]*Sa_sim[i_day,i_sim]^(2/3)/(1 + exp(-eta*(rho/rho_sim[i_day,i_sim] - 1)))
      ############## mass of foetus included in cost of maintenance ##########
      Cm_sim[i_day,i_sim] <- Sigma_M*(Sa_sim[i_day,i_sim]+ S_foetus_sim[i_day,i_sim] + Theta_F*F_sim[i_day,i_sim])^0.75
      starv_prop[i_day,i_sim] <- (1 - xi_m)*(rho_sim[i_day,i_sim] - rho_s)/((rho - rho_s) - xi_m*(rho_sim[i_day,i_sim] - rho_s))
      # female abandons lactation if rho_sim is close to the starvation threshold and can start feeding at a reduced rate
      if (rho_sim[i_day,i_sim] < (rho_s+0.05)) 
      {starv_prop[i_day,i_sim] <- 0
      if (R_sim[i_day] <- 0) R_sim[i_day] <- Rmean/2}
    }
    
  ###### calculate energy intake of current calf
    if (calf_age[i_day,i_sim] > 0){
      # if calf_age = 1, set total mass of calf to that of foetus and its reserves to rho_s
      if (calf_age[i_day,i_sim]==1) {
        F_C[i_day,i_sim] <- rho_s*S_foetus_sim[i_day,i_sim]
         Sa_C[i_day,i_sim] <- S_foetus_sim[i_day,i_sim] - F_C[i_day,i_sim]
         }
      
      rho_C_sim[i_day,i_sim] <- F_C[i_day,i_sim]/(F_C[i_day,i_sim] + Sa_C[i_day,i_sim])
      if (rho_C_sim[i_day,i_sim] <= 0) rho_C_sim[i_day,i_sim] <- 0.000001
      # check for starvation related death of calf
      if (rho_C_sim[i_day,i_sim] < rho_s) {
          is_calf_dead <- f.starve_death(rho_C_sim[i_day,i_sim],rho_s)
          if (is_calf_dead == 0) 
            {calf_deaths[i_day,i_sim] <- calf_age[i_day,i_sim]
              calf_age[(i_day+1):max_days[i_sim],i_sim] <- 0}
              }
      # determine demand for milk
      # check if mother is still alive
      if(age_days[i_day]>0){
      calf_demand[i_day,i_sim] <- 1/(1 + exp(-eta*(rho_C/rho_C_sim[i_day,i_sim] - 1)))
      calf_demand[i_day,i_sim] <- ifelse(calf_demand[i_day,i_sim]>1,1,calf_demand[i_day,i_sim])
      # determine milk intake
      Im_C[i_day,i_sim] <- phi_L*calf_demand[i_day,i_sim]*milk_prop[calf_age[i_day,i_sim]]*starv_prop[i_day,i_sim]*Sa_C[i_day,i_sim]^(2/3)
      }
      # determine calf's food intake
      Ir_C[i_day,i_sim] <- R_sim_calf[i_day]*assim_effic[calf_age[i_day,i_sim]]*Sa_C[i_day,i_sim]^(2/3)/(1 + exp(-eta*(rho/rho_C_sim[i_day,i_sim] - 1)))
      # determine cost of maintenance for calf
      Cm_C[i_day,i_sim] <- Sigma_M*(Sa_C[i_day,i_sim] + F_C[i_day,i_sim]*Theta_F)^0.75
      
      
      # determine calf's surplus and modify reserves accordingly
      I_C_total <- Ir_C[i_day,i_sim]+Im_C[i_day,i_sim]
      # is there enough energy to cover maintenance and growth?
      if (I_C_total < (Cm_C[i_day,i_sim] + Sigma_G*GR_C[i_day,i_sim])){
        # check if theres's any surplus. If not, set growth to zero
        if (I_C_total > 0){
        ####### apply Kappa rule and modify growth rates to fit if there's not enough energy
          available_for_calf_growth <- ifelse((I_C_total-Cm_C[i_day,i_sim])>I_C_total*(1-Kappa),(I_C_total-Cm_C[i_day,i_sim]),I_C_total*(1-Kappa))
          if(available_for_calf_growth < Sigma_G*GR_C[i_day,i_sim])
          GR_C[i_day,i_sim] <- available_for_calf_growth/Sigma_G 
            } else {GR_C[i_day,i_sim] <- 0}
          }
        surplus_C[i_day,i_sim] <- I_C_total - Cm_C[i_day,i_sim] - Sigma_G*GR_C[i_day,i_sim]
        energy.density_C <- ifelse(surplus_C[i_day,i_sim] > 0, epsi_plus_pups, epsi_minus)
        # modify calf's energy reserves
        F_C[(i_day+1),i_sim] <- F_C[i_day,i_sim] + surplus_C[i_day,i_sim]/energy.density_C
        # increase structural mass of calf
        Sa_C[(i_day+1),i_sim] <- Sa_C[i_day,i_sim] + GR_C[i_day,i_sim]
        # end of calf calculation
        #########################
        }
    
    # check if mother is still alive
    if(age_days[i_day]>0){
      # calculate cost of lactation (if any)
      CLact_sim[i_day, i_sim] <- Im_C[i_day,i_sim]/Sigma_L
      # determine mother's surplus and modify her reserves accordingly
      if (Ir_sim[i_day,i_sim] < (Cm_sim[i_day,i_sim] + CLact_sim[i_day,i_sim] + Sigma_G*(GR_sim[i_day,i_sim] + GR_foetus_sim[i_day,i_sim]))){
        surplus2 <- Ir_sim[i_day,i_sim] - CLact_sim[i_day,i_sim]
        # is there any spare energy once the costs of lactation have been subtracted?
          if (surplus2 > 0){
          ####### apply Kappa rule and modify growth rates to use available energy
          available_for_growth <- ifelse((surplus2-Cm_sim[i_day,i_sim])>surplus2*(1-Kappa),(surplus2-Cm_sim[i_day,i_sim]),surplus2*(1-Kappa))
            if (available_for_growth < Sigma_G*(GR_sim[i_day,i_sim] + GR_foetus_sim[i_day,i_sim])){
            growth_ratio <- GR_sim[i_day,i_sim]/(GR_sim[i_day,i_sim] + GR_foetus_sim[i_day,i_sim])
            GR_sim[i_day,i_sim] <- growth_ratio*available_for_growth/Sigma_G 
            GR_foetus_sim[i_day,i_sim] <- (1-growth_ratio)*available_for_growth/Sigma_G 
            }
          } else {GR_sim[i_day,i_sim]<- GR_foetus_sim[i_day,i_sim] <- 0}
      } 
    # end of modifications to maternal and foetal growth rates
    ###########################
      surplus_sim[i_day] <- Ir_sim[i_day,i_sim] - Cm_sim[i_day,i_sim] - CLact_sim[i_day,i_sim] - Sigma_G*(GR_sim[i_day,i_sim] + GR_foetus_sim[i_day,i_sim])
    # calculate female reserve mass in next time step
    energy.density <- ifelse(surplus_sim[i_day] > 0, epsi_plus, epsi_minus)
    F_sim[(i_day+1),i_sim] <- F_sim[i_day,i_sim] + surplus_sim[i_day]/energy.density
    if (F_sim[(i_day+1),i_sim] < 0) F_sim[(i_day+1),i_sim] <- 0
    # calculate structural mass of female and her foetus (if there is one) in next time step
    Sa_sim[(i_day+1),i_sim] <- Sa_sim[i_day,i_sim] + GR_sim[i_day,i_sim]
    # zero foetus size after its born or it dies
    if (i_day==calf_birth_day) S_foetus_sim[(i_day):max_days[i_sim],i_sim] <- 0
    S_foetus_sim[(i_day+1),i_sim] <- S_foetus_sim[i_day,i_sim] + GR_foetus_sim[i_day,i_sim]
    }
    # end of simulated day
    }

######## check if female has starved to death ###############
  # have female reserves ever dropped below rho_s?
  female_starve_days <- which(rho_sim[1:max(alive),i_sim]<rho_s)
  did_female_starve[i_sim] <- length(female_starve_days)
    if (length(female_starve_days>0)){
    # is_female_dead will contain days on which female is predicted to die
    is_female_dead <- 1
    for (i_starve in 1:length(female_starve_days)){
      # carry out Bernouilli trial on each female_starve_day to see if female dies
      is_female_dead[i_starve] <- f.starve_death(rho_sim[female_starve_days[i_starve],i_sim],rho_s)} 
      death_days <- which(is_female_dead==0)
        if (length(death_days)>0) {
          # death_day is earliest day on which death is predicted to occur
          death_day <- female_starve_days[min(death_days)]
          # remove all foetuses that were predicted to be conceived after death_day
          conceive_record[death_day:max_days[i_sim],i_sim] <- 0
          # remove all calves that were predicted to be born after death_day
          calf_age[death_day:max_days[i_sim],i_sim] <- rep(c(0),length(death_day:max_days[i_sim]))
          calf_deaths[death_day:max_days[i_sim],i_sim] <- 0
          # record day on which starvation death occurred and modify life_expectancy
          starve_death_day[i_sim] <- death_day
          life_expectancy[i_sim] <- death_day + age_day1
          } 
        } # end of female starvation calculations
      
  # calculate summary statistics
    life_expectancy[i_sim] <- life_expectancy[i_sim]/365
    time_as_adult[i_sim] <- length(which(conceive_record[,i_sim]>0))
    repro_success[i_sim] <- length(which(calf_age[,i_sim]==max_age_calf))/2
    calves_born[i_sim] <- length(which(calf_age[,i_sim]==1))
    if (time_as_adult[i_sim] > 0) {
      # check to see if last foetus can be born before female dies
      if((conceive_record[time_as_adult[i_sim],i_sim]+Tp)>max(alive)) time_as_adult[i_sim] <- time_as_adult[i_sim]-1
      birth_rate[i_sim] <- calves_born[i_sim]/time_as_adult[i_sim]}
      calf_starve_deaths[i_sim] <- length(which(calf_deaths[,i_sim]>0))
      calf_birth_days <- which(calf_age[,i_sim]==1)
      if (length(calf_birth_days)>1) {
          birth_intervals <- sapply(2:length(calf_birth_days), function(X) calf_birth_days[X]-calf_birth_days[X-1])
          ICI[i_sim] <- mean(birth_intervals/365)
          }
    }
# end of loop over multiple females


  
