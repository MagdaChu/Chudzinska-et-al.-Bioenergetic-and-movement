

########## SETUP INITIAL VECTORS & MATRICES ################

max_days <- rep(maximum_days,sim_number)
# years.disturb      <- ceiling(maximum_days/length(Dive.loss.vector)) + 1

# These matrices will hold daily values for all the variables for each simulated female
Sa_sim <- rho_sim <- rho_preg <- Cm_sim <- Ir_sim <- CGest_sim <- surplus_sim <- starv_prop <- matrix(c(0),nrow = maximum_days,ncol = sim_number)
# F_sim is total reserves, F_preg is threshold level of reserves to become pregnant, 
F_sim <- F_preg <- matrix(c(0),nrow = maximum_days,ncol = sim_number)
# preg_state indicates days on which a female can become pregnant
# conceive_record indicates day on which each female actually gives birth
preg_state <- conceive_record <- matrix(c(0),nrow = maximum_days,ncol = sim_number)

# set up matrices for calf related variables
calf_demand <- CGa_C <- matrix(c(0), nrow = maximum_days, ncol = sim_number)
calf_day_record <- calf_deaths <- calf_age <- matrix(c(0), nrow = maximum_days, ncol = sim_number)
Ir_C <- Im_C <- Cm_C <- surplus_C <- rho_C <- F_C <- Sa_C <- matrix(c(0), nrow = maximum_days, ncol = sim_number)

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
   
    R_sim <- rep(Rmean,maximum_days)
    if (stochastic_R) {R_sim <- Rmean*rbeta(maximum_days,a,b)/mu}
    if (one_age_class) {life_expectancy[i_sim] <- maximum_age} else {
      # determine life expectancy from age_day1 by sampling from cumulative survival curve
      check <- runif(1,0,cum_Surv[age_day1])
      life_expectancy[i_sim] <- min(which(cum_Surv<=check))
          }
      max_days[i_sim] <- life_expectancy[i_sim] - age_day1
    
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
        R_sim <- R_sim*intake_mod[,i_sim]
        }

    
    # fill vector with age on each simulated day
    age_days <- c(age_day1:(age_day1+max_days[i_sim]-1))

  # transfer values for core mass and growth costs for female from lifetime calculations
  Sa_sim[1:max_days[i_sim],i_sim] <- Sa[age_days]
  CGa_sim <- CGa[age_days]

  # calculate initial reserve level 
    F_sim[1,i_sim] <- rho_j*Sa_sim[1,i_sim]/(1-rho_j)
    # identify days on which female can become pregnant
    preg_state[1:max_days[i_sim],i_sim] <- preg_possible[1:max_days[i_sim]]
    
    # loop over simulation days
    #####################################################            
  for (i_day in 1:max_days[i_sim])  {
 
  F_preg[i_day,i_sim] <- F_neo + rho_s[i_day]*Sa_sim[i_day,i_sim]/(1-rho_s[i_day])
  rho_preg[i_day,i_sim] <- F_preg[i_day,i_sim]/(F_preg[i_day,i_sim] + Sa_sim[i_day,i_sim])
  # check if female can become pregnant
  if (preg_state[i_day,i_sim] > 0 & F_sim[i_day,i_sim] >= F_preg[i_day,i_sim]){
    ################ CODE FOR FEMALES THAT ARE ABLE TO BECOME PREGNANT ON i_day ############ 
    # reset pregnancy status to 0 for rest of fertile period
    preg_state[(i_day+1):(i_day+implant_period),i_sim] <- 0
    # check if ovum implants successfully
    if (rbinom(1,1,fert_success)==1){
    tmp  <- f.foetusandcalf_life(i_sim) 
    conceive <- i_day
    foetal_life <- tmp$life_expect_f
    conceive_record[conceive,i_sim] <- foetal_life
    # "birth_day" is day on which female gives birth OR day on which foetus dies
    birth_day <- conceive + foetal_life
      # check if there is enough time for female to give birth before end of simulation
      if (birth_day > max_days[i_sim]){
      CGest_sim[conceive:max_days[i_sim],i_sim] <- CGest[1:length(conceive:max_days[i_sim])]
      } else {
      CGest_sim[conceive:(birth_day-1),i_sim] <- CGest[1:foetal_life]
      # calf_days is life expectancy of calf (max value = 345 days)
      calf_days <- tmp$life_expect_c
      calf_day_record[birth_day, i_sim] <- calf_days
      # check if foetus dies during pregnancy (in which case calf_days was set to 0)
      if (calf_days>0){
        if ((calf_days+birth_day) > max_days[i_sim]) calf_days <- max_days[i_sim] - birth_day
        # insert calf ages, calf core mass and calf growth costs  
        calf_age[birth_day:(birth_day+calf_days-1),i_sim] <- c(1:calf_days)
        Sa_C[birth_day:(birth_day+calf_days-1),i_sim] <- Sa[1:calf_days]
        CGa_C[birth_day:(birth_day+calf_days-1),i_sim] <- CGa[1:calf_days]
        }
        # end of calculations for foetuses that are actually born
        ############################################################### 
          }
          # end of calculations for females that will give birth after simulation ends 
            }
            # end of calculations for females that fail to implant successfully
    }
    # end of all calculations relating to females who can become pregnant on i_day
  
  # calculate energy intake of mother
  rho_sim[i_day,i_sim] <- (F_sim[i_day,i_sim])/(F_sim[i_day,i_sim] + Sa_sim[i_day,i_sim])
  if (rho_sim[i_day,i_sim] <= 0) rho_sim[i_day,i_sim] <- 0.000001
  Ir_sim[i_day,i_sim] <- assim_effic[age_days[i_day]]*R_sim[i_day]*Sa_sim[i_day,i_sim]^(2/3)/(1 + exp(-eta*(rho_t[i_day]/rho_sim[i_day,i_sim] - 1)))
  Cm_sim[i_day,i_sim] <- Sigma_M*(Sa_sim[i_day,i_sim]+Theta_F*F_sim[i_day,i_sim])^0.75
  starv_prop[i_day,i_sim] <- (1 - xi_m)*(rho_sim[i_day,i_sim] - rho_s[i_day])/((rho_t[i_day] - rho_s[i_day]) - xi_m*(rho_sim[i_day,i_sim] - rho_s[i_day]))
  if (rho_sim[i_day,i_sim] < rho_s[i_day]) starv_prop[i_day,i_sim] <- 0
    
  ###### calculate milk intake of current calf
    if (calf_age[i_day,i_sim] > 0){
      # set calf's energy reserves to initial level if calf_age = 1
      if (calf_age[i_day,i_sim]==1) F_C[i_day,i_sim] <- F_start
      rho_C[i_day,i_sim] <- F_C[i_day,i_sim]/(F_C[i_day,i_sim] + Sa_C[i_day,i_sim])
      if (rho_C[i_day,i_sim] <= 0) rho_C[i_day,i_sim] <- 0.000001
      # check for starvation related death of calf
      if (rho_C[i_day,i_sim] < rho_s[i_day]) {
          is_calf_dead <- f.calf_starve_death(rho_C[i_day,i_sim],rho_s[i_day])
          if (is_calf_dead == 0) 
            {calf_deaths[i_day,i_sim] <- calf_age[i_day,i_sim]
              calf_age[(i_day+1):max_days[i_sim],i_sim] <- 0
             }
          }
      # determine demand for milk
      calf_demand[i_day,i_sim] <- 1/(1 + exp(-eta*(rho_t[i_day]/rho_C[i_day,i_sim] - 1)))
      calf_demand[i_day,i_sim] <- ifelse(calf_demand[i_day,i_sim]>1,1,calf_demand[i_day,i_sim])
      # determine milk intake
      Im_C[i_day,i_sim] <- phi_L*calf_demand[i_day,i_sim]*milk_prop[calf_age[i_day,i_sim]]*starv_prop[i_day,i_sim]*Sa_C[i_day,i_sim]^(2/3)
      
        # determine calf's food intake
        Ir_C[i_day,i_sim] <- R_sim[i_day]*assim_effic[calf_age[i_day,i_sim]]*Sa_C[i_day,i_sim]^(2/3)/(1 + exp(-eta*(rho_t[i_day]/rho_C[i_day,i_sim] - 1)))
        # determine cost of maintenance. 
        Cm_C[i_day,i_sim] <- Sigma_M*(Sa_C[i_day,i_sim] + F_C[i_day,i_sim]*Theta_F)^0.75
        # determine calf's surplus
        surplus_C[i_day,i_sim] <- Im_C[i_day,i_sim] + Ir_C[i_day,i_sim] - Cm_C[i_day,i_sim] - CGa_C[i_day,i_sim]
        efficC <- ifelse(surplus_C[i_day,i_sim] > 0, epsi_plus[i_day], epsi_minus[i_day])
        # modify calf's energy reserves
        F_C[(i_day+1),i_sim] <- F_C[i_day,i_sim] + surplus_C[i_day,i_sim]/efficC
        # end of lactation calculation
        }
    
    # determine mother's surplus and modify her reserves accordingly
    # calculate cost of lactation (if any)
    CLact_sim <- Im_C[i_day,i_sim]/Sigma_L
    surplus_sim[i_day]  <- Ir_sim[i_day,i_sim] - Cm_sim[i_day,i_sim] -  CLact_sim - CGest_sim[i_day,i_sim] - CGa_sim[i_day]
    effic <- ifelse(surplus_sim[i_day] > 0, epsi_plus[i_day], epsi_minus[i_day])
    F_sim[(i_day+1),i_sim] <- F_sim[i_day,i_sim] + surplus_sim[i_day]/effic
    if (i_day > max_days[i_sim]) F_sim[(i_day+1),i_sim] <- 0
    if (F_sim[(i_day+1),i_sim] < 0) F_sim[(i_day+1),i_sim] <- 0
  # end of simulated day
    }

######## check if female has starved to death ###############
  # have female reserves ever dropped below rho_s?
  female_starve_days <- which(rho_sim[1:max_days[i_sim],i_sim]<rho_s[1:max_days[i_sim]])
  did_female_starve[i_sim] <- length(female_starve_days>0)
    if (length(female_starve_days>0)){
    # is_female_dead will contain days on which female is predicted to die
    is_female_dead <- 1
    for (i_starve in 1:length(female_starve_days)){
      # carry out Bernouilli trial on each female_starve_day to see if female dies
      is_female_dead[i_starve] <- f.calf_starve_death(rho_sim[female_starve_days[i_starve],i_sim],rho_s[female_starve_days[i_starve]])} 
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
        }
      
  
  # calculate summary statistics
    time_as_adult[i_sim] <- (max_days[i_sim] - (min_age+365))
    repro_success[i_sim] <- length(which(calf_age[,i_sim]==Tl))/2
    calves_born[i_sim] <- length(which(calf_age[,i_sim]==1))
    if (time_as_adult[i_sim] > 0) birth_rate[i_sim] <- calves_born[i_sim]/ceiling(time_as_adult[i_sim]/365)
      calf_starve_deaths[i_sim] <- length(which(calf_deaths[,i_sim]>0))
      calf_birth_days <- which(calf_age[,i_sim]==1)
      if (length(calf_birth_days)>1) {
          birth_intervals <- sapply(2:length(calf_birth_days), function(X) calf_birth_days[X]-calf_birth_days[X-1])
          ICI[i_sim] <- mean(birth_intervals)/365
          # calculate survival rate for all calves in their first year
          ####### this needs modification to prevent calf_birth_days+365 > max_days[i_sim] 
          #calf_survivors1 <- length(which(calf_age[(calf_birth_days+365),i_sim]>0))
          #calf_survival1[i_sim] <- calf_survivors1/length(calf_birth_days)
          }
      #calf_survival[i_sim] <- 2*repro_success[i_sim]/calves_born[i_sim]
      }


# end of loop over multiple females


  
