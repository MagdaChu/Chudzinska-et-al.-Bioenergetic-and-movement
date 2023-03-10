

##################################################
###### SET KEY PARAMETERS

max_age <- 30

# set assumed reserve level of calves at weaning 
rho_j <- 0.3243

# steepness of assimilation response
#eta <- 15
eta <- 20
# relative cost of maintaining reserves
Theta_F <- 0.2
# shape parameter for starvation-induced reduction in milk delivery
# xi_m <- -2
xi_m <- -3
# shape parameter for calf foraging efficiency
#gamma <- 4
# age (days) at which resource foraging becomes 50% efficient
#Tr <- 150
#Tr <- 130
# set length of gestation
Tp <- 305
#Tp <- 335
# probability that implantation will occur
#fert_success <- 0.95
fert_success <- 1

# set foetal mortality.  This value is based on Murphy et al 2015 
foetal_mortality <- 0.197
# parameter defining strength of starvation related mortality
#mu_s <- 0.2
#mu_s <- 0.3

# field metabolic maintenance scalar
Sigma_M <- Sigma_M_full*0.294
# length of lactation. Assumed to be 8 months
Tl <- 250
# set age of female (in days) at start of simulation
age_day1 <- Tl+1
######################################################
# calculate maximum age in days 
maximum_age <- max_age*365
# set maximum duration of simulation
maximum_days <- maximum_age + age_day1

##################################################
# create vector of dates coded as Julian days, starting on 1 June + Tl, no leap years
birthday <- 152 # 01/06
day1 <- birthday + Tl - 365
julian_days <- c(day1:365)
remaining <- length(julian_days) + 1
julian_days[remaining:365] <- c(1:(day1-1))
# replicate the Julian days vectors
max_years <- ceiling(maximum_days/365) + 1
julian_days <- rep(julian_days,max_years)

# variations in target body condition over year

#rho_max <- 0.45
rho_max <- 0.35
#rho_min <- 0.3
rho_min <- 0.24
rho_annual <- rep(rho_max, 365)
# steady decline in rho from 1 March (Julian 60) to 1 July (Julian 182)
rho_annual[60:182] <- rho_max - (rho_max-rho_min)*c(1:123)/(182-59)
# steady rise in rho from 2 July (Julian 183) to 1 December (Julian 335)
#rho_annual[183:335] <- rho_min + (rho_max-rho_min)*c(1:153)/(335-182)
####### 9 March modification
# rho_min remains low until 15 September (Julian 258) then rises to 1 December (Julian 335)
rho_annual[183:258] <- rho_min
rho_annual[259:335] <- rho_min + (rho_max-rho_min)*c(1:77)/(335-258)

# annual cycle in catabolic and anabolic reserve conversion efficiencies in MJ/kg
epsi_minus_annual <- rep(20,365)
epsi_minus_annual[60:182] <- 35
epsi_plus_annual <- rep(55,365)
epsi_plus_annual[60:182] <- 28

# rescale rho and epsi's so that they match females time scale and calculate rho_s
rho_t <- rho_annual[julian_days[1:365]]
epsi_minus <- epsi_minus_annual[julian_days[1:365]]
epsi_plus <- epsi_plus_annual[julian_days[1:365]]

rho_t <- rep(rho_t,max_years)
epsi_minus <- rep(epsi_minus, max_years)
epsi_plus <- rep(epsi_plus, max_years)
# calculate starvation threshold
rho_s <- rho_t - 0.1
#rho_s <- rho_t - 0.11
#########################################################

#######################################################
######### GROWTH PARAMETERS FOR FEMALE ################
# length at infinity 
Linf <- 160
# length at birth 
Lb <- 70
# growth rate 
k <- 0.0015
# mass length scaling 
omega1 <- 5.9*10^(-5)
# mass length scaling exponent
omega2 <- 2.67

# energy efficiency per unit structural mass 
Sigma_G <- 30 

############## calculate length at age
ages <- c(1:(maximum_days+1))
La <- Linf - (Linf-Lb)*exp(-k*ages)
# calculate core weight at age 
Sa <- omega1*La^omega2
######### calculate growth costs for all age classes
CGa <- Sigma_G*(Sa[2] - Sa[1])
CGa[2:maximum_days] <- sapply(2:maximum_days, function(X) Sigma_G*(Sa[X+1] - Sa[X]))

#########################################################
############## LACTATION PARAMETERS #####################
# day of lactation when female contribution begins to decrease
#Tn <- 120
# shape parameter for decrease in milk supply with calf age
xi_c <- 0.5

### calculate change in % of calf's demand delivered by female over course of lactation
# extended beyond Tl to allow modelling of calf's growth post weaning if required
milk_days <- c(1:(Tl+365))
temp_prop <- (1 - (milk_days-Tn)/(Tl-Tn))/(1 - xi_c*(milk_days-Tn)/(Tl-Tn))
milk_prop <- ifelse(temp_prop>1,1,temp_prop)
milk_prop[(Tl+1):(Tl+365)] <- 0

# milk equivalent of R
phi_L <- 3.55
# efficiency of conversion of mother's reserves to milk
Sigma_L <- 0.86

# calculate feeding efficiency on each day over all possible ages
assim_effic <- sapply(1:maximum_days, function(X) X^gamma/(Tr^gamma + X^gamma))

##################### PREGNANCY PARAMETERS #############
# minimum age for reproduction
min_age <- 3
min_age <- min_age*365

# set period (in Julian days) during which implantation can occur, can be +/- 10 days either side of 31 July
implant_start <- birthday - Tp + 365 - 10
implant_period <- 20
# if rho is above threshold for pregnancy in this period, implantation can occur
starts <- which(julian_days==implant_start)
# "preg_possible" is 1 on days on which implantation can occur
# extends beyond max_age to allow modelling of calf survival beyond death of female
preg_possible <- rep(c(0),maximum_days)
for(i_preg in 1:length(starts)){preg_possible[starts[i_preg]:(starts[i_preg]+implant_period)] <- 1}
preg_possible[1:min_age] <- 0

###### calculate length and weight at age for foetus (assumes linear growth in length)
gest_t <- c(1:(Tp+1))
# NB L_foetus[Tp+1] is calculated solely for the purpose of determining growth costs
# L_foetus[Tp] is equivalent to La[0], i.e. length at birth in the cB growth curve
L_foetus <- gest_t*Lb/Tp
S_foetus <- omega1*L_foetus^omega2

# calculate maintenance cost of foetal tissue for a 10 year old female
CM_foetus <- sapply(1:Tp, function(X) Sigma_M*((Sa[(10*365+X)]+S_foetus[X])^0.75 - (Sa[(10*365+X)])^0.75))
# calculate growth cost of foetus
CG_foetus <- 0
CG_foetus[1] <- Sigma_G*S_foetus[2]
CG_foetus[2:Tp] <- sapply(2:Tp, function(X) Sigma_G*(S_foetus[X+1] - S_foetus[X]))
# calculate calf's total reserves at birth = minimum value of rho_s
F_start <- rho_s[365-age_day1]*S_foetus[Tp+1]/(1 - rho_s[365-age_day1])

# calculate daily and total cost of gestation
Cost_gestation <- sum(CG_foetus[1:Tp]) + sum(CM_foetus[1:Tp]) + max(epsi_plus)*F_start
# calculate extra cost of adding F_start and spread this over entire gestation period
extra_fat_cost <- Cost_gestation/(Cost_gestation - max(epsi_plus)*F_start)
CGest <- (CG_foetus + CM_foetus)*extra_fat_cost

# calculate foetal growth cost part of threshold reserves for start of pregnancy
#F_neo <- Cost_gestation/max(epsi_minus)

######################################################
######################################################
#F_neo <- 0.6*sum(CG_foetus[1:Tp])/max(epsi_minus)
F_neo <- sum(CG_foetus[1:Tp])/max(epsi_minus)
# F_neo <- (sum(CG_foetus[1:Tp])+max(epsi_plus)*F_start)/max(epsi_minus)
# calculate survival rates for foetus
daily_foetal_surv <- (1-foetal_mortality)^(1/Tp)
cum_Surv_foetus <- sapply(1:maximum_age,function(i) daily_foetal_surv^i)
cum_Surv_foetus[maximum_age] <- 0

# age-specific cumulative survival based on Sinclair et al. (2020) using Siler function in Barlow & Boveng
W_surv <- 0.8455
W_surv[2:5] <- 0.85
W_surv[6:(max_age-10)] <- 0.925
first_days <- seq(1,(maximum_age-1), 365)
cum_Surv <- 1
a2 <- 2.2
a1 <- 10
a3 <- 0.01
b1 <- 24
b3 <- 8
index <-(1:(30*365))/(30*365)
cum_Surv[1:length(index)] <- exp(-a2*index)*exp(-(a1/b1)*(1-exp(-b1*index)))*exp((a3/b3)*(1-exp(b3*index)))
cum_Surv[length(cum_Surv)] <- 0

f.foetusandcalf_life <- function(i){
  # determine foetal and calf life expectancy
  # with these values, ~20% of calves die from causes other than starvation before weaning
  check <- runif(2)
  life_expect_f <- min(which(cum_Surv_foetus<=check[1]))
  if (life_expect_f <= Tp) {life_expect_c <- 0} else 
    {life_expect_c <- min(which(cum_Surv<=check[2]))
     life_expect_f <- Tp}
  # model calf growth and energy expenditure up to age 1 year (less 20 day, to prevent any overlap in calves)
   if (life_expect_c > (365-20)) life_expect_c <- (365-20)
  # return values
  list (life_expect_f = life_expect_f, life_expect_c = life_expect_c)
  }

f.calf_starve_death <- function(rhoC,rho_s){
  # determine whether a calf dies (return 0) as a result of starvation on a particular day
  mort <- mu_s*(rho_s/rhoC - 1)
  rbinom(1,1,exp(-mort))
  }

##### parameters for beta distribution if modelling stochasticity in resources
mu <- 0.25
SD <- 0.075
b <- mu*(1-mu)^2/SD^2 + mu - 1
a <- mu*b/(1-mu)

