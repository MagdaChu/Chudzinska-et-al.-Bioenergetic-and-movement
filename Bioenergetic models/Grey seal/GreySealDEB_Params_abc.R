

##################################################
###### SET KEY PARAMETERS

# Kappa determines what proportion of energy intake is diverted to growth if total is insufficient to cover all costs
# Kappa = 0 gives complete priority to growth, Kappa = 1 gives complete priority to maintenance
#Kappa <- 0.7
# steepness of assimilation response
eta <- 40
# relative cost of maintaining reserves
Theta_F <- 0.2

# shape parameter for starvation-induced reduction in milk delivery
xi_m <- -2
# length of lactation from Hall & Russell
Tl <- 18

# duration of post-weaning fast
pw_fast <- 21
# shape parameter for calf foraging efficiency
#upsilon <- 0.95
# age (days) at which resource foraging becomes 50% efficient
#Tr <- 75

# length of gestation
Tp <- 240
# probability that implantation will occur
fert_success <- 1

foetal_mortality <- 0.0
# starvation threshold
#rho_s <- 0.1

# parameter defining strength of starvation related mortality
#mu_s <- 0.2
# field metabolic maintenance scalar
Sigma_M <- Sigma_M_full*0.294

# set age of female (in days) at start of simulation
age_day1 <- 366
# set maximum simulated age for calf (20 days less than 1 year to prevent overlap of simulated calves)
max_age_calf <- 365-20

# maximum age (years) of female
max_age <- 40
maximum_age <- max_age*365+2
# set maximum duration of simulation to allow birth and survival of last foetus to be determined
maximum_days <- maximum_age

# birthday_temp is mean pupping date for colony, in this case I've used 23 November for Farne Islands
# I've used 2021 because it isn't a leap year!
birthday_temp <- "23/11/2021"
# standard deviation for pupping date
sd_birthday <- 0
birthday_temp <- strptime(birthday_temp,format="%d/%m/%Y")
# convert birthday to day of year, +1 because 1 January = Julian day 0
mean_birthday <- birthday_temp$yday + 1
# set day of year when implantation occurs 
implant_day <- mean_birthday - Tp
if(implant_day<0) implant_day <- implant_day +365
# set day of pregnancy when female decides whether or not to continue
#decision_day <- 200
decision_day_of_year <- implant_day + decision_day - 1 
if (decision_day_of_year > 365) decision_day_of_year <- decision_day_of_year - 365

##################################################
# calculate feeding efficiency on each day over all possible ages
assim_effic <- sapply(1:maximum_age, function(X) X^upsilon/(Tr^upsilon + X^upsilon))
# assume no feeding during lactation and post-weaning fast
assim_effic[1:(pw_fast+Tl)] <- 0
# reduce assimilation efficiency during moult to take account of extra time spent hauled out
moult_reduction <- 0.5
moult_duration <- 20
moult_start <- as.integer(moult_duration/2)
implant_age <- seq((2*365-Tp),maximum_age,365)
moult_effect <- rep(1,maximum_age)
for (i in 1:length(implant_age)) {
  moult_effect[(implant_age[i]-moult_start):(implant_age[i]+moult_start)] <- (1-moult_reduction)
  }
assim_effic <- assim_effic*moult_effect

############# input rho and epsi values #########
# Target condition for adults
rho <- 0.5
# Target condition for pups calculated on assumption that all of mass gain by calves during lactation is reserves
rho_C <- 0.7
# set assumed reserve level on age_day1
# this value is derived from earlier simulations
rho_start <- 0.24
# Energy density of reserve tissue from Reilly et al
epsi <- 25.8
# Energy required to build up pup reserves
# value derived from Reilly et al (1996) - see report
epsi_plus_pups <- epsi + 2.665456
# Energy required to build up adult reserves, value from Hin et al.
epsi_plus <- epsi*55/40
# Energy provided by metabolizing reserves, value from Hin et al.
epsi_minus <- epsi*0.9


#######################################################
######### GROWTH PARAMETERS FOR FEMALE ################
# taken from McLaren 1993 La[t] = Linf*(1-exp(-k*(t-x0*365)))^b
Linf <- 184
k <- 0.182/365
x0 <- -0.4*365
b <- 0.27

# mass length scaling 
omega1 <- 1.5*10^(-4)
# mass length scaling exponent
omega2 <- 2.575

# energy efficiency per unit structural mass 
Sigma_G <- 30 

############## calculate length at age
ages <- c(1:(maximum_age+1))
La <- Linf*(1 - exp(-k*(ages - x0)))^b
#  calculate predicted length at birth
Lb <- Linf*(1 - exp(-k*(- x0)))^b
# calculate structural mass at age 
Sa <- omega1*La^omega2
######### calculate daily growth rate
GR <- Sa[1] - omega1*Lb^omega2
GR[2:maximum_age] <- sapply(2:maximum_age, function(X) (Sa[X] - Sa[X-1]))

#########################################################
############## LACTATION PARAMETERS #####################
# day of lactation when female contribution begins to decrease
Tn <- as.integer(0.95*Tl)
# shape parameter for decrease in milk supply with calf age
xi_c <- 0.99

### calculate change in % of calf's demand delivered by female over course of lactation
# extended beyond Tl to allow modelling of calf's growth post weaning if required
temp_prop <- rep(1,Tl)
milk_prop <- ifelse(temp_prop>1,1,temp_prop)
milk_prop[(Tl+1):(max_age_calf)] <- 0

# milk equivalent of R calculated from values in Reilly et al (1996)
phi_L <- (13.5*Tl + 31*epsi)/(Tl*Sa[9]^(2/3))
# efficiency of conversion of mother's reserves to milk
Sigma_L <- 0.86

##################### PREGNANCY PARAMETERS #############
# minimum age for reproduction
min_age <- 3
min_age <- min_age*365

###### calculate length and weight at age for foetus (assumes linear growth in length)
gest_t <- c(1:(Tp+1))
# NB L_foetus[Tp+1] is calculated solely for the purpose of determining growth costs
# L_foetus[Tp] is equivalent to La[0], i.e. length at birth in the cB growth curve
L_foetus <- gest_t*Lb/Tp
S_foetus <- omega1*L_foetus^omega2

# calculate daily growth rate of foetus
GR_foetus <- 0
GR_foetus[1] <- S_foetus[2]
GR_foetus[2:Tp] <- sapply(2:Tp, function(X) (S_foetus[X+1] - S_foetus[X]))

# calculate survival rates for foetus
daily_foetal_surv <- (1-foetal_mortality)^(1/Tp)
cum_Surv_foetus <- sapply(1:maximum_age,function(i) daily_foetal_surv^i)
cum_Surv_foetus[maximum_age] <- 0

# age-specific cumulative survival based on values in Thomas et al 2019 using Siler function in Barlow & Boveng
alpha1 <- 0.007
beta1 <- 0.01
alpha2 <- 1.3*10^(-4)
beta2 <- 0.1*10^(-6)
daily_mort <- alpha1*exp(-beta1*(1:maximum_age)) + alpha2*exp(beta2*(1:maximum_age))
cum_Surv <- cumprod(exp(-daily_mort))
cum_Surv[maximum_age] <- 0


f.foetusandcalf_life <- function(i){
  # determine foetal and calf life expectancy
  # with these values, ~20% of calves die from causes other than starvation before weaning
  check <- runif(2)
  life_expect_f <- min(which(cum_Surv_foetus<=check[1]))
  if (life_expect_f <= Tp) {life_expect_c <- 0} else 
    {life_expect_c <- min(which(cum_Surv<=check[2]))
     life_expect_f <- Tp}
  # calf growth and energy expenditure are only modelled up to max_age_calf
   if (life_expect_c > max_age_calf) life_expect_c <- max_age_calf
  # return values
  list (life_expect_f = life_expect_f, life_expect_c = life_expect_c)
  }

f.starve_death <- function(rho,rho_s){
  # determine whether an animal dies (return 0) as a result of starvation on a particular day
  mort <- mu_s*(rho_s/rho - 1)
  rbinom(1,1,exp(-mort))
  }

##### parameters for beta distribution if modelling stochasticity in resources
mu_beta <- 0.25
SD_beta <- 0.075
b_beta <- mu_beta*(1-mu_beta)^2/SD_beta^2 + mu_beta - 1
a_beta <- mu_beta*b_beta/(1-mu_beta)

###### parameter values for calculating P{pupping}
#b_pupping <- 0.17
# Isle of May values for PPM
a_pupping <- -12.2
b_pupping <- 0.13/1.4

# X is total weight of female on decision_day
f.prob_of_pupping <- function(X){
  p_pup <- exp(a_pupping+b_pupping*X)/(1+exp(a_pupping+b_pupping*X))
  # determine whether pregnancy will continue
  rbinom(1,1,p_pup)
  }

####### scalar to create cyclical variation in Rmean 
day <- c(0:(2*364+1))
# set amplitude of cycle (+/- 10% in this case) amplitude = 0 gives constant R
amplitude <- 0.1
scalar <- 1 + amplitude*sin(2*day*pi/365)
# we use mean_birthday and an offset term to determine when during the simulated year "scalar" has its mean value
# this coincides with day[1] and day[366]
# we want R to have its mean value (or as close as we can manage) on mean_birthday and then decrease
offset <- 184
# if we want the maximum value of R to occur on mean_birthday
# offset <- 92
# if we want R to be at its mean value on mean_birthday and to increase after this
#offset <- 0
start_day <- mean_birthday - offset
start_day <- ifelse (start_day > 0, start_day, start_day+365)
first_days <- c(start_day:365)
last_days <- c(1:(365-length(first_days)))
julian_scale <- rep(c(first_days,last_days),2)

birthday <- as.integer(rnorm(1,mean_birthday,sd_birthday))
#day 1 is first day of simulation
day1 <- birthday
if (day1 > 365) day1 <- day1 - 365
julian_days <- c(day1:365)
remaining <- length(julian_days) + 1
julian_days[remaining:365] <- c(1:(day1-1))
# replicate julian_days
max_years <- ceiling(maximum_days/365) + 1
julian_days <- rep(julian_days,max_years)

