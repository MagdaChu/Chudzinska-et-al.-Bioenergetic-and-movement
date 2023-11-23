
##################################################
###### SET KEY PARAMETERS
# Kappa determines what proportion of energy intake is diverted to growth if total is insufficient to cover all costs
# Kappa = 0 gives complete priority to growth, Kappa = 1 gives complete priority to maintenance
Kappa <- 0.2
KappaLact <- 0.6
# steepness of assimilation response (18) - grey seals have 20 but 2-25 was explored
eta <- 20
# relative cost of maintaining reserves (11)
Theta_F <- 0.2

# shape parameter for starvation-induced reduction in milk delivery (24,23 and 40) 
xi_m <- -2
# length of lactation (27 and 2)
Tl <- 23 # 24: Bowen 1992;Harkonen 1990 - 28 days; Moulbert et al 1993 - 24 days, Sauve et al 2014 - 33 days, Cordes et al 2013 18-23 days
# duration of post-weaning fast ()
pw_fast <- 15 # (Moulbert et al 1993)

# shape parameter for calf foraging efficiency (36A)
upsilon <- 1.8
# age (days) at which resource foraging becomes 50% efficient (20 and 36B)
Tr <- 65 
# length of gestation (1)
Tp <- 240
# probability that implantation will occur (part of 38?)
fert_success <- 1
# Foetal mortality arbitrarily set to 0.1 (JUST A STOCHASTIC EVENT, NOT RELATED TO CONDITIONS, set it to 0 for initial simulation)
foetal_mortality <- 0.1
# parameter defining strength of starvation related mortality (37 and 21)
mu_s <- 0.2

# field metabolic maintenance scalar (9 and 29, 0.294 = K (8))
Sigma_M <- 2*0.294
# set age of female (in days) at start of simulation
age_day1 <- 366 

# calculate maximum age in days 
max_age <- 30 # Hall 2019
maximum_age <- max_age*365
# set maximum duration of simulation
maximum_days <- maximum_age + age_day1
# set maximum simulated age for calf (20 days less than 1 year to prevent overlap of simulated calves)
max_age_calf <- 365-20

################# BIRTHING

# birthday_temp is mean pupping date for colony, in this case I've used 23 November for Farne Islands
# I've used 2021 because it isn't a leap year!
birthday_temp <- "17/06/2021" # (43) # Crdes et al 2013
# sd_birthday alloes birth day of individual females to be different
sd_birthday <- 0
birthday_temp <- strptime(birthday_temp,format="%d/%m/%Y")
# convert birthday to day of year, +1 because 1 January = Julian day 0
mean_birthday <- birthday_temp$yday + 1
# set day of year when implantation occurs 
implant_day <- mean_birthday - Tp
if(implant_day<0) implant_day <- implant_day +365
# set day of pregnancy when female decides whether or not to continue
decision_day <- 160
decision_day_of_year <- implant_day + decision_day - 1 
if(decision_day_of_year>365) decision_day_of_year <- decision_day_of_year-365


###############
# calculate feeding efficiency on each day over all possible ages (36), some seals may not reach max effiency until age 1+
assim_effic <- sapply(1:maximum_days, function(X) X^upsilon/(Tr^upsilon + X^upsilon))
# assume no feeding during lactation and post-weaning fast
assim_effic[1:(pw_fast+Tl)] <- 0
# reduce assimilation efficiency during moult to take account of extra time spent hauled out
moult_reduction <- 0.5
moult_duration <- 20
moult_start <- as.integer(moult_duration/2)
implant_age <- seq((2*365-Tp),maximum_days,365)
moult_effect <- rep(1,maximum_days)
for (i in 1:length(implant_age)) {
  moult_effect[(implant_age[i]-moult_start):(implant_age[i]+moult_start)] <- (1-moult_reduction)
}
assim_effic <- assim_effic*moult_effect

############# input rho and epsi values #########
# Target condition for adults
rho <- 0.43
# rho_C calculated on assumption that all of mass gain by calves during lactation is reserves
rho_C <- 0.55 # for pups (13) Jorgensen et al 2001 - 0.35, if we take the same assumption os for grey seals, that pups mainly increase in fat and not
# core mass, than they start at around 9kg when born and finish at about 20 when weaned (Harkonen and HeideJ 1990), 20-26kg (Moulbert et al 1993,2003, but this is Canada, but pups are bigger than so % is similar) meaning that 
# they gained 11 kg in fat which is 0.55 of total body mass
# set assumed reserve level on age_day1
# this value is a bit higher than value in Hall & McConnell (2007)
rho_start <- 0.20
# energy density of reserve tissue from Reilly et al (33 and 15)
epsi <- 25.8 # same as for grey now
# epsi_plus_pups is energy required to build up pup reserves
# value derived from Reilly et al (1996) - see report
epsi_plus_pups <- epsi + 2.665456 # (33 and 16) # additonal costs of accumulating reserves for pups
# equivalent value for adults, value from Hin et al.
epsi_plus <- epsi*55/40 # (33 and ???) # from appendix in Hin
# energy provided by metabolizing reserves, value from Hin et al.
epsi_minus <- epsi*0.9 # (33 and)

# input starvation threshold (14 and 32)
rho_s <- 0.1

#######################################################
######### GROWTH PARAMETERS FOR FEMALE most from 28 ################
# taken from McLaren 1993 La[t] = Linf*(1-exp(-k*(t-x0*365)))^b
Linf <- 140.5 # max length 4
#k <- 0.0012#0.52/365  #28 and 5 to get 0.0012 in HS table, Hall et al 2019
#x0 <- -0.35*365 # for b= 1: x0 = -1.65*365 # #28 bylo -0.35
#b <- 0.26 #1 #28 0.30
k=0.441/365 # as in Hall 2019 0.441
x0 <- -2.02*365 # as in Hall 2019

# mass length scaling (28 and 6)
omega1 <- 0.000033 
# mass length scaling exponent (28 and 7)
omega2 <- 2.875 

# energy efficiency per unit structural mass (10 and 30)
Sigma_G <- 25 

############## calculate length at age - 28
ages <- c(1:(maximum_days+1))
#La <- Linf*(1 - exp(-k*(ages - x0)))^b

La <- Linf*(1 - exp(-k*(ages - x0))) # as in Hall 2012
#  calculate predicted length at birth, same as L0 (3)
Lb <- Linf*(1 - exp(-k*(- x0))) # as in Hall 2012
# calculate core weight at age 
Sa <- omega1*La^omega2
######### calculate daily growth rate
GR <- Sa[1] - omega1*Lb^omega2
GR[2:maximum_days] <- sapply(2:maximum_days, function(X) (Sa[X] - Sa[X-1]))

#########################################################
############## LACTATION PARAMETERS #####################
# day of lactation when female contribution begins to decrease
Tc <- as.integer(0.70*Tl) # (39 and 22, called there Tc);  # from Bowen 1992, Mulebert and Bowen 0.6 (24/30)
# shape parameter for decrease in milk supply with calf age
xi_c <- 0.99 
# day of lactation when female starts eating again
Lact_feed <- 10 # Boness et al 1994, Thompson et a; 1994), Bowen et al 2001: 8th day
# at what proportion of R females feed at Lact_feed till the end of lactation
# Bowen et al 2001 says that 30% of energy given to pups comes from external
# food (not stored fat) after females start foraging'
R_prop_lactation <- 0.6 
### calculate change in % of calf's demand delivered by female over course of lactation
# extended beyond Tl to allow modelling of calf's growth post weaning if required
temp_prop <- rep(1,Tl)

# milk_prop <- ifelse(temp_prop>1,1,temp_prop)
# milk_prop[(Tl+1):(Tl+365)] <- 0
# plot(milk_prop) # this may have to be changed

# Magda's alternative approach where provision decrease, assuming that harbour seal pup may start getting food before (Sauve et al. 2014)
# this does not make a huge difference
milk_prop_decrease <- 0.9
milk_prop <- ifelse(temp_prop>1,1,temp_prop)
milk_prop[Tc:length(milk_prop)] <- milk_prop_decrease
milk_prop[(Tl+1):(Tl+365)] <- 0
#plot(milk_prop)

# milk equivalent of R calculated from values in Reilly et al (1996)
daily_pup_ee <- 8.3 # in MJ/day assuming that it is 3.7 * BMR as describer by Reilly 1996 so and average mass of 15kg so 70*15^0.75*0.004184 #13.5 MJ/d Reilly 1996 for grey seal
#pup_mass_gain <- 31 # in kg at the end of the lactation assuming rho_pup
# Bowen et al. 2001 says 0.64 kg/day
pup_mass_gain_day <- 0.55#0.5
pup_mass_gain <- pup_mass_gain_day*Tl # assuming 0.5kg gain per day (Dube, Hammill & Barrette 2003;Muelbert, Bowen & Iverson 2003, Bowen et al 1992: 0.8kg/day)
phi_L <- (daily_pup_ee*Tl + pup_mass_gain*epsi)/(Tl*mean(Sa[1:Tl])^(2/3)) # (25 and 41)
# efficiency of conversion of mother's reserves to milk
Sigma_L <- 0.86 # (17 and 34)

##################### PREGNANCY PARAMETERS #############
# minimum age for reproduction
min_age <- 4 # Harkonen and Jorgensen 1990 and iPCoD helpfile Table1
min_age <- min_age*365

###### calculate length and weight at age for foetus (assumes linear growth in length)
gest_t <- c(1:(Tp+1))
# NB L_foetus[Tp+1] is calculated solely for the purpose of determining growth costs
# L_foetus[Tp] is equivalent to La[0], i.e. length at birth in the cB growth curve
L_foetus <- gest_t*Lb/Tp # this is backcalculating from lenmgth at birth, 
S_foetus <- omega1*L_foetus^omega2 # daily growth rate

# calculate daily growth rate of foetus
GR_foetus <- 0
GR_foetus[1] <- S_foetus[2]
GR_foetus[2:Tp] <- sapply(2:Tp, function(X) (S_foetus[X+1] - S_foetus[X]))

# calculate survival rates for foetus
daily_foetal_surv <- (1-foetal_mortality)^(1/Tp)
cum_Surv_foetus <- sapply(1:maximum_age,function(i) daily_foetal_surv^i)
cum_Surv_foetus[maximum_age] <- 0

# age-specific cumulative survival based on values in iPCoD help file using Siler function in Barlow & Boveng
# THIS HAS TO BE REGION ADJUSTED AS IT DIFFERS FOR HARBOUR SEALS
# alpha1 <- 0.008 # changed from 0.007 in grey seals
# beta1 <- 0.01
# alpha2 <- 1.3*10^(-4)
# beta2 <- 0.1*10^(-6)
#see 'survival_tuning_harbour_seal' (pup surv 0.5, adult 0.92, juv 0.78)
alpha1 <- 0.0022
beta1 <- 0.0019
alpha2 <- 2.1*10^(-4)
beta2 <- 0.3*10^(-6)

#values as for ann_surv <- 0.5, ann_surv[2:First_birth] <- 0.86; ann_surv[(First_birth+1):max_age] <- 0.96

# alpha1 <- 0.005
# beta1 <- 0.0055
# alpha2 <- 1.3*10^(-4)
# beta2 <- 3.6*10^(-9)

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
  # model calf growth and energy expenditure up to max_age_calf
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
mu_beta <- 0.55#0.25
SD_beta <- 0.075
b_beta <- mu_beta*(1-mu_beta)^2/SD_beta^2 + mu_beta - 1
a_beta <- mu_beta*b_beta/(1-mu_beta)


###### parameter values for calculating P{pupping} (38)
# to get b_pupping I had to find a value which would result in p_pup close to 0.99
# assuming mass_at_weaning =60 so average mass at weaning
b_pupping <- 0.21
# female mass at weaning that results in 0.5 probability of pupping
# to get skipping_point I calculated proportion of skipping_point used by John in comparion to mass of grey seal at weaning (88.5/143=0.62) and used this 0.62 to calculate skipping_point
# for harbour assuming weight at weaning = 60
skipping_point <- 45 #88.5 #(alpha) # 88.5 for grey seals see Fig 2 in SMout et al 2019
# this function calculates apparent weight at the end of lactation and determines whether or not pregnancy
# will continue, if it doesn't, calf_days will be reset to current age of foetus
# X is total weight of female on decision_day (180 for grey seals)

### new equation




### old equation
f.prob_of_pupping <- function(X){
  a_pupping <- -skipping_point*b_pupping
  mass_at_weaning <- X*(1 - 0.17*decision_day/Tp)
  p_pup <- exp(a_pupping+b_pupping*mass_at_weaning)/(1+exp(a_pupping+b_pupping*mass_at_weaning))
  # determine whether pregnancy will continue
  rbinom(1,1,p_pup)
  p_pup
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


