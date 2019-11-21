# -------------------------------------------------------------
# Setting up R Environment for Sensitivity Analysis - Sobol'
# -------------------------------------------------------------

# install packages - only run once
install.packages("dplyr")
install.packages("fast") # FAST method
install.packages("forcats")
install.packages("ggplot2")
install.packages("hydroGOF")
install.packages("sensitivity") # Sobol' method
install.packages("tidyr")
install.packages("purrr")
install.packages("devtools")
install.packages("SWATdata")
install.packages("EcoHydRology")
install.packages("SWATmodel")
devtools::install_github("chrisschuerz/SWATplusR")
#devtools::install_github("chrisschuerz/SWATplusR", ref = "dev") # development version

# load packages
library(dplyr)
library(fast)
library(forcats)
library(ggplot2)
library(hydroGOF)
library(lubridate)
library(purrr)
library(sensitivity)
library(tidyr)
library(SWATplusR)
library(SWATdata)
library(EcoHydRology)
library(SWATmodel)


# ------------------------------------------------------------
# load the data ... insert your own data here
# ------------------------------------------------------------

# where to store demo ... copy-paste your own filepath
demo_path <- "C:/Users/echelsvi/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/SensitivityAnalysis_CloudCreek/" #forward slashes
demo_path2 <- "C:/Users/echelsvi/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/SensitivityAnalysis_CloudCreek/swatplus_rev57_demo"
# load SWAT+ demo 
path_plus <- load_demo(dataset = "project", version = "plus", path = demo_path, revision = 57)
demo_path3 <- "C:/Users/echelsvi/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/SensitivityAnalysis_CloudCreek/swatplus_rev57_demo"
# load the dataset, and limit time period
q_obs <- load_demo(dataset = "observation")
q_obs <- filter(q_obs, date >= ymd("2003-01-01"),
                date <= "2007-12-31")

# what does the data look like?
head(q_obs) #structure of tibble
dim(q_obs) #dimensions
plot(q_obs, type = "l") #time series plot 

run_swatplus(project_path = demo_path3, #path to SWAT project folder
             output =list(q_sim = define_output(file = "channel", #variables to extract from SWAT model runs
                                                variable = "flo_out",
                                                unit = 1)))#,
             parameter = par,
             start_date = "2000-01-01", #SWAT simulation start date
             end_date = "2007-12-31", #SWAT simulation end date
             years_skip = 3, #number of simulation years that are skipped b4 writing SWAT outputs
             n_thread = 4, #number of threads for parallel model run
             add_date = FALSE) #add a date column to every simulation output table?

# ------------------------------------------------------------
# set-up for sensitivity analysis
# ------------------------------------------------------------

# name the parameters (will be used for sensitivity analysis) ...insert your parameters
#   parameter name
#   change = relchg (alter param by fraction of initial param value), pctchg (alter param by % of initial param value), 
#             abschg (add an absolute value to the initial param value), absval (replace param by absolute value)
#   can add more | for additional param restraints (subbasins, land uses, soil types, etc.)
#             ex.) par_mod <- c("CN2.mgt|change = abschg | sub == 1:5 | luse %in% c('WWHT', 'CORN') " = 5)
par_names <- c("cn2.hru | change = abschg", 
               "lat_ttime.hru | change = absval",
               "lat_len.hru | change = absval",
               "k.sol | change = pctchg",
               "z.sol | change = pctchg",
               "esco.hru | change = absval",
               "epco.hru | change = absval")



# create function that returns a scalar variable for which the sensitivity is assessed
#     par = SWAT model parameters (vector or tibble)
#     obs = df of observed values to be used for NSE

swat_sobol <- function(par, obs) {
  
  #implement sampled parameter combos and simulate daily discharges for start-end period
  q_sim <- run_swatplus(project_path = proj_path, #path to SWAT project folder
                        output =list(q_sim = define_output(file = "channel", #variables to extract from SWAT model runs
                                                           variable = "flo_out",
                                                           unit = 1)),
                        parameter = par,
                        start_date = "2000-01-01", #SWAT simulation start date
                        end_date = "2007-12-31", #SWAT simulation end date
                        years_skip = 3, #number of simulation years that are skipped b4 writing SWAT outputs
                        n_thread = 4, #number of threads for parallel model run
                        add_date = FALSE) #add a date column to every simulation output table?
  
  #NSE to evaluate simulations for the start-end period
  nse_q <- map_dbl(q_sim$simulation$q_sim/8.64, ~ NSE(.x, obs))
  return(nse_q)
}

# use NSE to identify relevant parameters and/or parameters bounds
dotty_iter1 <- q_iter1$parameter$values %>%
  mutate(nse = nse_iter1) %>%
  filter(nse > -5) %>%
  gather(key = "par", value = "val", -nse)

ggplot(data = dotty_iter1) +
  geom_point(aes(x = val, y = nse)) +
  facet_wrap(.~par, ncol = 3, scales = "free") +
  theme_bw()


# list bounds of each parameter, based off dotty plots
par_bound <- tibble("cn2.hru | change = abschg" = c(-15, 10),
                    "lat_ttime.hru | change = absval" = c(0.5, 50),
                    "lat_len.hru | change = absval" = c(10, 100),
                    "k.sol | change = pctchg" = c(-50, 50),
                    "z.sol | change = pctchg" = c(-50, 50),
                    "esco.hru | change = absval" = c(0, 1),
                    "epco.hru | change = absval" = c(0, 1))
#View(par_bound)



# ------------------------------------------------------------
# Sobol' sensitivity analysis
# ------------------------------------------------------------

# create 2 random samples w/ same sample size for parameters 
n_par  <- 7 #number of parameters in sensitivity analysis
n_samp <- 500 #sample size
x1 <- data.frame(matrix(runif(n_par * n_samp), nrow = n_samp)) %>% #uniform distrbution vector with n_par*n_samp obs.
  set_names(., names(par_bound)) %>% #name x1 using parameter names
  map2_dfc(., par_bound, ~ (.x * (.y[2] - .y[1]) + .y[1]))
#View(x1)
x2 <- data.frame(matrix(runif(n_par * n_samp), nrow = n_samp)) %>%
  set_names(., names(par_bound)) %>%
  map2_dfc(., par_bound, ~ (.x * (.y[2] - .y[1]) + .y[1]))
#View(x2)


# model function
#swat_nse <- enter your model here

# run Sobol
sens_sobol <- sobol(model = swat_nse, #enter your model here
                    X1 = x1, #first random sample
                    X2 = x2,#second random sample
                    obs = q_obs$discharge, #observations
                    nboot = 100) #number of bootstrap replicates


# plot Sobol
plot_sobol <- sens_sobol$S %>%
  mutate(parameter = rownames(.)) %>%
  mutate(parameter = factor(parameter) %>% fct_reorder(., original))
ggplot(data = plot_sobol) +
  geom_pointrange(aes(x = parameter, y = original ,
                      ymin = `min. c.i.`, ymax = `max. c.i.`)) +
  coord_flip() +
  xlab("Parameter") +
  ylab("Sensitivity") +
  theme_bw()




# ------------------------------------------------------------
# the end
# ------------------------------------------------------------