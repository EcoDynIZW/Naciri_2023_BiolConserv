#==============================================================================#
#                                                                              #
#               Statistical tests - Effect of the time of the year             #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(lubridate)
library(RVAideMemoire)
Sys.setenv(LANG = "en")

# A. Count the roadkills -----------------------------------------------------
hy_carcasses <- read_delim("data/01_hy_carcasses.csv",
                           ";", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(date_obs = ymd(paste0(substr(date_obs, 1, 4), "-", substr(date_obs, 5, 6), "-", substr(date_obs, 7, 8))))

# ~ 1. Temporal ----------------------------------------------------------------

# ~~~ a. During each season ----------------------------------------------------

hy_carcasses_season <- hy_carcasses %>%
  mutate(year = year(date_obs),
         herds_position = ifelse(date_obs %within% interval(ymd(paste0(year, "-08-01")),
                                                            ymd(paste0(year, "-10-31"))),
                                 "north west",
                                 ifelse(date_obs %within% interval(ymd(paste0(year, "-05-11")),
                                                                   ymd(paste0(year, "-07-31"))),
                                        "transit",
                                        ifelse(date_obs %within% interval(ymd(paste0(year, "-11-01")),
                                                                          ymd(paste0(year, "-12-20"))),
                                               "transit", "south east")))) %>%
  filter(year != 2020, year != 2021,    # Remove the years 2020-2021 as they were incomplete because of covid-19
         year != 2023) %>%              # Remove the years 2023 as it was incomplete because of record termination
  count(herds_position)
sum(hy_carcasses_season$n)


# ~~~ b. During each month -----------------------------------------------------
hy_carcasses_month <- hy_carcasses %>%
  mutate(month = month(date_obs),
         year = year(date_obs)) %>%
  filter(year != 2020, year != 2021, # Remove the years 2020-2021 as they were incomplete because of covid-19
         year != 2023) %>%           # Remove the years 2023 as it was incomplete because of record termination
  count(month)

sum(hy_carcasses_month$n)


# ~ 2. Spatiotemporal ----------------------------------------------------------

hy_carcasses_season_spatial <- hy_carcasses %>%
  mutate(year = year(date_obs),
         herds_position = ifelse(date_obs %within% interval(ymd(paste0(year, "-08-01")),
                                                            ymd(paste0(year, "-10-31"))),
                                 "north west",
                                 ifelse(date_obs %within% interval(ymd(paste0(year, "-05-11")),
                                                                   ymd(paste0(year, "-07-31"))),
                                        "transit",
                                        ifelse(date_obs %within% interval(ymd(paste0(year, "-11-01")),
                                                                          ymd(paste0(year, "-12-20"))),
                                               "transit", "south east")))) %>%
  filter(!is.na(long),
         year != 2020, year != 2021,        # Remove the years 2020-2021 as they were incomplete because of covid-19
         year != 2023,                      # Remove the years 2023 as it was incomplete because of record termination
         age %in% c("adult", "unknown", "subadult"),
         herds_position %in% c("south east", "north west"))



# B. Do the tests ------------------------------------------------------------

# ~ a. Temporal pattern ------------------------------------------------------

# By month
dmultinom(x = hy_carcasses_month$n,
          prob = c(31/365, 28/365, 31/365, 30/365,
                   31/365, 30/365, 31/365, 31/365,
                   30/365, 31/365, 30/365, 31/365))

multinomial.theo.multcomp(x = hy_carcasses_month$n,
                                         p = c(31/365, 28/365, 31/365, 30/365,
                                               31/365, 30/365, 31/365, 31/365,
                                               30/365, 31/365, 30/365, 31/365),
                                         prop = FALSE,
                                         p.method = "fdr")

sum(hy_carcasses_month$n)

# By season
{ # Dry and wet seasons from the Serengeti IV book
  start_dry <- as.Date("1989-08-01")
  end_dry <- as.Date("1989-10-31")
  start_wet <- as.Date("1989-12-21")
  end_wet <- as.Date("1990-05-10")
  nbr_days <- 365
}

proba3 <- c(length(seq(start_dry, end_dry, by = "day"))/(length(seq(start_dry, end_dry, by = "day")) +
                                                           length(seq(start_wet, end_wet, by = "day"))), # North west
            # South east
            length(seq(start_wet, end_wet, by = "day"))/(length(seq(start_dry, end_dry, by = "day")) +
                                                           length(seq(start_wet, end_wet, by = "day"))))

dbinom(x = hy_carcasses_season$n[hy_carcasses_season$herds_position == "north west"],
       size = sum(hy_carcasses_season$n[1:2]),
       prob = proba3[1])
sum(hy_carcasses_season$n[1:2])


# ~ b. Spatiotemporal pattern --------------------------------------------------

# lat
wilcox.test(formula = lat ~ herds_position,
            data = hy_carcasses_season_spatial,
            alternative = "two.sided", exact = TRUE)
nrow(hy_carcasses_season_spatial)
# long
wilcox.test(formula = long ~ herds_position,
            data = hy_carcasses_season_spatial,
            alternative = "two.sided", exact = TRUE)
