#==============================================================================#
#                                                                              #
#               Statistical tests - Effect of the time of the year             #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(lubridate)
library(EMT)
library(RVAideMemoire)
Sys.setenv(LANG = "en")

# A. Count the roadkills -----------------------------------------------------
hy_carcasses <- read_delim("data/01_hy_carcasses.csv",
                           ";", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(date_obs = ymd(paste0(substr(date_obs, 1, 4), "-", substr(date_obs, 5, 6), "-", substr(date_obs, 7, 8))))

# ~ 1. Temporal ----------------------------------------------------------------

# ~~~ a. During each season ----------------------------------------------------

{ # Dry and wet seasons from the Serengeti IV book
  start_transit_1 <- as.Date("1989-05-11")
  end_transit_1 <- as.Date("1989-07-31")
  start_dry <- as.Date("1989-08-01")
  end_dry <- as.Date("1989-10-31")
  start_transit_2 <- as.Date("1989-11-01")
  end_transit_2 <- as.Date("1989-12-20")
  start_wet <- as.Date("1989-12-21")
  end_wet <- as.Date("1990-05-10")
  nbr_days <- 365
}

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

{ # Dry and wet seasons from the Serengeti IV book
  start_transit_1 <- as.Date("1989-05-11")
  end_transit_1 <- as.Date("1989-07-31")
  start_dry <- as.Date("1989-08-01")
  end_dry <- as.Date("1989-10-31")
  start_transit_2 <- as.Date("1989-11-01")
  end_transit_2 <- as.Date("1989-12-20")
  start_wet <- as.Date("1989-12-21")
  end_wet <- as.Date("1990-05-10")
  nbr_days <- 365
}

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
  filter(location_certainty_score >= 0.75,
         year != 2020, year != 2021,        # Remove the years 2020-2021 as they were incomplete because of covid-19
         year != 2023,                      # Remove the years 2023 as it was incomplete because of record termination
         age %in% c("adult", "unknown", "subadult"),
         herds_position %in% c("south east", "north west"))



# B. Do the tests ------------------------------------------------------------

# ~ a. Temporal pattern ------------------------------------------------------

# Season (with transit)
proba2 <- c(length(seq(start_dry, end_dry, by = "day"))/nbr_days, # North west
            # South east
            length(seq(start_wet, end_wet, by = "day"))/nbr_days,
            #Transit
            (nbr_days - length(c(seq(start_dry, end_dry, by = "day"),
                                 seq(start_wet, end_wet, by = "day"))))/nbr_days)

RVAideMemoire::multinomial.theo.multcomp(x = hy_carcasses_season$n,
                                         p = proba2,
                                         prop = FALSE,
                                         p.method = "fdr")

# Season (with transit but binomial)
proba3 <- c(length(seq(start_dry, end_dry, by = "day"))/(length(seq(start_dry, end_dry, by = "day")) +
                                                           length(seq(start_wet, end_wet, by = "day"))), # North west
            # South east
            length(seq(start_wet, end_wet, by = "day"))/(length(seq(start_dry, end_dry, by = "day")) +
                                                           length(seq(start_wet, end_wet, by = "day"))))

dbinom(x = hy_carcasses_season$n[hy_carcasses_season$herds_position == "north west"],
       size = sum(hy_carcasses_season$n[1:2]),
       prob = proba3[1])
sum(hy_carcasses_season$n[1:2])


# Month
dmultinom(x = hy_carcasses_month$n,
          prob = c(31/365, 28/365, 31/365, 30/365,
                   31/365, 30/365, 31/365, 31/365,
                   30/365, 31/365, 30/365, 31/365))

RVAideMemoire::multinomial.theo.multcomp(x = hy_carcasses_month$n,
                                         p = c(31/365, 28/365, 31/365, 30/365,
                                               31/365, 30/365, 31/365, 31/365,
                                               30/365, 31/365, 30/365, 31/365),
                                         prop = FALSE,
                                         p.method = "fdr")

sum(hy_carcasses_month$n)

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

# With transit
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
  filter(location_certainty_score >= 0.75,
         year != 2020, year != 2021, year != 2023,  # Remove years with incomplete seasons because of covid 19
         age %in% c("adult", "unknown", "subadult")) 
# lat
kruskal.test(formula = lat ~ herds_position,
             data = hy_carcasses_season_spatial)
FSA::dunnTest(lat ~ herds_position, data = hy_carcasses_season_spatial)

# long
kruskal.test(formula = long ~ herds_position,
             data = hy_carcasses_season_spatial)
FSA::dunnTest(long ~ herds_position, data = hy_carcasses_season_spatial)
