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


# ~~~ a. By month --------------------------------------------------------------
daily_monitoring <- read_csv("data/monitoring_sessions_processed_1989_2022.csv")

daily_monitoring_month <- daily_monitoring %>%
  mutate(year = year(date),
         month = month(date)) %>%
  group_by(year, month) %>%
  summarise(monitored = mean(monitored))


kruskal.test(formula = monitored ~ month, data = daily_monitoring_month)

res <- dunn.test::dunn.test(x = daily_monitoring_month$monitored, g = daily_monitoring_month$month,
                     method = "bh",
                     kw = TRUE, # Do the Kruskal Wallis test
                     label = TRUE)

# ~~~ b. By season -------------------------------------------------------------
daily_monitoring <- read_csv("data/monitoring_sessions_processed_1989_2022.csv")

daily_monitoring_season <- daily_monitoring %>%
  filter(herds_position_w_transit %in% c("north west", "south east")) %>%  # Remove the "transit season"
  group_by(year, herds_position_w_transit) %>%
  summarise(monitored = mean(monitored))

wilcox.test(formula = monitored ~ herds_position_w_transit,
            data = daily_monitoring_season,
            alternative = "two.sided", exact = TRUE)



