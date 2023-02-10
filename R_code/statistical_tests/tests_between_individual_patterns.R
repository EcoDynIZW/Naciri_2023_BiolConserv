#==============================================================================#
#                                                                              #
#            Statistical test of the effect of age, sex and rank               #
#                                                                              #
#==============================================================================#

library(tidyverse)
Sys.setenv(LANG = "en")

# A. Effect of age =============================================================

# Proportion of each age class in the study population (calculated from the three
# study clans over the whole study period)

proportions <- c(0.6608927, 0.1320453, 0.2070620)    
names(proportions) <- c("Adults", "Subadults", "Cubs")

# ~ 2. Count roadkills of each age class --------------------------------------

hy_carcasses <- read_delim("data/01_hy_carcasses.csv",
                           ";", escape_double = FALSE, trim_ws = TRUE)
{counts <- hy_carcasses %>%
    count(age, sex)
  nb.ad.F <- as.numeric(counts[1,3])
  nb.ad.M <- as.numeric(counts[2,3])
  nb.ad <- as.numeric(counts[1,3] + counts[2,3] + counts[3,3])
  nb.sub <- as.numeric(counts[7,3] + counts[8,3] + counts[9,3])
  nb.cub <- as.numeric(counts[4,3] + counts[5,3] + counts[6,3])
  nb.unkn <- as.numeric(counts[10,3] + counts[11,3])
}


# ~ 3. Run statistical tests ---------------------------------------------------

# Adults and subadults VS cubs (without unknowns as adults)
dbinom(x =  nb.ad + nb.sub,
       size = nb.ad + nb.sub + nb.cub,
       prob = proportions[1] + proportions[2])
# p < 0.01 (0.0026)
nb.ad + nb.sub + nb.cub

# (with unknwons as adults)
dbinom(x = nb.ad + nb.sub + nb.unkn,
       size = nb.ad + nb.sub + nb.cub + nb.unkn,
       prob = proportions[1] + proportions[2])
# p = 0.001
nb.ad + nb.sub + nb.cub + nb.unkn


# B. Effect of sex =============================================================

# ~ 1. Calculate the average adult and subadult sex ratio ---------------------

proportions <- c(0.36202532, 0.29886742, 0.06795470, 0.06409061, 0.20706196)
names(proportions) <- c("Adult females", "Adult males", "Subadults females", "Subadults males", "Cubs")

ad.sex.ratio <- unname(proportions[1]/(proportions[1] + proportions[2]))

# ~ 2. Count individuals of each sex -------------------------------------------
hy_carcasses <- read_delim("data/01_hy_carcasses.csv",
                           ";", escape_double = FALSE, trim_ws = TRUE)

ad.counts <- hy_carcasses %>%
  filter(age == "adult",
         sex %in% c("F", "M")) %>%
  count(sex)
nb.ad.F <- as.numeric(ad.counts[1, 2])
N <- nb.ad.F + as.numeric(ad.counts[2, 2])

# Adults With individuals of unknwon age considered adults --
ad.counts.bis <- hy_carcasses %>%
  filter(age %in% c("adult", "unknown"),
         sex %in% c("F", "M")) %>%
  count(sex)
nb.ad.F.bis <- as.numeric(ad.counts.bis[1, 2])
N.bis <- nb.ad.F + as.numeric(ad.counts.bis[2, 2])


# ~ 3. Statistical tests -------------------------------------------------------

# ----------- Adults --
dbinom(x = nb.ad.F, size = N, prob = ad.sex.ratio)
# p = 0.0058
N

# With unknown
dbinom(x = nb.ad.F.bis, size = N.bis, prob = ad.sex.ratio)
# p = 0.0086
N.bis


# C. Effect of rank ============================================================

# All the carcasses of females from the three studied clans have collision
# certainty scores = 1.

hy_carcasses_clan_members <- read_delim("data/02_hy_carcasses_clan_members.csv",
                                        ";", escape_double = FALSE, trim_ws = TRUE)

hy_carcasses_clan_members <- hy_carcasses_clan_members %>%
  filter(sex == "F",
         age_category %in% c("adult", "subadult")) # remove cubs


# ~ Statistical tests ----------------------------

wilcox.test(hy_carcasses_clan_members$standardized_rank, mu = 0, exact = TRUE)
nrow(hy_carcasses_clan_members)
