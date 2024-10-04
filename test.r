library(tidyr)
library(dplyr)

rm(list=ls())


# #Build simple data
# lbin <- 5
LW <- read.csv("data/ch_2023.csv")

LW <- LW %>%
  filter(Length < 180) %>%
  filter(Year < 2022) %>%
  # mutate(Length = scale(Length)) %>%
  select(Tagged, as.Smolt,	As.Adult.ballard, Year, ReleaseSite, Length, ReleaseWeek) %>%
  filter(ReleaseSite %in% c("LWBEAR", 'LWCEDR','LWISQH','LWBLRD', 'LWKENM')) %>%
  group_by(Tagged, as.Smolt,	As.Adult.ballard, Year, ReleaseSite, Length, ReleaseWeek) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(idx = row_number()) %>%
  mutate(y_i = factor(Year, levels = sort(unique(LW$Year)))
         # ,l_i = factor(ceiling(Length/lbin)*lbin,
         #              levels = sort(unique(ceiling(LW$Length/lbin)*lbin)))
  )
locs <- c('Tagged', 'as.Smolt',	'As.Adult.ballard')

# Wen <- read.csv('data/Wen.csv')
# Wen <- as.data.frame(do.call(cbind,Wen)) %>%
#   filter(ch_i.Trib==1)
# locs <- c('ch_i.Trib', 'ch_i.FirstTrap', 'ch_i.LwrWen', 'ch_i.MCN','ch_i.Bon', 'ch_i.TWX')

MR_settings <- list(
  data = na.omit(data.frame(LW)),
  locs = locs,
  frms = list(
    phi = 'state ~ poly(Length,2) + (1|ReleaseSite/y_i)', #Survival
    p = 'state ~  poly(ReleaseWeek,2) + (1|ReleaseSite:y_i)', #Detection probability for smolts at the locks.
    lam = 'state ~ (1|ReleaseSite/y_i)'  #joint probability of detection and survival for the  last location
  ),
  make_sd = TRUE
)

tmb_list <- new("tmb_list",
                MR_settings = MR_settings)
tmb_obj <- make_tmb_lists(tmb_list)
tmb_model <- build_tmb(tmb_obj)

tmb_model@TMB$opt
#
p <- plot(obj = tmb_model,
          type = "re",
          process = 'p',
          color = 'ReleaseSite',
          factor = NULL,
          var = NULL)
p <- plot(obj = tmb_model,
          type = "re",
          process = 'phi',
          color = 'ReleaseSite',
          factor = NULL,
          var = NULL)
p <- plot(obj = tmb_model,
          type = "re",
          process = 'lam',
          color = 'ReleaseSite',
          var = NULL)

AIC(tmb_model)


