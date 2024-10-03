library(tidyr)
library(dplyr)

rm(list=ls())


#Build simple data
lbin <- 5
x <- read.csv("data/ch_2023.csv")

LW <- x %>%
  filter(Length < 180) %>%
  select(Tagged, as.Smolt,	As.Adult.ballard, Year, ReleaseSite, Length, ReleaseWeek) %>%
  filter(ReleaseSite %in% c("LWBEAR", 'LWCEDR')) %>%
  group_by(Tagged, as.Smolt,	As.Adult.ballard, Year, ReleaseSite) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(idx = row_number()) %>%
  mutate(y_i = factor(Year, levels = sort(unique(LW$Year)))
         # ,l_i = factor(ceiling(Length/lbin)*lbin,
         #              levels = sort(unique(ceiling(LW$Length/lbin)*lbin)))
  )

MR_settings <- list(
  data = na.omit(data.frame(LW)),
  frms = list(
    phi = 'state ~ (1|ReleaseSite:y_i)', #Survival
    p = 'state ~  (1|y_i)', #Detection probability for smolts at the locks.
    lam = 'state ~  (1|ReleaseSite:y_i)'  #Dummy last location, The locks as adults
  ),
  make_sd = TRUE
)

tmb_list <- new("tmb_list",
                raw_data = na.omit(data.frame(LW)), #as.data.frame(do.call(cbind,Orange)),
                MR_settings = MR_settings)
tmb_obj <- make_tmb_lists(tmb_list)
x <- build_tmb(tmb_obj)

plot(obj = x,
     type = "re",
     process = 'phi',
     factor = NULL,
     var = NULL)
AIC(x)


