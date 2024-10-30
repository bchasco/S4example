library(tidyr)
library(dplyr)

rm(list=ls())


# #Build simple data
locs <- c('Tagged', 'as.Smolt',	'As.Adult.ballard')
states <- c("Alive", "Dead")

LW <- read.csv("data/combinedData.csv") %>%
  filter(Length < 180) %>%
  filter(Year < 2022) %>%
  mutate(Species = ifelse(Species==1,"Chinook","Coho")) %>%
  filter(Species == "Chinook") %>%
  mutate(temp_i = factor(round(LKWA10mTemp), levels = sort(unique(round(LKWA10mTemp))))) %>%
  mutate(rw_i = factor(ReleaseWeek, levels = sort(unique(ReleaseWeek)))) %>%
  select(Tagged, as.Smolt,	As.Adult.ballard, Year, Species, ReleaseSite, rw_i, temp_i) %>%
  filter(ReleaseSite %in% c("LWCEDR","LWBEAR", "LWBLRD", "LWISQH", "LWKENM")) %>%
  group_by(!!!syms(locs), ReleaseSite, rw_i, Year) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(idx = row_number()) %>%
  mutate(y_i = factor(Year, levels = sort(unique(Year))),
         init = 1, #everyone starts out alive
         init_state = 1,
         init_loc = 1,
         int = 1
  )


MR_settings <- list(
  data = na.omit(data.frame(LW)),
  locs = locs,
  states = states,
  model_type = "MS",
  frms = list(
    phi = 'int ~ (1|y_i)', #Survival
    p = 'int ~  (1|y_i)', #Detection probability
    lam = 'int ~ 1'  #joint last location
  ),
  make_sd = FALSE
)


tmb_list <- new("tmb_list",
                MR_settings = MR_settings)
tmb_obj <- make_tmb_lists(tmb_list)
tmb_model <- build_tmb(tmb_obj)

print(round(plogis(tmb_model@TMB$opt$par),2))
print(tmb_model@TMB$sd)
print(apply(tmb_model@TMB$rep$gamma,c(2,3,4),mean))
