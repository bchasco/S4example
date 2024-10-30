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
  group_by(!!!syms(locs),rw_i, ReleaseSite,Year) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(idx = row_number()) %>%
  mutate(y_i = factor(Year, levels = sort(unique(Year)))
  )

MR_settings <- list(
  data = na.omit(data.frame(LW)),
  locs = locs,
  states = states,
  model_type = "CJS",
  frms = list(
    # phi = 'state ~ 1', #Survival
    phi = 'state ~ poly(rw_i,2) + (1 | ReleaseSite:y_i)', #Survival
    p = 'state ~  (1 | y_i)', #Detection probability
    lam = 'state ~ (1 | y_i)'  #joint last location
  ),
  make_sd = TRUE
)


tmb_list <- new("tmb_list",
                MR_settings = MR_settings)
tmb_obj <- make_tmb_lists(tmb_list)
st <- Sys.time()
tmb_model <- build_tmb(tmb_obj)

print(tmb_model@TMB$rep$gamma[3,,])
print(tmb_model@TMB$rep$omega[3,,])
print(round(plogis(tmb_model@TMB$opt$par),2))
print(Sys.time()- st)

# AIC(tmb_model)
# source('plots.r')
#
