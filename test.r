library(tidyr)
library(dplyr)

rm(list=ls())


# #Build simple data
lbin <- 5
locs <- c('Tagged', 'as.Smolt',	'As.Adult.ballard')

LW <- read.csv("data/ch_2023.csv") %>%
  filter(Length < 180) %>%
  filter(Year < 2022) %>%
  select(Tagged, as.Smolt,	As.Adult.ballard, Year, ReleaseSite, Length, ReleaseWeek) %>%
  filter(ReleaseSite %in% c("LWBEAR", 'LWCEDR','LWISQH','LWBLRD', 'LWKENM')) %>%
  group_by(!!!syms(locs),Year,ReleaseSite,ReleaseWeek) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(idx = row_number()) %>%
  mutate(y_i = factor(Year, levels = sort(unique(LW$Year)))
  )


MR_settings <- list(
  data = na.omit(data.frame(LW)),
  locs = locs,
  frms = list(
    phi = 'state ~ (1|ReleaseSite) + (1|y_i)', #Survival
    p = 'state ~  (1|ReleaseSite:y_i)', #Detection probability for smolts at the locks.
    lam = 'state ~ (1|ReleaseSite:y_i)'  #joint probability of detection and survival for the  last location
  ),
  make_sd = TRUE
)


tmb_list <- new("tmb_list",
                MR_settings = MR_settings)
tmb_obj <- make_tmb_lists(tmb_list)
tmb_model <- build_tmb(tmb_obj)
AIC(tmb_model)

data.frame(rbind(tmb_model@phi.lm.form$phi.data,tmb_model@lam.lm.form$lam.data) %>% filter(loc != locs[1]) %>% group_by(ReleaseSite,loc, Year) %>% summarize(n = sum(n*state)),
           pred = tmb_model@TMB$sd$value[tmb_model@TMB$sd$value>0],
           sd = tmb_model@TMB$sd$sd[tmb_model@TMB$sd$value>0]) %>%
  ggplot(aes(x = Year, y = n)) +
  geom_point(size = 2) +
  geom_point(aes(x = Year, y = pred, color = ReleaseSite)) +
  geom_errorbar(aes(ymin = pred - 1.96 * sd, ymax = pred + 1.96 * sd, color = ReleaseSite)) +
  facet_grid(loc~ReleaseSite, scales = "free_y") +
  theme_classic()

# tmb_model@TMB$tmb.data$raw %>% group_by(Year,ReleaseSite) %>% summarize(n = sum(n*state))

response_plot(tmb_model,
              process = "p",
              var = 'ReleaseWeek')

p <- plot(obj = tmb_model,
        type = "re",
        process = 'phi',
        re_ef = "ReleaseSite:y_i",
        x = 'y_i',
        color = 'ReleaseSite',
        factor = NULL,
        var = NULL)
# #
p$p + ggplot2::facet_wrap(~ReleaseSite)
#
# p$p + ggplot2::facet_wrap(~ReleaseSite)
#
# p <- plot(obj = tmb_model,
#           type = "re",
#           process = 'lam',
#           # color = 'ReleaseSite',
#           var = NULL)
#
#
#
