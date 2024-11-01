x <- read.csv("data/ch_2023.csv")
dim(x)

LW <- x %>%
  select(Tagged, as.Smolt,	As.Adult.ballard, Year, ReleaseSite, Length + ReleaseWeek) %>%
  filter(ReleaseSite %in% c("LWBEAR", 'LWCEDR')) %>%
  group_by(Tagged, as.Smolt,	As.Adult.ballard, Year, ReleaseSite, Length, ReleaseWeek) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(idx = row_number()) #%>%
  # pivot_longer(cols = c('Tagged', 'as.Smolt',	'As.Adult.ballard'), names_to = "loc", values_to = "state")

save(LW,file = "data/LW.rda")
