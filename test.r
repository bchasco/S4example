library(tidyr)
library(dplyr)

rm(list=ls())
#Build simple data
LW$fYear <- factor(LW$Year, levels = sort(unique(LW$Year)))
#
# # Step 1: Create phi.loc by reshaping and filtering
# phi.loc <- LW %>%
#   mutate(id = row_number()) %>%
#   arrange(id) %>%
#   pivot_longer(cols = c(Tagged, as.Smolt, As.Adult.ballard),
#                names_to = 'loc',
#                values_to = 'state') %>%
#   filter(loc %in% c('as.Smolt'))
#
# # Step 2: Create model matrix and combine with 'id'
# phi.frm <- 'state ~ ReleaseSite + Year'
# phi.mat <- data.frame(model.matrix(as.formula(phi.frm), phi.loc)) %>%
#   mutate(id = phi.loc$id) %>%
#   group_by(id) %>%
#   summarise(y_values = list(pick(everything())))
#
# # Step 1: Create phi.loc by reshaping and filtering
# lam.loc <- LW %>%
#   mutate(id = row_number()) %>%
#   arrange(id) %>%
#   pivot_longer(cols = c(Tagged, as.Smolt, As.Adult.ballard),
#                names_to = 'loc',
#                values_to = 'state') %>%
#   filter(loc %in% c('As.Adult.ballard'))
#
# # Step 2: Create model matrix and combine with 'id'
# lam.frm <- 'state ~ Year'
#
# lam.mat <- data.frame(model.matrix(as.formula(lam.frm), lam.loc)) %>%
#   mutate(id = lam.loc$id) %>%
#   group_by(id) %>%
#   summarise(y_values = list(pick(everything())))

tmb_list <- new("tmb_list",
                raw_data = na.omit(data.frame(LW)), #as.data.frame(do.call(cbind,Orange)),
                phi.frm = 'Tagged ~ 1',
                phi.tmp.frm = 'state ~ ReleaseSite + poly(Length,1) + fYear',
                p.frm = 'Tagged ~ 1',
                lam.frm = 'Tagged ~ 1'
)

tmb_obj <- make_tmb_lists(tmb_list)
tmb_model <- build_tmb(tmb_obj)
marginal_plot(tmb_model,
     proc = 'phi',
     columns = 'Length')

# plot(tmb_model@TMB$rep$lam.re, type="l")

