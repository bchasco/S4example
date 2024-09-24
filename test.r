rm(list=ls())
#Build simple data
LW$fYear <- factor(LW$Year, levels = sort(unique(LW$Year)))
tmb_list <- new("tmb_list",
                raw_data = na.omit(data.frame(LW)), #as.data.frame(do.call(cbind,Orange)),
                phi.frm = 'Tagged ~ (1|fYear)',
                p.frm = 'Tagged ~ 1',
                lam.frm = 'Tagged ~ (1|fYear)'
)

tmb_obj <- make_tmb_lists(tmb_list)
tmb_model <- build_tmb(tmb_obj)
# marginal_plot(tmb_model,
#      proc = 'phi',
#      columns = 'Length')

tmb_model@TMB$sd
