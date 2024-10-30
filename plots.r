x <- tmb_model
data.frame(x@TMB$tmb.data$raw,est = x@TMB$rep$phi[,2]) %>%
  # filter(ReleaseSite == "LWBEAR") %>%
  filter(y_i %in% 2017:2021) %>%
  group_by(rw_i, ReleaseSite, temp_i, y_i) %>%
  summarize(est = mean(est)) %>%
  ggplot(aes(x = rw_i, y = temp_i, fill = est)) +
  geom_tile() +
  facet_grid(ReleaseSite~y_i) +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  theme_bw()

# data.frame(rbind(tmb_model@phi.lm.form$phi.data,tmb_model@lam.lm.form$lam.data) %>% filter(loc != locs[1]) %>% group_by(ReleaseSite,loc, Year) %>% summarize(n = sum(n*state)),
#            pred = tmb_model@TMB$sd$value[tmb_model@TMB$sd$value>0],
#            sd = tmb_model@TMB$sd$sd[tmb_model@TMB$sd$value>0]) %>%
#   ggplot(aes(x = Year, y = n)) +
#   geom_point(size = 2) +
#   geom_point(aes(x = Year, y = pred, color = ReleaseSite)) +
#   geom_errorbar(aes(ymin = pred - 1.96 * sd, ymax = pred + 1.96 * sd, color = ReleaseSite)) +
#   facet_grid(loc~ReleaseSite, scales = "free_y") +
#   theme_classic() +
#   continuous_
#
#
#
# response_plot(tmb_model,
#               process = "phi",
#               var = "ReleaseWeek",
#               factor = "Species")
# #
# p <- plot(obj = tmb_model,
#         type = "re",
#         process = 'phi',
#         re_ef = "rw_i:temp_i:Species",
#         x = 'temp_i',
#         color = "rw_i",
#         factor = NULL,
#         var = NULL)
# # #
# p$p + ggplot2::facet_grid(~Species)

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
