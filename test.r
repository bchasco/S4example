rm(list=ls())
#Build simple data
data("CO2")

# Orange$Tree <- factor(Orange$Tree, levels = sort(unique(Orange$Tree)))
# Step 2: Create a SurvivalModel object
iris$fSpecies <- factor(iris$Species)
tmb_list <- new("tmb_list",
                raw_data = na.omit(LW), #as.data.frame(do.call(cbind,Orange)),
                formula = 'Tagged ~ (1|Year)'
                )

tmb_obj <- make_tmb_lists(tmb_list)
tmb_model <- build_tmb(tmb_obj)
plot(tmb_model,
     columns = 'Solar.R')
