rm(list=ls())
#Build simple data
data("Orange")

# Orange$Tree <- factor(Orange$Tree, levels = sort(unique(Orange$Tree)))
# Step 2: Create a SurvivalModel object
tmb_list <- new("tmb_list",
                raw_data = iris, #as.data.frame(do.call(cbind,Orange)),
                vars = list(res = 'Sepal.Length',
                            pred_vars = c('Sepal.Width')),
                formula = 'Sepal.Length ~ poly(Sepal.Width,1) + poly(Petal.Width,6) + factor(Species)'
                )

tmb_obj <- make_tmb_lists(tmb_list)
tmb_model <- build_tmb(tmb_obj)
plot(tmb_model,
     columns = 'Petal.Width')
