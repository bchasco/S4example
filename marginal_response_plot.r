library(ggplot2)
# Assuming the original poly() transformation is stored in the S4 object

var <- "ReleaseWeek"  # Variable of interest
process <- 'p'  # Process to use (e.g., 'p')
my_factor <- NULL  # Factor variable for grouping
nr <- 20  # Number of rows for the prediction
obj <- tmb_model  # The TMB model object

v <- obj@MR_settings$data[, var]  # Extract the variable data

# Access the formula terms dynamically from the S4 slot
slot_name <- 'MR_settings'  # Slot name for settings
f <- slot(obj, slot_name)  # Extract slot
term_obj <- terms(formula(f$frms[process][[1]]))  # Extract terms

# Extract predictor variables
predictor_terms <- attr(term_obj, "term.labels")
var_form <- predictor_terms[grep(var, predictor_terms)]  # Get formula for the variable
var_transform <- eval(parse(text = var_form), obj@MR_settings$data)  # Apply transformation

# Extract relevant columns for factors in the model
slot_name <- paste0(process,".lm.form")
slot_list <- slot(obj, slot_name)
mod_factors <- names(attr(slot_list[[slot_name]]$X, "contrasts"))
mod_factor_cols <- grep(paste(mod_factors, collapse = "|"), names(slot_list[[paste0(process,'.mat')]][[1]]))

# Get the "assign" attribute
assign_attr <- attr(slot_list[[slot_name]]$X, "assign")
terms_labels <- attr(terms(formula(obj@MR_settings$frms[[process]])), "term.labels")

# Create a prediction data frame based on mean values
pred_df <- matrix(colMeans(slot_list[[slot_name]]$X),
                  nrow = nr,
                  ncol = ncol(slot_list[[slot_name]]$X),
                  byrow = TRUE)

# Set factor columns to 0 and generate a sequence for the variable
pred_df[, mod_factor_cols] <- 0
my_seq <- seq(min(v), max(v), length.out = nr)
pred_df[, grep(var, colnames(slot_list[[slot_name]]$X))] <- predict(var_transform, my_seq)

# Initialize data frames for predictions and standard errors
pred_var <- data.frame(var = my_seq)
pred_sd <- data.frame(var = my_seq)
names(pred_var) <- var
names(pred_sd) <- var

# Initialize counters
icnt <- 2
par_idx <- grep(paste0(process, ".list.b"), names(obj@TMB$opt$par))

# Loop over factors to compute predictions and standard errors
fact_loop <- c(1)
if(!is.null(my_factor))(
  fact_loop <- c(fact_loop,mod_factor_cols)
)
for (i in fact_loop) {
  tmp_df <- pred_df
  tmp_df[, i] <- 1  # Set the factor column to 1 for the current factor

  # Calculate predicted values
  pred_var[, icnt] <- as.data.frame(tmp_df %*% tmb_model@TMB$opt$par[names(tmb_model@TMB$opt$par) %in% paste0(process, '.list.b')])

  # Calculate variances and standard errors
  variances <- diag(tmp_df %*% obj@TMB$sd$cov.fixed[par_idx, par_idx] %*% t(tmp_df))
  standard_error <- sqrt(variances)
  pred_sd[, icnt] <- standard_error
  if(!is.null(my_factor)){
    names(pred_var)[icnt] <- unique(slot_list[[paste0(process,'.data')]][, my_factor])[icnt - 1, 1]
    names(pred_sd)[icnt] <- unique(slot_list[[paste0(process,'.data')]][, my_factor])[icnt - 1, 1]
  }

  icnt <- icnt + 1
}

# Reshape the predictions and standard errors using pivot_longer
pred_var <- pred_var %>%
  pivot_longer(cols = -all_of(var), names_to = "factor", values_to = "est")

pred_sd <- pred_sd %>%
  pivot_longer(cols = -all_of(var), names_to = "factor", values_to = "sd") %>%
  select(sd)
# Combine predictions and standard errors
pred_var <- cbind(pred_var, pred_sd)


# Plot the results with ggplot2
pred_var %>%
  ggplot(aes(x = !!sym(var), y = plogis(est), fill = factor)) +
  geom_line(aes(color = factor), size = 1.2) +
  ylab(process) +
  geom_ribbon(aes(ymin = plogis(est - 1.96 * sd), ymax = plogis(est + 1.96 * sd)), alpha = 0.3) +
  theme_classic()

