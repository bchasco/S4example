re_plot <- function(tmb_out,proc){
  slot_name <- paste0(proc, ".lm.form")  # Create the dynamic slot name
  # Extract the terms object
  f <- slot(tmb_out, slot_name)
  Zt_names <- rownames(f$reTrms$Zt)  # Names from the Zt matrix
  print(Zt_names)
  # Assuming the names in 'Zt_names' are formatted like 'ReleaseSite_fYear'
  # Split the names to create separate 'ReleaseSite' and 'fYear' columns
  split_names <- strsplit(Zt_names, ":")  # Adjust the separator as per your naming scheme
  print(split_names)
  # Create a data.frame
  re_df <- data.frame(do.call(rbind, split_names))  # Convert list of split names to a matrix/data.frame
  colnames(re_df) <- c("y_i")  # Nam
  print(re_df)
  f <- slot(tmb_out, 'TMB')
  print(names(f$rep))
  re_df$est <- f$rep[[paste0(proc,'.re')]]
  print(re_df)
  library(ggplot2)
  p <- ggplot(re_df,aes(x = as.numeric(y_i), y = plogis(2+est)))+ #, color = re_site)) +
    # ylim(50,130)+
    geom_line() +
    ylim(0,1) #+

  if(proc== "p"){
    p <- p +
      ylab("detection probability")
  }
  if(proc== "phi"){
    p <- p +
      ylab("survival probability")
  }
  print(p)
  # scale_fill_viridis_c(option = "magma", direction = 1) +
  #   facet_grid(~ReleaseSite)
}
