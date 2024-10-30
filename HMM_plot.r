library(ggplot2)
library(ggforce)
library(dplyr)
library(reshape2)

# Simulate some example data (use your actual data)
# x <- tmb_model@TMB$rep$gamma
x <- reshape2::melt(apply(tmb_model@TMB$rep$gamma,c(2,3,4),mean))

# Process the data for plotting
g <- x %>%
  filter(value >=0, value <1) %>%
  mutate(loc_id = as.numeric(factor(loc, labels = tmb_model@MR_settings$locs))) %>%
  filter(paste(init_state,next_state)!="Dead Dead",
         init_state!="Dead") %>%
  mutate(x1 = loc_id -1,
         x2 = ifelse(next_state == "Dead",loc_id-1,loc_id  - 0.1),
         y1 = ifelse(init_state == "Sub",1.0,2),
         y2 = ifelse(next_state == "Dead",ifelse(init_state=="Yr",2.5,0.5),
                     ifelse(init_state=="Sub",ifelse(next_state=="Yr",1.9,1),2))) %>%
  mutate(lx1 = loc_id -1,
         lx2 = ifelse(next_state == "Dead",loc_id-1.1,ifelse(init_state==next_state,loc_id,loc_id - 0.3)),
         ly1 = ifelse(init_state == "Sub",1.0,2),
         ly2 = ifelse(next_state == "Dead",ifelse(init_state=="Yr",2.65,0.75),
                      ifelse(init_state=="Sub",ifelse(next_state=="Yr",2.15,1.05),2))) %>%
  mutate(rot = ifelse(next_state == "Dead",90,ifelse(next_state == init_state,0,60)),
         angle = atan2(abs(y2 - y1), abs(x2 - x1)) * 180 / pi) %>%
  filter(loc_id>1,
         value !=0) %>%
  mutate(label = paste0(round(value,2)))

# Plot
p <- g %>%
  ggplot(aes(x = x1, y = y1)) +
  xlim(0, 8) +
  scale_y_continuous(
    breaks = c(0.5, 1, 2, 2.5),
    labels = c("Dead", "Alive", "Yearling", "Dead")
  ) +
  scale_x_continuous(
    breaks = seq_along(locs),
    labels = tmb_model@MR_settings$locs
  ) +
  ylab("State") +
  xlab("Detection location") +

  # Add arrows for transitions
  geom_segment(aes(xend = x2, yend = y2),
               size = 0.5,
               linejoin = "mitre",
               arrow = arrow(type = "closed")) +

  # Add dynamic text labels along the arrows
  geom_text(aes(x = (lx1 + lx2) / 2,
                y = (ly1 + ly2) / 2,
                label = label,
                angle   = angle),  # Use dynamic angle for rotation
            size = 4, hjust = 0.5, vjust = -0.5) +

  theme_classic()

# Print and save the plot
print(p)
# ggsave(p, filename = "output/HMM_plot_wo_barge.png", width = 10, height = 5)
