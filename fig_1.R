constant_ld<-as_tibble(fread("/Users/jessegarcia/Documents/SLiM_ParallelRProject/data2/constant_dfe_and_constant_selection_ld.csv", data.table = FALSE))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


selection_levels=c("0", "-1e-04", "-0.001", "-0.01", "-0.1", "DFE")
selection_levels=levels=rev(selection_levels)

recombination_levels=c("r=1e-06", "r=1e-07", "r=1e-08" ,"r=1e-09")
dominance_levels=c("h=0",  "h=0.5")



cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


selection_levels=c("0", "-1e-04", "-0.001", "-0.01", "-0.1", "DFE")
selection_levels=rev(selection_levels)


## Changing the order of dominance levels for this
recombination_levels=c("r=1e-06", "r=1e-07", "r=1e-08" ,"r=1e-09")
dominance_levels=c("h=0.5", "h=0" )



summarized_constant_ld=constant_ld %>%
  mutate(distance_breaks=cut_interval(x=dist, length=1000, labels=FALSE)) %>%
  group_by(recombination_rate, selection_coefficient, dominance_coefficient,seed, distance_breaks) %>%
  summarise(mean_r_2=mean(r_2)) %>% 
  ungroup() %>%
  mutate(selection_coefficient = replace(selection_coefficient, selection_coefficient == -999, "DFE")) %>%
  mutate(selection_coefficient = factor(selection_coefficient, levels=levels) ) %>%
  mutate(dominance_coefficient = glue("h={dominance_coefficient}")) %>%
  mutate(recombination_rate=glue("r={recombination_rate}")) %>%
  mutate(dominance_coefficient=factor(dominance_coefficient,levels= dominance_levels)) %>%
  mutate(recombination_rate=factor(recombination_rate, levels=recombination_levels)) %>%
  group_by(recombination_rate, selection_coefficient, dominance_coefficient, distance_breaks) %>%
  mutate(min_r2=min(mean_r_2), max_r2=max(mean_r_2)) %>%
  ungroup()




  
yaxis_set_zero_plot<-summarized_constant_ld %>%
  mutate(recombination_rate=fct_rev(recombination_rate)) %>%
  ggplot(aes(x=distance_breaks, y=mean_r_2, color=factor(selection_coefficient)) )+  
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), aes(fill=selection_coefficient)) +
  theme_bw() +
  labs(x="Distance (kbp)", y=bquote("Mean"~r^2 ), colour="Selection Coefficient" ) +
  facet_wrap(dominance_coefficient ~ recombination_rate, scales="free_y", nrow=2) +
  scale_colour_manual(values=cbPalette) + 
  scale_fill_manual(values=cbPalette) +
  theme(strip.background =element_rect(fill="white"), text = element_text(size = 22))  +
  scale_x_continuous(breaks = scales::pretty_breaks(n=3)) + 
  guides(fill=FALSE)  + scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))


ggsave(filename="fig_1_yaxis_0.png", plot=yaxis_set_zero_plot, width=20, height=12)
