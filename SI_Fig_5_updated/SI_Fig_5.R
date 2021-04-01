library(data.table)
library(tidyverse)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

colors=c("0"="#0072B2", "-1e-04"="#F0E442", "-0.001"="#009E73", "-0.01"= "#56B4E9", "-0.1"="#E69F00", "DFE"= "#999999")

selection_levels=c("0", "-1e-04", "-0.001", "-0.01", "-0.1", "DFE")
selection_levels=rev(selection_levels)


## Changing the order of dominance levels for this
recombination_levels=c("r=1e-06", "r=1e-07", "r=1e-08" ,"r=1e-09")
dominance_levels=c("h=0.5", "h=0" )

library(data.table)
library(tidyverse)
library(glue)

## This file can be found on the paper's data dryad dataset.
doubletons_ld<-as_tibble(fread("/Users/jessegarcia/Documents/SLiM_ParallelRProject copy/data2/constant_population_constant_selection_constant_dfe_doubletons_all_replicates_ld.csv", data.table = FALSE))

set.seed(1)
subsampled_doubletons<-doubletons_ld %>%
  filter(selection_coefficient != -0.1) %>%
  group_by(recombination_rate, selection_coefficient, dominance_coefficient, seed) %>%
  nest() %>%
  ungroup() %>%
  group_by(recombination_rate, selection_coefficient, dominance_coefficient) %>%
  sample_n(150) %>%
  mutate(group=rep(c(1,2,3), each=50)) %>%
  unnest(data)


subsampled_doubletons_summary<-subsampled_doubletons %>%
  group_by(recombination_rate, selection_coefficient, dominance_coefficient,group) %>%
  summarise(mean_r_2=mean(r_2)) %>% 
  ungroup() %>%
  mutate(selection_coefficient = replace(selection_coefficient, selection_coefficient == -999, "DFE")) %>%
  mutate(selection_coefficient = factor(selection_coefficient, levels=selection_levels) ) %>%
  mutate(dominance_coefficient = glue("h={dominance_coefficient}")) %>%
  mutate(recombination_rate=glue("r={recombination_rate}")) %>%
  mutate(dominance_coefficient=factor(dominance_coefficient,levels= dominance_levels)) %>%
  mutate(recombination_rate=factor(recombination_rate, levels=recombination_levels)) 


subsampled_doubletons_summary_error_bars<-subsampled_doubletons_summary %>%
  group_by(recombination_rate,selection_coefficient,dominance_coefficient) %>%
  summarise(max_r_2=max(mean_r_2), min_r_2=min(mean_r_2), middle_r_2= sort(mean_r_2, decreasing = FALSE)[2] , mean_of_mean_r2=mean(mean_r_2)) %>% 
  ungroup()





si_fig_5 <- subsampled_doubletons_summary_error_bars %>%
  mutate(recombination_rate=fct_rev(recombination_rate)) %>%
  ggplot(aes(x=recombination_rate, y=mean_of_mean_r2, colour=selection_coefficient)) +
  geom_errorbar(aes(x=recombination_rate,ymin = min_r_2, ymax=max_r_2) , size=2 , width=.1,alpha=0.5)  +
  facet_wrap(~dominance_coefficient, nrow = 2) +
  geom_line(aes(group=selection_coefficient)) +
  scale_colour_manual(values=colors) + 
  guides(fill=FALSE)  +
  theme_bw() +
    theme(strip.background =element_rect(fill="white"), text = element_text(size = 18)) +
  labs(x="Recombination rate",y=bquote("Mean"~r^2 ), colour="Selection Coefficient") +
  geom_line(aes(group=selection_coefficient), size=2)

si_fig_5
ggsave(filename="../figures/si_figure_5_mean_r2.tiff", plot=si_fig_5, width=20, height=12)
