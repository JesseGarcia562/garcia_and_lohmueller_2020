selection_levels <- c("0", "-1e-04", "-0.001", "-0.01", "-0.1", "DFE")
selection_levels <- rev(selection_levels)

recombination_levels <- c("r=1e-06", "r=1e-07", "r=1e-08" ,"r=1e-09")
dominance_levels <- c("h=0.5",  "h=0")

doubletons_ld <- as_tibble(fread("constant_population_constant_selection_constant_dfe_doubletons_all_replicates_ld.csv", data.table = FALSE)) %>%
  mutate(d=case_when(
    r_2 == 4e-04 ~ -0.0004 , 
    r_2 == 0.2399 ~ 0.0096, 
    r_2 == 1 ~ 0.0196
  ))

colors <- c("0"="#0072B2", "-1e-04"="#F0E442", "-0.001"="#009E73", "-0.01"= "#56B4E9", "-0.1"="#E69F00", "DFE"= "#999999")


set.seed(1)
subsampled_doubletons <- doubletons_ld %>%
  filter(selection_coefficient != -0.1) %>%
  group_by(recombination_rate, selection_coefficient, dominance_coefficient, seed) %>%
  nest() %>%
  ungroup() %>%
  group_by(recombination_rate, selection_coefficient, dominance_coefficient) %>%
  sample_n(150) %>%
  mutate(group=rep(c(1,2,3), each=50)) %>%
  unnest(data) %>% 
  ungroup()


subsampled_doubletons_summary <- subsampled_doubletons %>%
  group_by(recombination_rate, selection_coefficient, dominance_coefficient,group) %>%
  summarise(mean_d=mean(d))



subsampled_doubletons_summary <- subsampled_doubletons_summary %>%
  group_by(recombination_rate, selection_coefficient, dominance_coefficient) %>%
summarise(max_d=max(mean_d), min_d=min(mean_d), middle_d=sort(mean_d, decreasing = FALSE)[2] , mean_of_mean_d=mean(mean_d)) %>% 
  ungroup() %>%
  mutate(selection_coefficient = replace(selection_coefficient, selection_coefficient == -999, "DFE")) %>%
  mutate(selection_coefficient = factor(selection_coefficient, levels=levels) ) %>%
  mutate(dominance_coefficient = glue("h={dominance_coefficient}")) %>%
  mutate(recombination_rate=glue("r={recombination_rate}")) %>%
  mutate(dominance_coefficient=factor(dominance_coefficient,levels= dominance_levels)) %>%
  mutate(recombination_rate=factor(recombination_rate, levels=recombination_levels)) 


dprime_mean_plot <- subsampled_doubletons_summary %>%
  mutate(recombination_rate = fct_rev(recombination_rate)) %>%
  filter(selection_coefficient != -0.1) %>%
  ggplot(aes(x=recombination_rate, y=mean_of_mean_d_prime, colour=selection_coefficient)) +
  geom_point() +
    geom_errorbar(aes(x=recombination_rate,ymin = min_d_prime, ymax=max_d_prime) , size=3 , width=.2,alpha=0.5)  +
  facet_wrap(~dominance_coefficient , nrow=2)  +
  scale_colour_manual(values=colors) + 
  guides(fill=FALSE)  +
  theme_bw() +
    theme(strip.background =element_rect(fill="white"), text = element_text(size = 22)) +
  labs(x="Recombination rate",y=bquote("Mean D'" ), colour="Selection Coefficient") +
  geom_line(aes(group=interaction(selection_coefficient)), size=2, alpha=1)

dprime_mean_plot

