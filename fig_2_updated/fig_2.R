c("1.5", "", "" ,"6" ,"","","10")


## Recombination of 1e-09
simulation_df <- read_rds("configuration_gravel_simulations_low_recomb_ac_1_2.rds") %>% mutate_if(is.numeric, replace_na,  0) %>%
  mutate(distance=abs(pos_1-pos_2)) %>%
  mutate(variation_type=case_when(
    variation_type == 1 ~ "Nonsynonymous",
    variation_type == 2 ~ "Synonymous"
  )) %>%
  mutate(variation_type=as.factor(variation_type)) %>% filter(recombination_rate==1e-09 )

long_simulation_df_genotype_config <- gather(simulation_df, 
                                 "0|0,0|0", "0|0,0|1", "1|0,0|0", "0|0,1|0", "0|1,0|0", "1|0,0|1", "1|0,1|0", "0|1,1|0", "0|1,0|1", "1|1,0|0" ,"0|0,1|1" ,"1|1,1|1" , 
                                 value="count", 
                                 key="genotype") 


long_simulation_df_genotype_config <- annotate_genotypes(long_simulation_df_genotype_config, type="phased")


annotated_genotypes_recomb_1e9 <- long_simulation_df_genotype_config


## Recombination of 1e-08

simulation_df <- read_rds("configuration_gravel_simulations_average_recomb_ac_1_2.rds") %>% mutate_if(is.numeric, replace_na,  0) %>%
  mutate(distance=abs(pos_1-pos_2)) %>%
  mutate(variation_type=case_when(
    variation_type == 1 ~ "Nonsynonymous",
    variation_type == 2 ~ "Synonymous"
  )) %>%
  mutate(variation_type=as.factor(variation_type)) %>% filter(recombination_rate==1e-08 )

long_simulation_df_genotype_config <- gather(simulation_df, 
                                 "0|0,0|0", "0|0,0|1", "1|0,0|0", "0|0,1|0", "0|1,0|0", "1|0,0|1", "1|0,1|0", "0|1,1|0", "0|1,0|1", "1|1,0|0" ,"0|0,1|1" ,"1|1,1|1" , 
                                 value="count", 
                                 key="genotype") 


long_simulation_df_genotype_config <- annotate_genotypes(long_simulation_df_genotype_config, type="phased")


annotated_genotypes_recomb_1e8 <- long_simulation_df_genotype_config





## Simulations

#Recombination 1e-09
simulations_doubletons_hr_recomb_1e9 <- plot_rh_unphased(annotated_genotypes_recomb_1e9, window_size = 1500, plot_type = "line_plot" , allele_count_to_analyze = 2)[[3]]
simulations_singletons_hr_recomb_1e9 <- plot_rh_unphased(annotated_genotypes_recomb_1e9, window_size = 1500, plot_type = "line_plot" , allele_count_to_analyze = 1)[[3]]
#Recombination 1e-08
simulations_doubletons_hr_recomb_1e8 <- plot_rh_unphased(annotated_genotypes_recomb_1e8, window_size = 1500, plot_type = "line_plot" , allele_count_to_analyze = 2)[[3]]
simulations_singletons_hr_recomb_1e8 <- plot_rh_unphased(annotated_genotypes_recomb_1e8, window_size = 1500, plot_type = "line_plot" , allele_count_to_analyze = 1)[[3]]



# extract a legend that is laid out horizontally
legend_b <- get_legend(
  simulations_doubletons_hr_recomb_1e8 + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

hr_plots <- plot_grid(
simulations_singletons_hr_recomb_1e9 + theme(legend.position = "none") + scale_x_discrete(labels=c("1.5", "", "" ,"6" ,"","","10")) + xlab(NULL)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()), 
simulations_doubletons_hr_recomb_1e9 + theme(legend.position = "none") + scale_x_discrete(labels=c("1.5", "", "" ,"6" ,"","","10")) + xlab(NULL) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()),
simulations_singletons_hr_recomb_1e8 + theme(legend.position = "none") + scale_x_discrete(labels=c("1.5", "", "" ,"6" ,"","","10")), 
simulations_doubletons_hr_recomb_1e8 + theme(legend.position = "none") + scale_x_discrete(labels=c("1.5", "", "" ,"6" ,"","","10")),
labels = c('A', 'B', 'C', 'D'), label_size = 12, nrow=2, ncol=2)



plot_grid(hr_plots,legend_b, ncol=1, rel_heights = c(1, .1))





hr_plots <- plot_grid(
simulations_singletons_hr_recomb_1e9 + theme(legend.position = "none") + scale_x_discrete(labels=c("1.5", "", "" ,"6" ,"","","10"))  + xlab(NULL)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + labs(title="r=1e-09"), 
simulations_singletons_hr_recomb_1e8 + theme(legend.position = "none") + scale_x_discrete(labels=c("1.5", "", "" ,"6" ,"","","10")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + ylab(" ")  + labs(title="r=1e-08"),
simulations_doubletons_hr_recomb_1e9 + theme(legend.position = "none") + scale_x_discrete(labels=c("1.5", "", "" ,"6" ,"","","10")) + xlab("Physical distance (kbp)")  ,
simulations_doubletons_hr_recomb_1e8 + theme(legend.position = "none") + scale_x_discrete(labels=c("1.5", "", "" ,"6" ,"","","10"))+ xlab("Physical distance (kbp)")  + ylab(" "),
labels = c('A', 'C', 'B', 'D'), label_size = 12, nrow=2, ncol=2)



plot_grid(hr_plots,legend_b, ncol=1, rel_heights = c(1, .1))
