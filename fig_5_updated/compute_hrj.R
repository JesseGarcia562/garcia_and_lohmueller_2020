library(tidyverse)

high_coverage_hg38_genotype_configuration <- read_csv("/Users/jessegarcia/Documents/SLiM_ParallelRProject copy/data2/high_coverage_hg38_genotype_configuration.csv")

long_low_genotype_config <- read_csv("long_low_genotype_config.csv")


## Configuration of doubletons in r=1e-09 Gravel Simulations
simulation_df<-read_rds("configuration_gravel_simulations_low_recomb_ac_1_2.rds") %>% mutate_if(is.numeric, replace_na,  0) %>%
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


annotated_genotypes <- long_simulation_df_genotype_config


## We need to put all the data here as well. 

## Simulations
simulations_doubletons_hr <- plot_rh_unphased(annotated_genotypes, window_size = 1500, plot_type = "line_plot" , allele_count_to_analyze = 2)

## High coverage hg38

highcoverage_hg38_doubletons_hr <- plot_rh_unphased(high_coverage_hg38_genotype_configuration %>% mutate(recombination_rate="unknown")%>% 
                                                      mutate(variation_type=case_when(
                                                        variation_type == "nonsynonymous_SNV" ~ "Nonsynonymous", 
                                                        variation_type == "synonymous_SNV" ~ "Synonymous" 
                                                      )) , window_size = 1500, plot_type = "line_plot" , allele_count_to_analyze = 2)

## low coverage hg19
low_coverage_doubletons_hr <- plot_rh_unphased(long_low_genotype_config %>% mutate(recombination_rate="unknown")%>% 
                                                 mutate(variation_type=case_when(
                                                   variation_type == "nonsynonymous_SNV" ~ "Nonsynonymous", 
                                                   variation_type == "synonymous_SNV" ~ "Synonymous" 
                                                 )) , window_size = 1500, plot_type = "line_plot" , allele_count_to_analyze = 2)

# extract a legend that is laid out horizontally
legend_b <- get_legend(
  low_coverage_doubletons_hr[[3]]+ 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

hr_plots <- plot_grid(
  simulations_doubletons_hr[[3]] + theme(legend.position = "none") + scale_x_discrete(labels=c("[0,1.5e+03]", "", "" ,"(4.5e+03,6e+03]" ,"","","(9e+03,1.0e+04]")) + theme(axis.text.x = element_text(size=14, angle =45 , hjust=1 ) ) + xlab(NULL), 
  low_coverage_doubletons_hr[[3]] + theme(legend.position = "none") + ylab(NULL)+ scale_x_discrete(labels=c("[0,1.5e+03]", "", "" ,"(4.5e+03,6e+03]" ,"","","(9e+03,1.0e+04]"))  + theme(axis.text.x = element_text(size=14 , angle =45, hjust=1)  ),
  highcoverage_hg38_doubletons_hr[[3]] + theme(legend.position = "none") + ylab(NULL)+ scale_x_discrete(labels=c("[0,1.5e+03]", "", "" ,"(4.5e+03,6e+03]" ,"","","(9e+03,1.0e+04]"))+ theme(axis.text.x = element_text(size=14, angle =45, hjust=1)) + xlab(NULL), 
  labels = c('A', 'B', 'C'), label_size = 12, nrow=1, ncol=3)

### just empirical

## High coverage hg38

highcoverage_hg38_singletons_hr <- plot_rh_unphased(high_coverage_hg38_genotype_configuration %>% mutate(recombination_rate="unknown")%>% 
                                                      mutate(variation_type=case_when(
                                                        variation_type == "nonsynonymous_SNV" ~ "Nonsynonymous", 
                                                        variation_type == "synonymous_SNV" ~ "Synonymous" 
                                                      )) , window_size = 1500, plot_type = "line_plot" , allele_count_to_analyze = 1)[[3]]



## low coverage hg19
low_coverage_singletons_hr <- plot_rh_unphased(long_low_genotype_config %>% mutate(recombination_rate="unknown")%>% 
                                                 mutate(variation_type=case_when(
                                                   variation_type == "nonsynonymous_SNV" ~ "Nonsynonymous", 
                                                   variation_type == "synonymous_SNV" ~ "Synonymous" 
                                                 )) , window_size = 1500, plot_type = "line_plot" , allele_count_to_analyze = 1)[[3]]


hr_plots <- plot_grid(low_coverage_doubletons_hr[[3]]+ xlab(NULL)  + theme(legend.position = "none") + 
                        theme(axis.title.x=element_blank(),
                              axis.text.x=element_blank(),
                              axis.ticks.x=element_blank(), text=element_text(size=22)) , 
                      low_coverage_singletons_hr + scale_x_discrete(labels=c("1.5", "", "" ,"6" ,"","","10")) + theme(legend.position = "none", text=element_text(size=22)), 
                      ncol=1,
                      labels=c("A", "B"),
                      label_size = 14,
                      align="hv") 


plot_grid(hr_plots, legend_b, ncol=1,  rel_heights = c(1, .09))
©