
library(tidyverse)
library(cowplot)

annotate_genotypes<-function(long_genotype_config, type){
  
  
  if (type == "unphased"){
  long_genotype_config<-long_genotype_config %>%
  mutate(words=case_when(
    genotype == "0/0,0/0" ~ "Homozygous reference (0/0,0/0)", 
    genotype == "0/0,0/1" ~ "Single heterozygote (0/1,0/0 or 0/0,0/1)", 
    genotype == "0/1,0/0" ~ "Single heterozygote (0/1,0/0 or 0/0,0/1)",
    genotype == "0/1,0/1" ~ "Double Heterozygote (0/1,0/1)"
  ))
  
  return(long_genotype_config)
  }
  
    
  if (type == "phased"){
  long_genotype_config<-long_genotype_config %>%
  mutate(words=case_when(
    genotype == "0|0,0|0" ~ "Homozygous reference (0|0,0|0)", 
    genotype == "0|0,0|1" ~ "Single heterozygote (0|0,0|1 & 1|0,0|0 & 0|1,0|0 & 0|0,1|0)",
    genotype == "1|0,0|0" ~ "Single heterozygote (0|0,0|1 & 1|0,0|0 & 0|1,0|0 & 0|0,1|0)",
    genotype == "0|1,0|0" ~ "Single heterozygote (0|0,0|1 & 1|0,0|0 & 0|1,0|0 & 0|0,1|0)",
    genotype == "0|0,1|0" ~ "Single heterozygote (0|0,0|1 & 1|0,0|0 & 0|1,0|0 & 0|0,1|0)",
    genotype == "0|1,0|1" ~ "Double heterozygote (0|1,0|1 & 0|1,1|0 & 1|0,1|0 & 1|0,0|1)",
    genotype == "0|1,1|0" ~ "Double heterozygote (0|1,0|1 & 0|1,1|0 & 1|0,1|0 & 1|0,0|1)",
    genotype == "1|0,1|0" ~ "Double heterozygote (0|1,0|1 & 0|1,1|0 & 1|0,1|0 & 1|0,0|1)",
    genotype == "1|0,0|1" ~ "Double heterozygote (0|1,0|1 & 0|1,1|0 & 1|0,1|0 & 1|0,0|1)"
    
  ))
  
  return(long_genotype_config)
  }
  
  
    if (type == "phased_coupling"){
  long_genotype_config<-long_genotype_config %>%
  mutate(words=case_when(
    genotype == "0|0,0|0" ~ "Homozygous reference (0|0,0|0)", 
    genotype == "0|0,0|1" ~ "Single heterozygote (0|0,0|1 & 1|0,0|0 & 0|1,0|0 & 0|0,1|0)",
    genotype == "1|0,0|0" ~ "Single heterozygote (0|0,0|1 & 1|0,0|0 & 0|1,0|0 & 0|0,1|0)",
    genotype == "0|1,0|0" ~ "Single heterozygote (0|0,0|1 & 1|0,0|0 & 0|1,0|0 & 0|0,1|0)",
    genotype == "0|0,1|0" ~ "Single heterozygote (0|0,0|1 & 1|0,0|0 & 0|1,0|0 & 0|0,1|0)",
    genotype == "0|1,0|1" ~ "Double heterozygote Coupling (0|1,0|1 & 1|0,1|0)",
    genotype == "0|1,1|0" ~ "Double heterozygote Repulsion (0|1,1|0 & 1|0,0|1)",
    genotype == "1|0,1|0" ~ "Double heterozygote Coupling (0|1,0|1 & 1|0,1|0)",
    genotype == "1|0,0|1" ~ "Double heterozygote Repulsion (0|1,1|0 & 1|0,0|1)"
    
  ))
  
  return(long_genotype_config)
  }
  
  
  
  
}




plot_rh_unphased<-function(annotated_genotypes, window_size, plot_type, allele_count_to_analyze){
  unique_words<-unique(annotated_genotypes$words) 

unique_words<-unique_words[!is.na(unique_words)]


testthat::expect_true(unique_words[3] == "Double Heterozygote (0/1,0/1)" | unique_words[3] == "Double heterozygote (0|1,0|1 & 0|1,1|0 & 1|0,1|0 & 1|0,0|1)")
if (plot_type == "bar_plot"){
doubleton_plots<-unique_words %>% map(~ { 
annotated_genotypes %>%
  filter(allele_count == allele_count_to_analyze, !is.na(words), words ==.x)  %>%
      mutate(distance_breaks=cut_interval(distance, length=window_size)) %>%
  group_by(distance_breaks, variation_type, recombination_rate, words) %>%
  summarise(mean_count=mean(count)) %>% 
    ungroup() %>%
  ggplot(aes(x=distance_breaks, y=mean_count, fill=variation_type)) +
  geom_col(position="dodge2") +
  facet_grid(~words,scales="free") + 
    theme_bw()
} ) } else if (plot_type == "line_plot") {
  
  doubleton_plots<-unique_words %>% map(~ { 
annotated_genotypes %>%
  filter(allele_count == allele_count_to_analyze, !is.na(words), words ==.x)  %>%
      mutate(distance_breaks=cut_interval(distance, length=window_size)) %>%
  group_by(distance_breaks, variation_type, recombination_rate, words) %>%
  summarise(mean.count=mean(count),
            sd.count = sd(count, na.rm = TRUE),
    n.count = n()) %>%
mutate(se.count = sd.count / sqrt(n.count),
 lower.ci.count = mean.count - se.count,
 upper.ci.count = mean.count + se.count) %>%
          ungroup() %>%
  ggplot(aes(x=distance_breaks, y=mean.count, colour=variation_type, group=variation_type)) +
    geom_errorbar(aes(ymin=lower.ci.count, ymax=upper.ci.count), width=.2) +
  ylab(bquote('Mean '*H[R]^.(allele_count_to_analyze))) +
  labs( x="Physical distance (kbp)", colour="Variation") +
  geom_point(size=2) +
  geom_line(size=1.5) +
#  facet_grid(~words,scales="free") + 
  theme_bw() +
  theme(text=element_text(size=16)) +
  scale_color_manual(values=c("Synonymous"="dodgerblue1", "Nonsynonymous"="darkorchid2"))
}  )
}


  
return(doubleton_plots)
}



## Recombination of 1e-09
# path on computer "/Users/jessegarcia/Documents/SLiM_ParallelRProject copy/data2/configuration_gravel_simulations_low_recomb_ac_1_2.rds"
simulation_df <- read_rds("configuration_gravel_simulations_low_recomb_ac_1_2.rds") %>% mutate_if(is.numeric, replace_na,  0) %>%
  mutate(distance=abs(pos_1-pos_2)) %>%
  mutate(variation_type=case_when(
    variation_type == 1 ~ "Nonsynonymous",
    variation_type == 2 ~ "Synonymous"
  )) %>%
  mutate(variation_type=as.factor(variation_type)) %>% filter(recombination_rate==1e-09 )

revisions_simulation_df <- read_rds("configuration_gravel_simulations_second_round_revisions_low_recomb_ac_1_2.rds") %>% mutate_if(is.numeric, replace_na,  0) %>%
  mutate(distance=abs(pos_1-pos_2)) %>%
  mutate(variation_type=case_when(
    variation_type == 1 ~ "Nonsynonymous",
    variation_type == 2 ~ "Synonymous"
  )) %>%
  mutate(variation_type=as.factor(variation_type)) %>% filter(recombination_rate==1e-09 )

revisions_simulation_second_df <- read_rds("configuration_gravel_simulations_second_round_revisions_second_simulations_low_recomb_ac_1_2.rds") %>% mutate_if(is.numeric, replace_na,  0) %>%
  mutate(distance=abs(pos_1-pos_2)) %>%
  mutate(variation_type=case_when(
    variation_type == 1 ~ "Nonsynonymous",
    variation_type == 2 ~ "Synonymous"
  )) %>%
  mutate(variation_type=as.factor(variation_type)) %>% filter(recombination_rate==1e-09 )

simulation_df <- bind_rows(simulation_df,revisions_simulation_df ,revisions_simulation_second_df)
#simulation_df<- revisions_simulation_df
long_simulation_df_genotype_config <- gather(simulation_df, 
                                 "0|0,0|0", "0|0,0|1", "1|0,0|0", "0|0,1|0", "0|1,0|0", "1|0,0|1", "1|0,1|0", "0|1,1|0", "0|1,0|1", "1|1,0|0" ,"0|0,1|1" ,"1|1,1|1" , 
                                 value="count", 
                                 key="genotype") 


long_simulation_df_genotype_config <- annotate_genotypes(long_simulation_df_genotype_config, type="phased")


annotated_genotypes_recomb_1e9 <- long_simulation_df_genotype_config


## Recombination of 1e-08
# path on computer "/Users/jessegarcia/Desktop/SLiM_ParallelRProject/data2/configuration_gravel_simulations_average_recomb_ac_1_2.rds"
simulation_df <- read_rds("configuration_gravel_simulations_average_recomb_ac_1_2.rds") %>% mutate_if(is.numeric, replace_na,  0) %>%
  mutate(distance=abs(pos_1-pos_2)) %>%
  mutate(variation_type=case_when(
    variation_type == 1 ~ "Nonsynonymous",
    variation_type == 2 ~ "Synonymous"
  )) %>%
  mutate(variation_type=as.factor(variation_type)) %>% filter(recombination_rate==1e-08 )

revisions_simulation_df <- read_rds("configuration_gravel_simulations_second_round_revisions_average_recomb_ac_1_2.rds") %>% mutate_if(is.numeric, replace_na,  0) %>%
  mutate(distance=abs(pos_1-pos_2)) %>%
  mutate(variation_type=case_when(
    variation_type == 1 ~ "Nonsynonymous",
    variation_type == 2 ~ "Synonymous"
  )) %>%
  mutate(variation_type=as.factor(variation_type)) %>% filter(recombination_rate==1e-08 )

revisions_simulation_second_df <- read_rds("configuration_gravel_simulations_second_round_revisions_second_simulations_average_recomb_ac_1_2.rds") %>% mutate_if(is.numeric, replace_na,  0) %>%
  mutate(distance=abs(pos_1-pos_2)) %>%
  mutate(variation_type=case_when(
    variation_type == 1 ~ "Nonsynonymous",
    variation_type == 2 ~ "Synonymous"
  )) %>%
  mutate(variation_type=as.factor(variation_type)) %>% filter(recombination_rate==1e-08 )

simulation_df <- bind_rows(simulation_df,revisions_simulation_df,revisions_simulation_second_df )
#simulation_df <- revisions_simulation_df
long_simulation_df_genotype_config <- gather(simulation_df, 
                                 "0|0,0|0", "0|0,0|1", "1|0,0|0", "0|0,1|0", "0|1,0|0", "1|0,0|1", "1|0,1|0", "0|1,1|0", "0|1,0|1", "1|1,0|0" ,"0|0,1|1" ,"1|1,1|1" , 
                                 value="count", 
                                 key="genotype") 


long_simulation_df_genotype_config <- annotate_genotypes(long_simulation_df_genotype_config, type="phased")


annotated_genotypes_recomb_1e8 <- long_simulation_df_genotype_config





## Simulations

#Recombination 1e-09
simulations_doubletons_hr_recomb_1e9 <- plot_rh_unphased(annotated_genotypes_recomb_1e9, window_size = 1500, plot_type = "line_plot" , allele_count_to_analyze = 2)[[3]] + theme(text = element_text(size = 20)) 
simulations_singletons_hr_recomb_1e9 <- plot_rh_unphased(annotated_genotypes_recomb_1e9, window_size = 1500, plot_type = "line_plot" , allele_count_to_analyze = 1)[[3]] + theme(text = element_text(size = 20)) 
#Recombination 1e-08
simulations_doubletons_hr_recomb_1e8 <- plot_rh_unphased(annotated_genotypes_recomb_1e8, window_size = 1500, plot_type = "line_plot" , allele_count_to_analyze = 2)[[3]] + theme(text = element_text(size = 20)) 
simulations_singletons_hr_recomb_1e8 <- plot_rh_unphased(annotated_genotypes_recomb_1e8, window_size = 1500, plot_type = "line_plot" , allele_count_to_analyze = 1)[[3]] + theme(text = element_text(size = 20)) 



# extract a legend that is laid out horizontally
legend_b <- get_legend(
  simulations_doubletons_hr_recomb_1e8 + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)




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



simulated_hr_plots<-plot_grid(hr_plots,legend_b, ncol=1, rel_heights = c(1, .1))
simulated_hr_plots
ggsave(filename = "figure_2_simulated_hr_plots.png", width=20, height=12)
