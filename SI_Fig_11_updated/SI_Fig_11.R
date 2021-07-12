##### HRJ coorelation for r=1e-09 allele count = 2 
library(tidyverse)
library(cowplot)
theme_set((theme_bw(base_size=22)))

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


simulation_df<-read_rds("/Users/jessegarcia/Documents/SLiM_ParallelRProject copy/data2/configuration_gravel_simulations_low_recomb_ac_1_2.rds") %>% mutate_if(is.numeric, replace_na,  0) %>%
  mutate(distance=abs(pos_1-pos_2)) %>%
  mutate(variation_type=case_when(
    variation_type == 1 ~ "Nonsynonymous",
    variation_type == 2 ~ "Synonymous"
  )) %>%
  mutate(variation_type=as.factor(variation_type)) %>%
  filter(recombination_rate==1e-09, allele_count == 2)

dprime_simulation_df<-simulation_df %>%
  mutate(pAB=(`0|1,0|1` + `1|0,1|0`)/100) %>%
  mutate(pA=(`0|1,0|1` + `0|1,0|0` + `1|0,0|0` + `1|0,1|0` +`1|0,0|1`+`0|1,1|0` +2*`1|1,0|0` + 2*`1|1,1|1`)/100) %>%
  mutate(pB= (`0|0,0|1` + `0|0,1|0` + `0|1,0|1` + `1|0,1|0` + 2*`0|0,1|1` + `1|0,0|1` + `0|1,1|0` + 2*`1|1,1|1`)/100) %>%
  mutate(d=pAB-pA*pB) %>%
  mutate(r_squared = d**2/(pA*(1-pA)*pB*(1-pB))) %>%
  mutate(d_prime = case_when(
    d < 0 ~ d / pmin(pA*pB , (1-pA)*(1-pB)),
    d > 0 ~ d / pmin(pA*(1-pB), (1-pA)*pB )
    ))


long_simulation_df_genotype_config<-gather(simulation_df, 
                                 "0|0,0|0", "0|0,0|1", "1|0,0|0", "0|0,1|0", "0|1,0|0", "1|0,0|1", "1|0,1|0", "0|1,1|0", "0|1,0|1", "1|1,0|0" ,"0|0,1|1" ,"1|1,1|1" , 
                                 value="count", 
                                 key="genotype") 


long_simulation_df_genotype_config<-annotate_genotypes(long_simulation_df_genotype_config, type="phased")


annotated_genotypes<-long_simulation_df_genotype_config

window_size =1500
allele_count_to_analyze = 2


unique_words<-unique(annotated_genotypes$words) 

unique_words<-unique_words[!is.na(unique_words)]


 doubleton_plots<-unique_words %>% map(~ { 
annotated_genotypes %>%
  filter(allele_count == allele_count_to_analyze, !is.na(words), words ==.x)  %>%
      mutate(distance_breaks=cut_interval(distance, length=window_size)) %>%
  group_by(distance_breaks, variation_type, recombination_rate, words) %>%
  summarise(mean_count=mean(count)) %>% 
    ungroup() 
}  )
 
hr_j<-doubleton_plots[[3]]

dprime_simulation_df<-dprime_simulation_df %>%
  filter(allele_count == allele_count_to_analyze)  %>%
      mutate(distance_breaks=cut_interval(distance, length=window_size)) %>%
  group_by(distance_breaks, variation_type, recombination_rate) %>%
  summarise(mean_dprime=mean(d_prime), mean_d=mean(d), mean_rsquare=mean(r_squared)) %>% 
    ungroup() 

cor.test(hr_j$mean_count, dprime_simulation_df$mean_dprime, method = 'spearman')


sim_plot_1<-bind_cols(hr_j,dprime_simulation_df) %>%
  ggplot(aes(x=mean_dprime, y=mean_count, color=variation_type...2)) +
  geom_point(size=2) + 
  labs(y=bquote('Mean '*H[R]^.(allele_count_to_analyze)), color="Variation", subtitle ="r=1e-09", x="Mean D'")+ 
  scale_color_manual(values=c("Synonymous"="dodgerblue1", "Nonsynonymous"="darkorchid2")) +
  theme(text=element_text(size=22))


##### HRJ for correlation r=1e-08 allele count =2

simulation_df<-read_rds("/Users/jessegarcia/Documents/SLiM_ParallelRProject/data2/configuration_gravel_simulations_average_recomb_ac_1_2.rds")  %>% mutate_if(is.numeric, replace_na,  0) %>%
  mutate(distance=abs(pos_1-pos_2)) %>%
  mutate(variation_type=case_when(
    variation_type == 1 ~ "Nonsynonymous",
    variation_type == 2 ~ "Synonymous"
  )) %>%
  mutate(variation_type=as.factor(variation_type)) %>%
  filter(allele_count==2)

dprime_simulation_df<-simulation_df %>%
  mutate(pAB=(`0|1,0|1` + `1|0,1|0`)/100) %>%
  mutate(pA=(`0|1,0|1` + `0|1,0|0` + `1|0,0|0` + `1|0,1|0` +`1|0,0|1`+`0|1,1|0` +2*`1|1,0|0` + 2*`1|1,1|1`)/100) %>%
  mutate(pB= (`0|0,0|1` + `0|0,1|0` + `0|1,0|1` + `1|0,1|0` + 2*`0|0,1|1` + `1|0,0|1` + `0|1,1|0` + 2*`1|1,1|1`)/100) %>%
  mutate(d=pAB-pA*pB) %>%
  mutate(r_squared = d**2/(pA*(1-pA)*pB*(1-pB))) %>%
  mutate(d_prime = case_when(
    d < 0 ~ d / pmin(pA*pB , (1-pA)*(1-pB)),
    d > 0 ~ d / pmin(pA*(1-pB), (1-pA)*pB )
    ))


long_simulation_df_genotype_config<-gather(simulation_df, 
                                 "0|0,0|0", "0|0,0|1", "1|0,0|0", "0|0,1|0", "0|1,0|0", "1|0,0|1", "1|0,1|0", "0|1,1|0", "0|1,0|1", "1|1,0|0" ,"0|0,1|1" ,"1|1,1|1" , 
                                 value="count", 
                                 key="genotype") 


long_simulation_df_genotype_config<-annotate_genotypes(long_simulation_df_genotype_config, type="phased")


annotated_genotypes<-long_simulation_df_genotype_config

window_size =1500
allele_count_to_analyze = 2


unique_words<-unique(annotated_genotypes$words) 

unique_words<-unique_words[!is.na(unique_words)]


 doubleton_plots<-unique_words %>% map(~ { 
annotated_genotypes %>%
  filter(allele_count == allele_count_to_analyze, !is.na(words), words ==.x)  %>%
      mutate(distance_breaks=cut_interval(distance, length=window_size)) %>%
  group_by(distance_breaks, variation_type, recombination_rate, words) %>%
  summarise(mean_count=mean(count)) %>% 
    ungroup() 
}  )
 
hr_j<-doubleton_plots[[3]]

dprime_simulation_df<-dprime_simulation_df %>%
  filter(allele_count == allele_count_to_analyze)  %>%
      mutate(distance_breaks=cut_interval(distance, length=window_size)) %>%
  group_by(distance_breaks, variation_type, recombination_rate) %>%
  summarise(mean_dprime=mean(d_prime), mean_d=mean(d), mean_rsquare=mean(r_squared)) %>% 
    ungroup() 

cor.test(hr_j$mean_count, dprime_simulation_df$mean_dprime, method = 'spearman')


sim_plot_2<-bind_cols(hr_j,dprime_simulation_df) %>%
  ggplot(aes(x=mean_dprime, y=mean_count, color=variation_type...2)) +
  geom_point(size=2) + 
  labs(y=bquote('Mean '*H[R]^.(allele_count_to_analyze)), color="Variation", subtitle ="r=1e-08", x="Mean D'")+ 
  scale_color_manual(values=c("Synonymous"="dodgerblue1", "Nonsynonymous"="darkorchid2"))+
  theme(text=element_text(size=22))





    prow_1<-plot_grid( sim_plot_1+ theme(legend.position="none"),
           sim_plot_2 + theme(legend.position = "none") + labs(y=""),labels=c("B","D"),
           align = 'vh',
           hjust = -1,
           nrow = 1
           )

### HRJ  for correlation r=1e-09 allele count =1
    
simulation_df<-read_rds("/Users/jessegarcia/Documents/SLiM_ParallelRProject copy/data2/configuration_gravel_simulations_low_recomb_ac_1_2.rds") %>% mutate_if(is.numeric, replace_na,  0) %>%
  mutate(distance=abs(pos_1-pos_2)) %>%
  mutate(variation_type=case_when(
    variation_type == 1 ~ "Nonsynonymous",
    variation_type == 2 ~ "Synonymous"
  )) %>%
  mutate(variation_type=as.factor(variation_type)) %>%
  filter(recombination_rate==1e-09, allele_count == 1)

dprime_simulation_df<-simulation_df %>%
  mutate(pAB=(`0|1,0|1` + `1|0,1|0`)/100) %>%
  mutate(pA=(`0|1,0|1` + `0|1,0|0` + `1|0,0|0` + `1|0,1|0` +`1|0,0|1`+`0|1,1|0` +2*`1|1,0|0` + 2*`1|1,1|1`)/100) %>%
  mutate(pB= (`0|0,0|1` + `0|0,1|0` + `0|1,0|1` + `1|0,1|0` + 2*`0|0,1|1` + `1|0,0|1` + `0|1,1|0` + 2*`1|1,1|1`)/100) %>%
  mutate(d=pAB-pA*pB) %>%
  mutate(r_squared = d**2/(pA*(1-pA)*pB*(1-pB))) %>%
  mutate(d_prime = case_when(
    d < 0 ~ d / pmin(pA*pB , (1-pA)*(1-pB)),
    d > 0 ~ d / pmin(pA*(1-pB), (1-pA)*pB )
    ))


long_simulation_df_genotype_config<-gather(simulation_df, 
                                 "0|0,0|0", "0|0,0|1", "1|0,0|0", "0|0,1|0", "0|1,0|0", "1|0,0|1", "1|0,1|0", "0|1,1|0", "0|1,0|1", "1|1,0|0" ,"0|0,1|1" ,"1|1,1|1" , 
                                 value="count", 
                                 key="genotype") 


long_simulation_df_genotype_config<-annotate_genotypes(long_simulation_df_genotype_config, type="phased")


annotated_genotypes<-long_simulation_df_genotype_config

window_size =1500
allele_count_to_analyze = 1


unique_words<-unique(annotated_genotypes$words) 

unique_words<-unique_words[!is.na(unique_words)]


 doubleton_plots<-unique_words %>% map(~ { 
annotated_genotypes %>%
  filter(allele_count == allele_count_to_analyze, !is.na(words), words ==.x)  %>%
      mutate(distance_breaks=cut_interval(distance, length=window_size)) %>%
  group_by(distance_breaks, variation_type, recombination_rate, words) %>%
  summarise(mean_count=mean(count)) %>% 
    ungroup() 
}  )
 
hr_j<-doubleton_plots[[3]]

dprime_simulation_df<-dprime_simulation_df %>%
  filter(allele_count == allele_count_to_analyze)  %>%
      mutate(distance_breaks=cut_interval(distance, length=window_size)) %>%
  group_by(distance_breaks, variation_type, recombination_rate) %>%
  summarise(mean_dprime=mean(d_prime), mean_d=mean(d), mean_rsquare=mean(r_squared)) %>% 
    ungroup() 

cor.test(hr_j$mean_count, dprime_simulation_df$mean_dprime, method = 'spearman')


sim_plot_3<-bind_cols(hr_j,dprime_simulation_df) %>%
  ggplot(aes(x=mean_dprime, y=mean_count, color=variation_type...2)) +
  geom_point(size=2) + 
  labs(y=bquote('Mean '*H[R]^.(allele_count_to_analyze)), color="Variation", subtitle ="r=1e-09", x="Mean D'")+ 
  scale_color_manual(values=c("Synonymous"="dodgerblue1", "Nonsynonymous"="darkorchid2")) +
  theme(text=element_text(size=22))


##### HRJ for correlation r=1e-08 allele count =1

simulation_df<-read_rds("/Users/jessegarcia/Documents/SLiM_ParallelRProject/data2/configuration_gravel_simulations_average_recomb_ac_1_2.rds")  %>% mutate_if(is.numeric, replace_na,  0) %>%
  mutate(distance=abs(pos_1-pos_2)) %>%
  mutate(variation_type=case_when(
    variation_type == 1 ~ "Nonsynonymous",
    variation_type == 2 ~ "Synonymous"
  )) %>%
  mutate(variation_type=as.factor(variation_type)) %>%
  filter(allele_count==1)

dprime_simulation_df<-simulation_df %>%
  mutate(pAB=(`0|1,0|1` + `1|0,1|0`)/100) %>%
  mutate(pA=(`0|1,0|1` + `0|1,0|0` + `1|0,0|0` + `1|0,1|0` +`1|0,0|1`+`0|1,1|0` +2*`1|1,0|0` + 2*`1|1,1|1`)/100) %>%
  mutate(pB= (`0|0,0|1` + `0|0,1|0` + `0|1,0|1` + `1|0,1|0` + 2*`0|0,1|1` + `1|0,0|1` + `0|1,1|0` + 2*`1|1,1|1`)/100) %>%
  mutate(d=pAB-pA*pB) %>%
  mutate(r_squared = d**2/(pA*(1-pA)*pB*(1-pB))) %>%
  mutate(d_prime = case_when(
    d < 0 ~ d / pmin(pA*pB , (1-pA)*(1-pB)),
    d > 0 ~ d / pmin(pA*(1-pB), (1-pA)*pB )
    ))


long_simulation_df_genotype_config<-gather(simulation_df, 
                                 "0|0,0|0", "0|0,0|1", "1|0,0|0", "0|0,1|0", "0|1,0|0", "1|0,0|1", "1|0,1|0", "0|1,1|0", "0|1,0|1", "1|1,0|0" ,"0|0,1|1" ,"1|1,1|1" , 
                                 value="count", 
                                 key="genotype") 


long_simulation_df_genotype_config<-annotate_genotypes(long_simulation_df_genotype_config, type="phased")


annotated_genotypes<-long_simulation_df_genotype_config

window_size =1500
allele_count_to_analyze = 1


unique_words<-unique(annotated_genotypes$words) 

unique_words<-unique_words[!is.na(unique_words)]


 doubleton_plots<-unique_words %>% map(~ { 
annotated_genotypes %>%
  filter(allele_count == allele_count_to_analyze, !is.na(words), words ==.x)  %>%
      mutate(distance_breaks=cut_interval(distance, length=window_size)) %>%
  group_by(distance_breaks, variation_type, recombination_rate, words) %>%
  summarise(mean_count=mean(count)) %>% 
    ungroup() 
}  )
 
hr_j<-doubleton_plots[[3]]

dprime_simulation_df<-dprime_simulation_df %>%
  filter(allele_count == allele_count_to_analyze)  %>%
      mutate(distance_breaks=cut_interval(distance, length=window_size)) %>%
  group_by(distance_breaks, variation_type, recombination_rate) %>%
  summarise(mean_dprime=mean(d_prime), mean_d=mean(d), mean_rsquare=mean(r_squared)) %>% 
    ungroup() 

cor.test(hr_j$mean_count, dprime_simulation_df$mean_dprime, method = 'spearman')


sim_plot_4<-bind_cols(hr_j,dprime_simulation_df) %>%
  ggplot(aes(x=mean_dprime, y=mean_count, color=variation_type...2)) +
  geom_point(size=2) + 
  labs(y=bquote('Mean '*H[R]^.(allele_count_to_analyze)), color="Variation", subtitle ="r=1e-08", x="Mean D'")+ 
  scale_color_manual(values=c("Synonymous"="dodgerblue1", "Nonsynonymous"="darkorchid2"))+
  theme(text=element_text(size=22))





    prow_2<-plot_grid( sim_plot_3+ theme(legend.position="none"),
           sim_plot_4 + theme(legend.position = "none") + labs(y=""), labels=c("A", "C"),
           align = 'vh',
           hjust = -1,
           nrow = 1
           )

    
    
    

# extract the legend from one of the plots
# (clearly the whole thing only makes sense if all plots
# have the same legend, so we can arbitrarily pick one.)
legend <- get_legend(sim_plot_1 + theme(legend.position = "bottom"))

# add the legend to the row we made earlier. Give it one-third of the width
# of one plot (via rel_widths).
si_fig_11 <- plot_grid( prow_2,prow_1 , legend, ncol = 1, rel_heights = c(1,1, .1))
si_fig_11


ggsave(filename="figures/si_figure_11_correlations_dprime_hrj.tiff", plot=si_fig_11, width=20, height=12)
