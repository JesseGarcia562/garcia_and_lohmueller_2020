library(tidyverse)

##### HRJ coorelation for r=1e-09 allele count = 2 

theme_set((theme_bw(base_size=22)))

simulation_df<-read_rds("configuration_gravel_simulations_low_recomb_ac_1_2.rds") %>% mutate_if(is.numeric, replace_na,  0) %>%
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

simulation_df<-read_rds("configuration_gravel_simulations_average_recomb_ac_1_2.rds")  %>% mutate_if(is.numeric, replace_na,  0) %>%
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





low_coverage_hg19<-"cofiguration_low_coverage_ac_1_2.rds"

long_low_genotype_config <- read_csv("long_low_genotype_config.csv") %>% mutate(recombination_rate="unknown")%>% 
  mutate(variation_type=case_when(
  variation_type == "nonsynonymous_SNV" ~ "Nonsynonymous", 
  variation_type == "synonymous_SNV" ~ "Synonymous" 
))%>%
  filter(allele_count==2)



low_coverage_hg19<-read_rds(low_coverage_hg19) %>% mutate_if(is.numeric, replace_na, 0) %>% mutate(data_type="Low Coverage (hg19)") %>%
  mutate(pAB=(`0|1,0|1` + `1|0,1|0`)/100) %>%
  mutate(pA=(`0|1,0|1` + `0|1,0|0` + `1|0,0|0` + `1|0,1|0` +`1|0,0|1`+`0|1,1|0` +2*`1|1,0|0` + 2*`1|1,1|1`)/100) %>%
  mutate(pB= (`0|0,0|1` + `0|0,1|0` + `0|1,0|1` + `1|0,1|0` + 2*`0|0,1|1` + `1|0,0|1` + `0|1,1|0` + 2*`1|1,1|1`)/100) %>%
  mutate(d=pAB-pA*pB) %>%
  mutate(r_squared = d**2/(pA*(1-pA)*pB*(1-pB))) %>%
  mutate(d_prime = case_when(
    d < 0 ~ d / pmin(pA*pB , (1-pA)*(1-pB)),
    d > 0 ~ d / pmin(pA*(1-pB), (1-pA)*pB )
    )) %>% 
  mutate(distance=abs(pos_2-pos_1)) %>%
  mutate(recombination_rate="unknown") %>%
  filter(allele_count==2)

## low coverage hg19


long_low_genotype_config<-annotate_genotypes(long_low_genotype_config, type="phased")


annotated_genotypes<-long_low_genotype_config

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

 
 
 dprime_low_coverage_df<-low_coverage_hg19 %>%
  filter(allele_count == allele_count_to_analyze)  %>%
      mutate(distance_breaks=cut_interval(distance, length=window_size)) %>%
  group_by(distance_breaks, variation_type, recombination_rate) %>%
  summarise(mean_dprime=mean(d_prime, na.rm=T), mean_d=mean(d), mean_rsquare=mean(r_squared, na.rm=T), count=length(d)) %>% 
    ungroup()

    
    prow_1<-plot_grid( sim_plot_1+ theme(legend.position="none"),
           sim_plot_2 + theme(legend.position = "none") + labs(y=""),
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
si_fig_10 <- plot_grid( prow_1, legend,ncol = 1, rel_heights = c(1, .1))
si_fig_10


ggsave(filename="figures/si_figure_10_correlations_dprime_hrj", plot=si_fig_10, width=20, height=12)
