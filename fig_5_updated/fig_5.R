library(tidyverse)
library(cowplot)


high_coverage_hg38_genotype_configuration <- read_csv("high_coverage_hg38_genotype_configuration.csv")

long_low_genotype_config <- read_csv("long_low_genotype_config.csv")

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
          summarise(mean_count=mean(count)) %>% 
          ungroup() %>%
          ggplot(aes(x=distance_breaks, y=mean_count, colour=variation_type, group=variation_type)) +
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

## High coverage hg38

highcoverage_hg38_doubletons_hr<-plot_rh_unphased(high_coverage_hg38_genotype_configuration %>% mutate(recombination_rate="unknown")%>% 
  mutate(variation_type=case_when(
  variation_type == "nonsynonymous_SNV" ~ "Nonsynonymous", 
  variation_type == "synonymous_SNV" ~ "Synonymous" 
)) , window_size = 1500, plot_type = "line_plot" , allele_count_to_analyze = 2)[[3]]




## low coverage hg19
low_coverage_doubletons_hr<-plot_rh_unphased(long_low_genotype_config %>% mutate(recombination_rate="unknown")%>% 
  mutate(variation_type=case_when(
  variation_type == "nonsynonymous_SNV" ~ "Nonsynonymous", 
  variation_type == "synonymous_SNV" ~ "Synonymous" 
)) , window_size = 1500, plot_type = "line_plot" , allele_count_to_analyze = 2)[[3]]




# extract a legend that is laid out horizontally
legend_b <- get_legend(
  low_coverage_doubletons_hr + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)




### just empirical

## High coverage hg38

highcoverage_hg38_singletons_hr<-plot_rh_unphased(high_coverage_hg38_genotype_configuration %>% mutate(recombination_rate="unknown")%>% 
  mutate(variation_type=case_when(
  variation_type == "nonsynonymous_SNV" ~ "Nonsynonymous", 
  variation_type == "synonymous_SNV" ~ "Synonymous" 
)) , window_size = 1500, plot_type = "line_plot" , allele_count_to_analyze = 1)[[3]]




## low coverage hg19
low_coverage_singletons_hr<-plot_rh_unphased(long_low_genotype_config %>% mutate(recombination_rate="unknown")%>% 
  mutate(variation_type=case_when(
  variation_type == "nonsynonymous_SNV" ~ "Nonsynonymous", 
  variation_type == "synonymous_SNV" ~ "Synonymous" 
)) , window_size = 1500, plot_type = "line_plot" , allele_count_to_analyze = 1)[[3]]



hr_plots<-plot_grid(
low_coverage_doubletons_hr + theme(legend.position = "none") + xlab(NULL)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + labs(title="Low coverage 1KGP"),
highcoverage_hg38_doubletons_hr + theme(legend.position = "none")+ scale_x_discrete(labels=c("1.5", "", "" ,"6" ,"","","10"))+ theme(axis.text.x = element_text(size=14, angle =45, hjust=1))  +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + ylab(" ")  + labs(title="High coverage NYGC 1KGP"), 
low_coverage_singletons_hr+ scale_x_discrete(labels=c("1.5", "", "" ,"6" ,"","","10"))  + theme(legend.position = "none"),
highcoverage_hg38_singletons_hr+ scale_x_discrete(labels=c("1.5", "", "" ,"6" ,"","","10"))  + theme(legend.position = "none") + ylab(" "),

labels = c('A', 'C', 'B', 'D'), label_size = 12, nrow=2, ncol=2)



plot_grid(hr_plots,legend_b, ncol=1, rel_heights = c(1, .1))
