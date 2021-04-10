library(tidyverse)
library(cowplot)


annotate_genotypes<-function(long_genotype_config){

    long_genotype_config<-long_genotype_config %>%
      mutate(words=case_when(
        genotype == "0/0,0/0" ~ "Homozygous reference (0/0,0/0)", 
        genotype == "0/0,0/1" ~ "Single heterozygote (0/1,0/0 or 0/0,0/1)", 
        genotype == "0/1,0/0" ~ "Single heterozygote (0/1,0/0 or 0/0,0/1)",
        genotype == "0/1,0/1" ~ "Double Heterozygote (0/1,0/1)"
      ))
}

compute_hrj_unphased<-function(genotype_counts, window_size, allele_count_to_analyze){
  
  long_genotype_configuration <- gather(genotype_counts, `0/0,0/0` ,`0/0,0/1` ,`0/1,0/0`, `0/1,0/1`,`0/0,1/1` ,`1/1,0/0` ,`1/1,1/1`  , `1/1,0/1`, `0/1,1/1` , 
                                        value="count", 
                                        key="genotype")
  
  annotated_genotype_configuration <- annotate_genotypes(long_genotype_configuration)
  

  double_heterozygote <- "Double Heterozygote (0/1,0/1)"

  hrj <- annotated_genotype_configuration %>%
      filter(allele_count == allele_count_to_analyze, !is.na(words), words == double_heterozygote)  %>%
      mutate(distance_breaks=cut_interval(distance, length=window_size)) %>%
      group_by(distance_breaks, variation_type, words) %>%
      summarise(hrj=mean(count), n=n()) %>% 
      ungroup() 

  return(hrj)
}


## Example 
high_coverage_hg38_genotype_counts <- read_csv("masked_high_coverage_hg38_genotype_counts.csv") 
highcoverage_hg38_doubletons_hr <- compute_hrj_unphased(high_coverage_hg38_genotype_counts  , window_size = 1500, allele_count_to_analyze = 2)


