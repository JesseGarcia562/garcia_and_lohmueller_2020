library(glue)  
library(tidyverse)
library(furrr)
library(glue)
library(pryr)
library(MatchIt)
plan(multiprocess)

chromosome<-1:22 
allele_count<-1:5

source("functions_for_appending_bvalues.R")
source("create_matched_ld_table_functions.R")
df<-crossing(chromosome, allele_count)
SGETaskID<-parse_integer(Sys.getenv("SGE_TASK_ID"))
population<-c("MXL", "YRI", "JPT", "CHB", "CEU")
population<-population[SGETaskID]



## Step 06
future_map2(.x = df$chromosome, .y=df$allele_count, ~ createVCF(population_code = population, chromosome = .x, variation_type = "nonsynonymous_SNV", allele_count=.y), .progress=T ) 
future_map2(.x = df$chromosome, .y=df$allele_count, ~ createVCF(population = population, chromosome = .x, variation_type = "synonymous_SNV", allele_count=.y), .progress=T ) 





## Step 07
future_map2(.x = df$chromosome, .y=df$allele_count, ~ compute_ld(population_code = population, chromosome = .x, variation_type = "nonsynonymous_SNV", allele_count=.y), .progress=T ) 

future_map2(.x = df$chromosome, .y=df$allele_count, ~ compute_ld(population_code = population, chromosome = .x, variation_type = "synonymous_SNV", allele_count=.y), .progress=T ) 

## Step 08

future_map2(.x = df$chromosome, .y=df$allele_count, ~ compute_genetic_map(population_code = population, chromosome = .x, variation_type = "nonsynonymous_SNV", allele_count=.y), .progress=T ) 


future_map2(.x = df$chromosome, .y=df$allele_count, ~ compute_genetic_map(population_code = population, chromosome = .x, variation_type = "synonymous_SNV", allele_count=.y), .progress=T ) 

## Step 09 

future_map2(.x = df$chromosome, .y=df$allele_count, ~ interpolate_with_genetic_map(population_code = population, chromosome = .x, variation_type = "nonsynonymous_SNV", allele_count=.y), .progress=T ) 


future_map2(.x = df$chromosome, .y=df$allele_count, ~ interpolate_with_genetic_map(population_code = population, chromosome = .x, variation_type = "synonymous_SNV", allele_count=.y), .progress=T ) 


## Step 10
nonsynymous_ld<-future_map2_dfr(.x = df$chromosome, .y=df$allele_count, ~ combine_ld_annotated_with_genetic_distance(population_code = population, chromosome = .x, variation_type = "nonsynonymous_SNV", allele_count=.y), .progress=T ) 
synonsymous_ld<-future_map2_dfr(.x = df$chromosome, .y=df$allele_count, ~ combine_ld_annotated_with_genetic_distance(population_code = population, chromosome = .x, variation_type = "synonymous_SNV", allele_count=.y), .progress=T ) 
ld_data<-bind_rows(synonsymous_ld, nonsynymous_ld) %>% write_rds(glue("../../data2/ld_{population}_annotated_low_coverage.rds"))

## Step 11
append_bvalues(hg19_bvalue_bed="../../data2/b_values_hg19.bed", empirical_ld_table = glue("/u/home/j/jessegar/project-klohmuel/SLiM_ParallelRProject/data2/ld_{population}_annotated_low_coverage.rds")) %>% 
write_rds(glue("../../data2/ld_{population}_annotated_bvalues_genetic_distance.rds"))


## Step 12

physical_bvalue=match_on_b_values_ac_chromosome_physical_distance(population_code=population, number_of_b_value_breaks=5, number_of_physical_distance_breaks=200)


physical=match_on_ac_chromosome_physical_distance(population_code=population,  number_of_physical_distance_breaks=200)
