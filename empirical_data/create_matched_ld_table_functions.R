library(tidyverse)
library(furrr)
library(cowplot)
library(ggplot2)
library(MatchIt)
library(glue)


## Step 06 
createVCF<-function( variation_type, allele_count, chromosome,population_code){

	
    vcf_input_path<-glue("/u/home/j/jessegar/project-klohmuel/SLiM_ParallelRProject/data/04_chr{chromosome}_{population_code}_PolarizedAndACFix.vcf")
	output_vcf<-glue("../../data/filtered_vcf/{variation_type}_ac_{allele_count}_chromosome{chromosome}_population_{population_code}_10000K.vcf")
	
	command_to_filter<-glue("grep -E '#|AC={allele_count};AF=.*;.*ExonicFunc.refGene={variation_type}' {vcf_input_path} > {output_vcf}")
	system(command_to_filter)

	command_to_filter
}

## Step 07
compute_ld<-function( vcftools_path = "/u/local/apps/vcftools/0.1.9/bin/vcftools",variation_type, allele_count, chromosome,population_code){

	
	vcf_input_path<-glue("../../data/filtered_vcf/{variation_type}_ac_{allele_count}_chromosome{chromosome}_population_{population_code}_10000K.vcf")
	
	command_to_compute<-glue("{vcftools_path} --vcf {vcf_input_path} --hap-r2 --ld-window-bp 100000 --out ../../data/ld/{variation_type}_ac_{allele_count}_chromosome{chromosome}_population_{population_code}")
	system(command_to_compute)

	command_to_compute
}

## Step 08
compute_genetic_map<-function( python="/usr/bin/python", 
interpolate_genetic_distance_script = "interpolate_genetic_distance_LD.py",
variation_type, 
allele_count, 
chromosome,
population_code){

	
	vcftools_ld_input_path<-glue(" ../../data/ld/{variation_type}_ac_{allele_count}_chromosome{chromosome}_population_{population_code}.hap.ld")
	genetic_map_path<-glue("/u/home/j/jessegar/project-klohmuel/LD1000GenomePolarizedVCF/GeneticMap/chr{chromosome}_average_noncarrier.gmap")
	output<-glue("../../data/ld//{variation_type}_ac_{allele_count}_chromosome{chromosome}_population_{population_code}_map.txt")
	command_to_compute<-glue("{python} {interpolate_genetic_distance_script} --genetic_map {genetic_map_path} --coordinates {vcftools_ld_input_path} --outfile {output}")
	system(command_to_compute)

	command_to_compute
}


## Step 09
interpolate_with_genetic_map<-function( python="/usr/bin/python", 
convert_physical_to_genetic_script = "convert_physical_to_genetic.py",
variation_type, 
allele_count, 
chromosome,
population_code){

	
	vcftools_ld_input_path<-glue(" ../../data/ld/{variation_type}_ac_{allele_count}_chromosome{chromosome}_population_{population_code}.hap.ld")
	genetic_map_path<-glue("/u/home/j/jessegar/project-klohmuel/LD1000GenomePolarizedVCF/GeneticMap/chr{chromosome}_average_noncarrier.gmap")
	specific_genetic_map_path<-glue("../../data/ld/{variation_type}_ac_{allele_count}_chromosome{chromosome}_population_{population_code}_map.txt")
	output<-glue("../../data/ld/{variation_type}_ac_{allele_count}_chromosome{chromosome}_population_{population_code}_genetic_distance_ld.txt")
	command_to_compute<-glue("{python} {convert_physical_to_genetic_script} --genetic_map {specific_genetic_map_path} --VCFtoolsout {vcftools_ld_input_path} --outfile {output}")
	system(command_to_compute)

	command_to_compute
}



## Step 10
combine_ld_annotated_with_genetic_distance<-function(
variation_type, 
allele_count, 
chromosome,
population_code){

	
	vcftools_ld_genetic_map_input_path<-glue("../../data/ld/{variation_type}_ac_{allele_count}_chromosome{chromosome}_population_{population_code}_genetic_distance_ld.txt")
	genetic_map_ld<-read_tsv(vcftools_ld_genetic_map_input_path, col_names=F) %>% mutate(allele_count = allele_count, chromosome = chromosome, variation_type=variation_type, population_code=population_code)
	return(genetic_map_ld)
}


## Step 11
source("functions_for_appending_bvalues.R")




## Step 12
prepping_ld_table<-function(population_code){
ld_table<-glue("../../data2/ld_{population_code}_annotated_bvalues_genetic_distance.rds")
ACFilteredLDBValues<-read_rds(ld_table) %>% dplyr::rename(genetic_distance=X5,physical_distance=X4) 
if(ACFilteredLDBValues$variation_type %>% unique() %>% length() !=  2){
	  stop("Data is not formatted Correctly. There are more than just nonsynonymous and synonymous SNVs.")
}
ACFilteredLDBValues<-ACFilteredLDBValues %>% mutate(Group=ifelse(variation_type == "nonsynonymous_SNV", F, T))
ACFilteredLDBValues$mean_b_value<-(ACFilteredLDBValues$BValuePOS1 + ACFilteredLDBValues$BValuePOS2)/2
ACFilteredLDBValues<-na.omit(ACFilteredLDBValues)
return(ACFilteredLDBValues)
}

  
match_on_b_values_ac_chromosome_physical_distance<-function(population_code, number_of_b_value_breaks, number_of_physical_distance_breaks){

ld_table<-prepping_ld_table(population_code=population_code) %>% filter(physical_distance < 10000)
ld_table$physical_distance_breaks<-cut(ld_table$physical_distance, breaks=number_of_physical_distance_breaks,include.lowest = T)
ld_table$b_value_breaks<-cut(ld_table$mean_b_value, breaks=number_of_b_value_breaks,include.lowest = T)

#Assigning all observations with the same AC, Bvalue Bin, and PhysicalDistance Bin the same label
ld_table_temp<-ld_table %>% mutate(label=group_indices(. ,allele_count,b_value_breaks,Chromosome, physical_distance_breaks))

#Matching Data based on Group Label NEEDS TO BE A DATAFRAME TO WORK

matching<-matchit(Group~label, as.data.frame(ld_table_temp), ratio=1.0, exact="label", replace=F)

#Matching Indices
controlRows<-matching$match.matrix[,1]

#Remove Unmatched Synonymous SNPS
matches<-controlRows[!is.na(controlRows)]

#Names of matches is the synonymous and the actual value is the matching nonsynonymous pair
synonymousPairs<-ld_table_temp[names(matches),]
synonymousPairs<-merge(synonymousPairs, ld_table, sort=F)

nonsynonymousPairs<-ld_table_temp[matches,]
nonsynonymousPairs<-merge(nonsynonymousPairs, ld_table, sort=F)

allMatches<-rbind(synonymousPairs, nonsynonymousPairs)


allMatches %>% 
write_rds(glue("../../data2/{population_code}_b_value_breaks_{number_of_b_value_breaks}_physical_breaks_{number_of_physical_distance_breaks}.rds" ) ) 
return(as_tibble(allMatches))

}

match_on_ac_chromosome_physical_distance<-function(population_code, number_of_physical_distance_breaks){

ld_table<-prepping_ld_table(population_code=population_code) %>% filter(physical_distance < 10000)
ld_table$physical_distance_breaks<-cut(ld_table$physical_distance, breaks=number_of_physical_distance_breaks,include.lowest = T)

#Assigning all observations with the same AC, Bvalue Bin, and PhysicalDistance Bin the same label
ld_table_temp<-ld_table %>% mutate(label=group_indices(. ,allele_count,Chromosome, physical_distance_breaks))

#Matching Data based on Group Label NEEDS TO BE A DATAFRAME TO WORK

matching<-matchit(Group~label, as.data.frame(ld_table_temp), ratio=1.0, exact="label", replace=F)

#Matching Indices
controlRows<-matching$match.matrix[,1]

#Remove Unmatched Synonymous SNPS
matches<-controlRows[!is.na(controlRows)]

#Names of matches is the synonymous and the actual value is the matching nonsynonymous pair
synonymousPairs<-ld_table_temp[names(matches),]
synonymousPairs<-merge(synonymousPairs, ld_table, sort=F)

nonsynonymousPairs<-ld_table_temp[matches,]
nonsynonymousPairs<-merge(nonsynonymousPairs, ld_table, sort=F)

allMatches<-rbind(synonymousPairs, nonsynonymousPairs)


allMatches %>% 
write_rds(glue("../../data2/{population_code}_physical_breaks_{number_of_physical_distance_breaks}.rds" ) ) 
return(as_tibble(allMatches)) 

}
