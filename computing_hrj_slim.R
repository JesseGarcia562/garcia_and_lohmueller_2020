library(data.table)
library(tidyverse)

library(reticulate)
library(Rcpp)
library(glue)

read_vcf<-function( vcf_input_path,allele_count,variation_type){

       
        vcf_input_path<-vcf_input_path
        raw_lines=read_lines(vcf_input_path)
        last_line=length(raw_lines)
        header=str_detect(raw_lines, pattern="CHR") %>% which()
	variants_with_annotation<-str_detect(raw_lines, glue("MT={variation_type};AC={allele_count};DP=1000")) %>% which()
	
	# If no variants with annotation
	if ( is_empty(variants_with_annotation) ) {
	lines_to_read<-c(header[1], header[1] + 1)
	lines_to_read=raw_lines[lines_to_read]
	} else {
	
	lines_to_read<-c(header[1], variants_with_annotation)
	lines_to_read=raw_lines[lines_to_read]
        	
	}
	
	read_tsv(lines_to_read, col_names=TRUE)
}




find_pairs_of_snps<-function(vcf, distance_limit){
possible_combinations<-combn(vcf$POS, m=2, simplify=FALSE)
distance<-possible_combinations %>% map(~ abs(.x[1] - .x[2]))
distance<-unlist(distance)
combinations<-possible_combinations[distance < distance_limit]
return(combinations)
}


get_genotypes<-function(combinations, vcf){
genotypes<-combinations %>% map(~ vcf %>% filter(POS ==.x[1] | POS==.x[2]))  

return(genotypes)
}


get_genotype_distribution<-function(genotype, allele_count){



genotype_distribution<-genotype %>% summarise_at(vars(starts_with("i")), list(~ paste(., collapse=","))) %>% select(-ID, -INFO) %>%

	unlist(., use.names=FALSE)  %>%
	   table()




}


turn_genotype_distribution_to_df<-function(genotype_distribution){

	column_names<-names(genotype_distribution)
	values_for_df<-unname(genotype_distribution)

	names(values_for_df)<-column_names

	## Named vectors must be turned to  lists, then tibbles
	return(values_for_df %>% as.list() %>% as_tibble()  )
}



combine_combinations_genotype_distribution_df<-function(combinations, genotype_distribution_df, variation){

	
	combination<-c(combinations, variation)
	df<-tibble(
		   pos_1=as.integer(combinations[1]),
		   pos_2=as.integer(combinations[2]),
		   variation_type=variation
		   )
	return( bind_cols(df, genotype_distribution_df )) 	
} 

get_configuration_table<-function(vcf_input_path,recombinationRate,variation_type, allele_count, chromosome, population, distance_limit){

 # browser()
  vcf<-read_vcf(vcf_input_path=vcf_input_path,allele_count=allele_count,variation_type=variation_type ) 

combinations<-find_pairs_of_snps(vcf=vcf, distance_limit=distance_limit)
genotypes<-get_genotypes(combinations=combinations, vcf=vcf)

genotype_distribution<-genotypes %>% map(~get_genotype_distribution(genotype=.x, allele_count=allele_count))

genotype_distribution_df<- genotype_distribution %>%
	map(~ turn_genotype_distribution_to_df(.x))

configuration_table<-map2_df(.x=combinations, 
			     .y=genotype_distribution_df , 
			     ~combine_combinations_genotype_distribution_df(combinations=.x, genotype_distribution_df=.y, variation=variation_type))

return(configuration_table %>% mutate(chromosome=chromosome, allele_count=allele_count,recombination_rate=recombinationRate))
}

safely_get_configuration_table<-safely(get_configuration_table)



## Usage for slim

ns_config_doubletons<-safely_get_configuration_table(variation_type = 1, allele_count=2, chromosome=1,recombinationRate=.y, population="sim", distance_limit=100000,vcf_input_path=.x)

