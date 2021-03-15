# Garcia and Lohmueller 2020
Scripts for Garcia and Lohmueller 2020


# Computing Hrj with SLIM VCFs

To compute Hrj for VCFs we used the R Based script titled: "computing_hrj_slim.R"

Their are seven arguments needed for the function get_configuration_table(): 

variation_type (integer) = This is the mutation type that was used in slim that is to be analyzed. In our simulations we used three mutation types (1=nonsynonymous, 2=synonymous and 3=intergenic). To study nonsynonymous mutations, this argument should be 1. 

allele_count (integer) = This is the allele count in the sample to be analyzed. In order to study doubletons (variants that show up twice in a sample) the argument should be 2. 

chromosome (integer) = This value indicates the chromosome number of the slim simulation. This value should be 1 for slim 3.

recombinationRate (float) = To annotate the data this argument records the recombination rate of the simulation. In simulations with a recombination rate of 1e-08 crossovers per base pair per generation this value will be equal to 1e-08. 

population (character) = This argument annotates the data with the population that the genetic information comes from. A usual parameter is "sim". 

distance_limit (integer) = This argument limits the base pairs allowed between two loci in the computation. If Loci 1 is at base pair 50 and Loci 2 is at base pair 10,500, the distance between them is 10,500 - 50 = 10,450 base pairs. 

vcf_input_path (character) = This argument indicates the path towards the VCF to be analyzed. 


# Computing Hrj with other data
