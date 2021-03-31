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

## Input data

To compute Hrj with data types that are not slim VCF's, first one must create a data frame with 9 columns.

pos_1 (double) = The position (bp) of the left most variant.
pos_2 (double) = The position (bp) of the right most variant.
variation_type (character) = Defines the annotation of both variants. Either "nonsynonymous_SNV" or "synonymous_SNV".
chromosome (double) = The chromosome number that both variants exist on. 
allele_count (double) = The number of times the variants exists in the sample. For doubletons this would be 2, and for singletons this would be 1.
data_type (character) = Metadata defining where the data and genotype call are coming from. 
`0/0,0/0` (double) = The number of individuals in the sample that are homozygous for the ancestral allele at pos_1 and pos_2
`0/0,0/1` (double) = The number of individuals in the sample that are homozygous for the ancestral allele at pos_1 and heterozygous for the ancestral allele at pos_2
`0/1,0/0` (double) = The number of individuals in the sample that are homozygous for the ancestral allele at pos_1 and heterozygous for the ancestral allele at pos_2
`0/1,0/1` (double) = The number of individuals in the sample that are heterozygous for the ancestral allele at pos_1 and pos_2
`0/0,1/1` (double) = The number of individuals in the sample that are homozygous for the ancestral allele at pos_1 and homozygous for the derived allele at pos_2
`1/1,0/0` (double) = The number of individuals in the sample that are homozygous for the derived allele at pos_1 and homozygous for the ancestral allele at pos_2
`1/1,1/1` (double) = The number of individuals in the sample that are homozygous for the derived allele at pos_1 and pos_2
`1/1,0/1` (double) = The number of individuals in the sample that are homozygous for the derived allele at pos_1 and heterozygous for the ancestral allele at pos_2
`0/1,1/1` (double) = The number of individuals in the sample that are heterozygous for the derived allele at pos_1 and homozygous for the ancestral allele at pos_2



# Simulation scripts

In the folder simulation_scripts, scripts to simulate Model 2 and Model 3 exist. Each row of the csv contains the script needed for a simulation replicate.  
