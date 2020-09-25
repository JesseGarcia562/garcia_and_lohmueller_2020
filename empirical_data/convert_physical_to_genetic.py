import os
import sys
import argparse
import csv
from bisect import bisect

def parse_args():
	"""
	Parse command-line arguments
	"""
	parser = argparse.ArgumentParser(description="This script converts LD decay output from VCFtools from physical distance to genetic distance")
	parser.add_argument(
			"--genetic_map", required=True,
			help="REQUIRED. Genetic map. This is the output from the Python script interpolate_genetic_distance_LD.py. The first column is the position. The second column is the cM, defined as the genetic distance between the current SNP and the previous SNP.")
	parser.add_argument(
			"--VCFtoolsout", required=True,
			help="REQUIRED. Output from VCF tools.")
	parser.add_argument(
			"--outfile", required=True,
			help="REQUIRED. Name of the output file")	

	args = parser.parse_args()
	return args


def main():
	args = parse_args()

	outfile = open(args.outfile, "w")

	genetic_map_positions = []
	genetic_map_cM = []
	genetic_map_positions_set = set()

	with open(args.genetic_map, "r") as f:
		for line in f:
			line = line.rstrip("\n")
			line = line.split("\t")
			genetic_map_positions_set.add(int(line[0]))
			genetic_map_positions.append(int(line[0]))
			genetic_map_cM.append(float(line[1]))

	with open(args.VCFtoolsout, "r") as f:
		for line in f:
			if not line.startswith('CHR'):
				line = line.rstrip("\n")
				line = line.split("\t")

				SNP1 = int(line[1])
				SNP2 = int(line[2])
				if SNP1 in genetic_map_positions_set and SNP2 in genetic_map_positions_set:
					SNP1_index = genetic_map_positions.index(SNP1)
					SNP2_index = genetic_map_positions.index(SNP2)
					physical_distance = genetic_map_positions[SNP1_index:SNP2_index+1]
					cM = genetic_map_cM[SNP1_index:SNP2_index+1]
					total_physical_distance = 0
					total_cM = 0
					for i in range(1, len(physical_distance)):
						total_physical_distance += physical_distance[i] - physical_distance[i-1]
						total_cM += cM[i]
					toprint = [str(SNP1), str(SNP2), str(total_cM), str(total_physical_distance), str(line[4]), str(line[5]), str(line[6])]
					print >>outfile, "\t".join(x for x in toprint)

main()
