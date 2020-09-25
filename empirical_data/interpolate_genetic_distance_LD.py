import os
import sys
import argparse
import csv
from bisect import bisect

def parse_args():
	"""
	Parse command-line arguments
	"""
	parser = argparse.ArgumentParser(description="This script takes in the genetic map and obtains the genetic distance for a position if that position is not in the genetic map by interpolation.")
	parser.add_argument(
			"--genetic_map", required=True,
			help="REQUIRED. Genetic map. This is the DECODE map. There are 4 columns. The columns are (1) chr, (2) RSID, (3) position, (3) cM defined as the distance between the current SNP and the previous SNP.")
	parser.add_argument(
			"--coordinates", required=True,
			help="REQUIRED. Coordinates. For this script, this is the output from VCFtool where the header is (1) chr, (2) pos, (3) pos, (5)r2 ")
	parser.add_argument(
			"--outfile", required=True,
			help="REQUIRED. Name of the output file")

	args = parser.parse_args()
	return args

def interpolate(known_cM_positions, unknown_cM_positions, known_cM, outfile):
	"""
	This function takes in:
	1. a list of positions where cM is known (known_cM_positions)
	2. a list of positions where cM is not known (unknown_cM_positions)
	3. a list of cM (known_cM)(same index as known_cM_positions)
	"""

	# for each position where cM is now known, find the index where it will be inserted
	for unknown_cM_position in unknown_cM_positions:
		# need to check so that will only consider the coordinates that do not fall into the NA 
		if unknown_cM_position > known_cM_positions[0] and unknown_cM_position < known_cM_positions[-1]:
			index_at_insert = bisect(known_cM_positions, unknown_cM_position)
			cM_before = known_cM[index_at_insert-1]
			cM_after = known_cM[index_at_insert]

			# compute (p_C - p_A)/(p_B - p_A)
			diff_physical_dist = float(unknown_cM_position - known_cM_positions[index_at_insert-1])/(known_cM_positions[index_at_insert] - known_cM_positions[index_at_insert-1])

			cM_interpolate = known_cM[index_at_insert]*diff_physical_dist
			known_cM[index_at_insert] = cM_after - cM_interpolate # need to update because cM is defined as the distance between the current SNP and previous SNP. Need to update so that the meaning of cM is consisten. 

			known_cM_positions.insert(index_at_insert, unknown_cM_position)

			known_cM.insert(index_at_insert, cM_interpolate)

			# print unknown_cM_position, cM_interpolate

	for i in range(len(known_cM_positions)):
		toprint = [str(known_cM_positions[i]), str(known_cM[i])]
		print >>outfile, "\t".join(x for x in toprint)
def main():
	args = parse_args()

	known_cM_positions = []
	known_cM = []
	unknown_cM_positions_set = set() # use set here because I want to curate the VCFtool output to get a list of SNPs that are unique

	with open(args.genetic_map, "r") as f:
		for line in f:
			line = line.rstrip("\n")
			line = line.split("\t")
			if line[3] != 'NA':
				known_cM_positions.append(int(line[2]))
				known_cM.append(float(line[3]))

	with open(args.coordinates, "r") as f:
		for line in f:
			if not line.startswith('CHR'):
				line = line.rstrip("\n")
				line = line.split("\t")
				unknown_cM_positions_set.add(int(line[1]))
				unknown_cM_positions_set.add(int(line[2]))

	unknown_cM_positions = list(unknown_cM_positions_set) # convert back to a list

	outfile = open(args.outfile, "w")

	interpolate(known_cM_positions, unknown_cM_positions, known_cM, outfile)

main()
