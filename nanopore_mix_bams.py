#!/usr/bin/env python


import argparse
import pysam
import csv
import sys
import random
import os.path


def parse_commandline():
	parser=argparse.ArgumentParser()
	parser.add_argument('--bam1', help='first bam file',required=True)
	parser.add_argument('--bam2', help='second bam file',required=True)
	parser.add_argument('--out_prefix',help='output prefix',required=True)
	parser.add_argument('--frac', help='mixing fraction (p <= bam 1/total)',default=0.5)
	parser.add_argument('--nreads', help='number of reads to generate',default=1000000)
	parser.add_argument('--threads',help='number of threads',required=False, default=1)
	parser.add_argument('--seed',help='random seed',required=False, default=1)
	args=parser.parse_args()
	print(args, file=sys.stderr)
	return args

def main(args):

	random.seed(int(args.seed))
	
	bam_file_1 = pysam.AlignmentFile(args.bam1, 'rb', threads=int(args.threads))
	bam_file_2 = pysam.AlignmentFile(args.bam2, 'rb', threads=int(args.threads))
	out_file = pysam.AlignmentFile(args.out_prefix + "_" + str(args.frac) + "_" + str(args.nreads) + "_" + str(args.seed) + ".bam", 'wb', threads=int(args.threads), template=bam_file_1)
		
	#random.seed(0)
	
	with open(args.out_prefix + "_" + str(args.frac) + "_" + str(args.nreads) + "_" + str(args.seed) + ".mix.txt",'w') as w:
		csv_w = csv.writer(w, delimiter='\t', quoting=csv.QUOTE_NONE)
		i = 0
		for (read_1,read_2) in zip(bam_file_1.fetch(until_eof=True), bam_file_2.fetch(until_eof=True)):
			if i > int(args.nreads):
				break
				
			rand_num = random.random()
			
			if read_1.is_unmapped or read_2.is_unmapped:
				continue
				
			if read_1.is_supplementary or read_1.is_secondary or read_2.is_supplementary or read_2.is_secondary:
				continue
	
			if rand_num <= float(args.frac):
				csv_w.writerow([read_1.query_name, os.path.basename(args.bam1)])
				out_file.write(read_1)
				i += 1
			else:
				csv_w.writerow([read_2.query_name, os.path.basename(args.bam2)])
				out_file.write(read_2)
				i += 1
	
				
if __name__ == '__main__':
	args=parse_commandline()
	main(args)
