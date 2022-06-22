#!/usr/bin/env python


import argparse
import pysam
import csv
import sys
from modbampy import ModBam
import pandas

def parse_commandline():
	parser=argparse.ArgumentParser()
	parser.add_argument('--bam', help='bam file',required=True)
	parser.add_argument('--fasta', help='fasta file',required=True)
	parser.add_argument('--out',help='output table',required=True)
	parser.add_argument('--nthreads',help='number of threads',required=False, default=1)
	args=parser.parse_args()
	print(args, file=sys.stderr)
	return args

def main(args):
	
	reference_fasta = pysam.FastaFile(args.fasta)
	
	#modbam needs to get all contig lengths first
	tmp_bam = pysam.AlignmentFile(args.bam, "rb")
	
	rlist = []
	for i in range(tmp_bam.nreferences):
		rname = tmp_bam.get_reference_name(i)
		rlength = tmp_bam.get_reference_length(rname)
		rlist.append([rname,rlength])
	with open(args.out,'w') as w:
		csv_w = csv.writer(w, delimiter='\t', quoting=csv.QUOTE_NONE)
		
		for rname,rlength in rlist:
			#thanks to ONT for their nice modbampy library! https://github.com/epi2me-labs/modbam2bed/issues/2
			with ModBam(args.bam, rname, 0, int(rlength)) as bam:
			
				
				for read in bam.reads():
					read_data = []
				
					length_pos_mod = 0 #number of candidate sites
					for pos_mod in read.mod_sites():
						qname, rpos, qpos, strand, mod_strand, cbase, mbase, score = pos_mod
						#if (strand == "+" and reference_fasta.fetch(rname,rpos,rpos+2) == 'CG'):                
						#	if score >= 0.8 * 255:
						#		read_data.append([rname, rpos, rpos+1, 1, strand])
						#	else:
						#		read_data.append([rname, rpos, rpos+1, 0, strand])
						#	length_pos_mod += 1
						#if (strand == "-" and reference_fasta.fetch(rname,rpos-1,rpos+1) == 'GC'):
						#	if score >= 0.8 * 255:
						#		read_data.append([rname, rpos-1, rpos, 1, strand])
						#	else:
						#		read_data.append([rname, rpos-1, rpos, 0, strand])
						if (score <= 0.2 * 255):
							read_data.append([rname, rpos, rpos+1, 0, score, strand, reference_fasta.fetch(rname,rpos,rpos+2), reference_fasta.fetch(rname,rpos-1,rpos+1)])
						if (score >= 0.8 * 255):
							read_data.append([rname, rpos, rpos+1, 1, score, strand, reference_fasta.fetch(rname,rpos,rpos+2), reference_fasta.fetch(rname,rpos-1,rpos+1)])
							length_pos_mod += 1
						
		
					if len(read_data) == 0:
						continue
				
					for row in read_data:
						csv_w.writerow(row + [qname])
	


				
if __name__ == '__main__':
	args=parse_commandline()
	main(args)
