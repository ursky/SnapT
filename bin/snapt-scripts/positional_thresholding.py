#!/usr/bin/env python
import sys
# usage: positional_thresholding.py assembly.fa nc_rna.gff 
# This script computes the optimal cut-off for the minimum distance of a ncRNA to the contig edge

def load_contigs(filename):
	contigs={}
	for line in open(filename):
		if line[0]==">":
			contig=line[1:-1]
			contigs[contig]=0
		else:
			contigs[contig]+=len(line.strip())
	return contigs


def load_transcripts(filename, lengths):
	distances=[]
	for line in open(filename):
		if line[0]=="#":
			continue
		cut=line.strip().split("\t")
		if len(cut)!=9:
			continue
		contig=cut[0]
		st = int(cut[3])
		fi = int(cut[4])
	
		distance_to_end = min(st, lengths[contig]-fi)
		#print lengths[contig], st, fi, distance_to_end
		distances.append(distance_to_end)
	return distances


def filter_transcripts(filename, lengths, min_dist):
	for line in open(filename):
		line=line.strip()
		if line[0]=="#":
			print line
			continue
		cut=line.split("\t")
		if len(cut)!=9:
			print line
			continue
		contig=cut[0]
		st = int(cut[3])
		fi = int(cut[4])
		
		distance_to_end = min(st, lengths[contig]-fi)
		if distance_to_end > min_dist:
			print line
			

def compute_statistics(contig_lengths, positions, bin_size):
	contig_length_bins=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
	for contig in contig_lengths:
		length=contig_lengths[contig]
		n_bins = int(0.5*length/bin_size)
		for i in range(n_bins):
			if i>=len(contig_length_bins):
				break
			contig_length_bins[i]+=1

	trans_position_bins=[]
	for i in range(len(contig_length_bins)):
		trans_position_bins.append(0)
	for transcript in positions:
		bin_pos = int(transcript/bin_size)
		if bin_pos<len(trans_position_bins):
			trans_position_bins[bin_pos]+=1

	bin_representation=[]
	for i, n_transcripts in enumerate(trans_position_bins):
		n_contigs = contig_length_bins[i]
		bin_representation.append(float(n_transcripts)/n_contigs)
	
	return bin_representation


def optimize_cut_off(contig_lengths, positions):
	cut_off=10
	while True:
		# computing overrepresentation at cut off cut_off
		representation = compute_statistics(contig_lengths, positions, cut_off)
		ratio1 = representation[0]/representation[1]
		ratio2 = representation[1]/representation[2]
		ratio3 = representation[2]/representation[3]
		ratio4 = representation[3]/representation[4]

		if ratio2<1.5 and ratio2>0.75 and ratio3<1.5 and ratio3>0.75:
			break
		cut_off+=10
		if cut_off>1000:
			break
	return cut_off+10
		

contig_lengths = load_contigs(sys.argv[1])
positions = load_transcripts(sys.argv[2], contig_lengths)
cut_off = optimize_cut_off(contig_lengths, positions)
filter_transcripts(sys.argv[2], contig_lengths, cut_off)

