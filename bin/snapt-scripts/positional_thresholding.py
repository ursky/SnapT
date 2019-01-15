#!/usr/bin/env python
import sys
# usage: positional_thresholding.py assembly.fa nc_rna.gff out_file.gff
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


def filter_transcripts(filename, lengths, min_dist, outfile):
	print "Filtering out transcripts that are less than "+str(min_dist)+"bp away from a contig edge..."
	f = open(outfile, "w")
	ct=0
	for line in open(filename):
		if line[0]=="#":
			f.write(line)
			continue
		cut=line.split("\t")
		if len(cut)!=9:
			f.write(line)
			continue
		contig=cut[0]
		st = int(cut[3])
		fi = int(cut[4])
		
		distance_to_end = min(st, lengths[contig]-fi)
		if distance_to_end > min_dist:
			f.write(line)
		else:
			ct+=1
	print "Removed "+str(ct)+" transcripts that were too close to contig edges."
	f.close()
			

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
		print "Computing overrepresentation at cut off: "+ str(cut_off)
		representation = compute_statistics(contig_lengths, positions, cut_off)
		print representation[0], representation[1], representation[2], representation[3], representation[4]
		if representation[2]==0 or representation[3]==0:
			cut_off+=10
			continue
		ratio2 = representation[1]/representation[2]
		ratio3 = representation[2]/representation[3]

		if ratio2<1.5 and ratio2>0.75 and ratio3<1.5 and ratio3>0.75:
			break
		cut_off+=10
		if cut_off>500:
			break
	print "Looks like using a minimum distance of "+str(cut_off+10)+" is safe!"
	return cut_off+10
		

contig_lengths = load_contigs(sys.argv[1])
positions = load_transcripts(sys.argv[2], contig_lengths)
cut_off = optimize_cut_off(contig_lengths, positions)
filter_transcripts(sys.argv[2], contig_lengths, cut_off, sys.argv[3])

