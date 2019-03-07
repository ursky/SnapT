#!/usr/bin/env python2
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

def load_contigs_from_gff(filename):
	contigs={}
	for line in open(filename):
		if line[0]=="#":
			continue
		cut=line.strip().split("\t")
		if len(cut)!=9:
			continue
		contig=cut[0]
		if "_len" not in contig:
			print "could not determine length of contig from name: " + contig
			continue
		cut_c=contig.split("_")
		length=int(cut_c[3])
		if contig not in contigs:
			contigs[contig]=length
	return contigs
		

def load_transcripts(filename, lengths):
	distances={}
	for line in open(filename):
		if line[0]=="#":
			continue
		cut=line.strip().split("\t")
		if len(cut)!=9:
			continue
		if cut[2]=="antisense transcript":
			continue
		contig=cut[0]
		st = int(cut[3])
		fi = int(cut[4])
	
		distance_to_end = min(st, lengths[contig]-fi)
		#print lengths[contig], st, fi, distance_to_end
		if contig not in distances:
			distances[contig]=[]
		distances[contig].append(distance_to_end)
	return distances


def filter_transcripts(filename, lengths, min_dist, min_length, outfile):
	print "Filtering out transcripts that are less than "+str(min_dist)+"bp away from a contig edge..."
	f = open(outfile, "w")
	ct_end=0
	ct_good=0
	ct_short=0
	ct_exons=0
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
		
		length = lengths[contig]
		distance_to_end = min(st, length-fi)
		
		if cut[2]=="EXON" or cut[2]=="Exon" or cut[2]=="exon":
			ct_exons+=1
		elif length<min_length:
			ct_short+=1
		elif distance_to_end<min_dist:
			ct_end+=1
		else:
			f.write(line)
			ct_good+=1
	print "Removed "+str(ct_short)+" transcripts that were on contigs that were too short (<"+str(min_length)+"), and "+str(ct_end)+" transcripts that were too close to contig edges. "+str(ct_good)+" good ncRNA transcripts remaining!"
	f.close()
			

def compute_statistics(contig_lengths, distances_to_edge, cutoff):
	contig_length_bins=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
	n_bins = len(contig_length_bins)


	for contig in distances_to_edge:
		for distance in distances_to_edge[contig]:
			if distance<cutoff:
				continue
			relative_position = n_bins*2*distance/contig_lengths[contig]
			contig_length_bins[relative_position]+=1
	for i in range(n_bins):
		contig_length_bins[i]=str(contig_length_bins[i])
	return contig_length_bins




def optimize_cut_off(contig_lengths, distances_to_edge):
	cut_off=0
	while True:
		if cut_off>300:
			break
		print "Computing overrepresentation at cut off: "+ str(cut_off)
		representation = compute_statistics(contig_lengths, distances_to_edge, cut_off)
		print " ".join(representation)
		if representation[2]==0 or representation[3]==0:
			cut_off+=10
			continue
		
		ratio = float(representation[0])/float(representation[1])
		if ratio<1.25:
			break
		else:
			cut_off+=10
	print "Looks like using a minimum distance of "+str(cut_off+10)+" is safe!"
	return cut_off+10
		

#MAIN
genome_file=sys.argv[1]
annotation_file=sys.argv[2]
out_file=sys.argv[3]
min_length = int(sys.argv[4])

contig_lengths = load_contigs(genome_file)
distances_to_edge = load_transcripts(annotation_file, contig_lengths)
cut_off = optimize_cut_off(contig_lengths, distances_to_edge)
filter_transcripts(annotation_file, contig_lengths, cut_off, min_length, out_file)





