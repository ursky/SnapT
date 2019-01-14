#!/usr/bin/env python
import sys

def load_gff_coordinates(filename):
	coordinates={}
	# load coordinates from gff file:
	for line in open(sys.argv[1]):
		if line[0]=="#":
			continue
		cut=line.strip().split("\t")
		
		contig = cut[0]
		st = int(cut[3])
		fi = int(cut[4])
		strand = cut[6]
	
		if contig not in coordinates:
			coordinates[contig]={}
			coordinates[contig]["+"]=[]
			coordinates[contig]["-"]=[]

		coordinates[contig][strand].append((st,fi))
	return coordinates


def intersect(st, fi, strand_coordinates, dist):
	diagnosis="intergenic"
	for i,coordinate in enumerate(strand_coordinates):
		# middle of gene
		if st>=coordinate[0] and fi<=coordinate[1]:
			diagnosis="gene"
			break
		# start of gene
		elif fi<=coordinate[1] and fi>=coordinate[0]:
			diagnosis="gene"
			break
		# end of gene
		elif st<=coordinate[1] and st>=coordinate[0]:
			diagnosis="gene"
			break
		# complete span of gene
		elif st<coordinate[0] and fi>coordinate[1]:
			diagnosis="gene"
			break
		# intergenic on start of contig
		elif i==0 and fi+dist<coordinate[0]:
			diagnosis="intergenic"
			break
		# intergenic between genes
		elif fi+dist<coordinate[0] and st-dist>strand_coordinates[i-1][1]:
			diagnosis="intergenic"
			break
		# intergenic on end of contig
		elif i+1==len(strand_coordinates) and st-dist>strand_coordinates[i-1][1]:
			diagnosis="intergenic"
			break
		else:
			diagnosis="ambiguous"
	return diagnosis



def load_transcripts_gff(filename, orfs, dist):
	other_strand={"+":"-", "-":"+"}
	intergenic=""
	antisense=""

	for line in open(filename):
		if line[0]=="#":
			continue
		cut=line.strip().split("\t")
		if cut[2]!="transcript":
			continue

		contig = cut[0]
		st = int(cut[3])
		fi = int(cut[4])
		strand = cut[6]
	
		if contig not in orfs:
			continue
		if orfs[contig][strand]==[]:
			continue

		# check if this transcript is intergenic on the strand that its expressed on
		strand_coordinates = orfs[contig][strand]
		diagnosis = intersect(st, fi, strand_coordinates, dist)

		# if the transcript is intergenic, check if there is a gene on the other strand
		if diagnosis=="intergenic":
			strand_coordinates = orfs[contig][other_strand[strand]]
			opposite_strand_diagnosis = intersect(st, fi, strand_coordinates, 0)
			if opposite_strand_diagnosis=="gene":
				diagnosis = "antisense"

		if diagnosis=="intergenic":
			cut[2]="intergenic transcript"
			intergenic+="\t".join(cut)+"\n"
		if diagnosis=="antisense":
			cut[2]="antisense transcript"
			antisense+="\t".join(cut)+"\n"
	return intergenic, antisense


MIN_DISTANCE_TO_GENE=30
ORFs = load_gff_coordinates(sys.argv[1])
intergenic, antisense = load_transcripts_gff(sys.argv[2], ORFs, MIN_DISTANCE_TO_GENE)


f1 = open(sys.argv[3], "w")
f2 = open(sys.argv[4], "w")

f1.write(intergenic)
f2.write(antisense)

f1.close()
f2.close()



