#!/usr/bin/env python2
import sys
#usage: intersect_gff.py annotation.gff transcripts.gff nc_transcripts.gff

def load_gff_coordinates(filename):
	coordinates={}
	# load coordinates from gff file:
	for line in open(sys.argv[1]):
		if line[0]=="#":
			continue
		cut=line.strip().split("\t")
		if len(cut) != 9:
			continue

		contig = cut[0]
		st = int(cut[3])
		fi = int(cut[4])
		strand = cut[6]
		if strand=="1":
			strand="+"
		if strand=="0":
			strand="-"
		gene = cut[8].split(";")[0]
	
		if contig not in coordinates:
			coordinates[contig]={}
			coordinates[contig]["+"]=[]
			coordinates[contig]["-"]=[]
		if strand==".":
			coordinates[contig]["+"].append((st,fi,gene))
			coordinates[contig]["-"].append((st,fi,gene))
		else:
			coordinates[contig][strand].append((st,fi,gene))
	return coordinates


def intersect(st, fi, strand_coordinates, dist, overlap):
	diagnosis="intergenic"
	coordinate=(0,1,"NA")
	for i,coordinate in enumerate(strand_coordinates):
		# middle of gene
		if st-overlap>=coordinate[0] and fi+overlap<=coordinate[1]:
			diagnosis="gene"
			break
		# start of gene
		elif fi+overlap<=coordinate[1] and fi-overlap>=coordinate[0]:
			diagnosis="gene"
			break
		# end of gene
		elif st+overlap<=coordinate[1] and st-overlap>=coordinate[0]:
			diagnosis="gene"
			break
		# complete span of gene
		elif st<coordinate[0] and fi>coordinate[1]:
			# small peptide encoded in ncRNA
			if coordinate[1]-coordinate[0]<100 and (fi-st)/(coordinate[1]-coordinate[0])>3:
				# besides the peptide, the transcript doesnt intersect with anything on the strand
				if i>0:
					if st-strand_coordinates[i-1][1]>dist:
						diagnosis="intergenic"
					else:
						diagnosis="ambiguous"
				
				if i+1<len(strand_coordinates):
					if strand_coordinates[i+1][0]-fi>dist:
						diagnosis="intergenic"
					else:
						diagnosis="ambiguous"
				if len(strand_coordinates)==1:
					diagnosis="intergenic"
				break
			# likely the peptide or gene itself
			else:
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
	return diagnosis, coordinate[2]



def load_transcripts_gff(filename, orfs, dist, overlap):
	other_strand={"+":"-", "-":"+"}
	intergenic=""
	antisense=""

	for line in open(filename):
		if line[0]=="#":
			continue
		cut=line.strip().split("\t")
		if "transcript" not in cut[2]:
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
		diagnosis,gene = intersect(st, fi, strand_coordinates, dist, 0)

		# if the transcript is intergenic, check if there is a gene on the other strand
		if diagnosis=="intergenic":
			strand_coordinates = orfs[contig][other_strand[strand]]
			opposite_strand_diagnosis,opposite_gene = intersect(st, fi, strand_coordinates, dist, overlap)
			if opposite_strand_diagnosis=="gene":
				diagnosis = "antisense"

		if diagnosis=="intergenic":
			cut[2]="intergenic transcript"
			intergenic+="\t".join(cut)+"\n"
			print "\t".join(cut)
		if diagnosis=="antisense":
			cut[2]="antisense transcript"
			cut[8]=cut[8]+' antisense_to_gene "'+opposite_gene+'";'
			antisense+="\t".join(cut)+"\n"
			print "\t".join(cut)
	return intergenic, antisense


MIN_DISTANCE_TO_GENE=30
MIN_OVERLAP_WITH_GENE=10
ORFs = load_gff_coordinates(sys.argv[1])
intergenic, antisense = load_transcripts_gff(sys.argv[2], ORFs, MIN_DISTANCE_TO_GENE, MIN_OVERLAP_WITH_GENE)



