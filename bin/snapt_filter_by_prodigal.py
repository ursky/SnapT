#!/usr/bin/env python2
import sys


# load prodigal annotation of transcripts
blacklist={}
for line in open(sys.argv[1]):
	if line[0]=="#":
		continue
	cut = line.strip().split("\t")
	# important to look at only the relevant strand (the strand of the ncRNA)
	if cut[6]!="+":
		continue
	transcript = cut[0].split("(")[0]
	positions = transcript.split(":")[1].split("-")
	rna_length = int(positions[1])-int(positions[0])
	orf_length = int(cut[4]) - int(cut[3])
	if orf_length > rna_length/2:
		blacklist[transcript]=None

# only print non-blacklisted transcripts
for line in open(sys.argv[2]):
	if line[0]=="#":
		print line.strip()
		continue
	cut=line.strip().split("\t")
	transcript=cut[0]+":"+str(int(cut[3])-1)+"-"+cut[4]
	if transcript not in blacklist:
		print line.strip()


	
