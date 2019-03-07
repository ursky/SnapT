#!/usr/bin/env python2
import sys

# load transcripts with rfam hits
skip={}
for line in open(sys.argv[1]):
	if line[0]=="#":
		continue
	cut=line.strip().split()
	if len(cut)<15:
		continue
	if "tRNA" in line or "RNase P" in line:
		skip[cut[2].split("(")[0]]=None

inter_skipped=0
inter_passed=0
anti_skipped=0
anti_passed=0

# only print out gff lines without hits
for line in open(sys.argv[2]):
	line=line.strip()
	if line[0]=="#":
		print line
	cut=line.split("\t")
	if len(cut)!=9:
		continue
	name=cut[0] + ":" + str(int(cut[3])-1) + "-" + cut[4]
	if name not in skip:
		print line
		if cut[2]=="antisense transcript":
			anti_passed+=1
		if cut[2]=="intergenic transcript":
			inter_passed+=1
	else:
		if cut[2]=="antisense transcript":
			anti_skipped+=1
		if cut[2]=="intergenic transcript":
			inter_skipped+=1

sys.stderr.write("Removed "+str(anti_skipped) + " antisense transcripts\n")
sys.stderr.write("Removed "+str(inter_skipped) + " intergenic transcripts\n")
sys.stderr.write("Kept "+str(anti_passed) + " antisense transcripts\n")
sys.stderr.write("Kept "+str(inter_passed) + " intergenic transcripts\n")
	





