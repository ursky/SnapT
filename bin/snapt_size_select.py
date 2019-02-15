#!/usr/bin/env python2
import sys

filename=sys.argv[1]
min_length=int(sys.argv[2])
max_length=int(sys.argv[3])

for line in open(filename):
	if line[0]=="#":
		print line.strip()
		continue
	cut=line.strip().split("\t")
	if len(cut)!=9:
		continue
	length = int(cut[4])-int(cut[3])
	if length<max_length and length>min_length:
		print line.strip()

