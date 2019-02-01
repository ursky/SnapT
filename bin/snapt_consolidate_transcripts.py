#!/usr/bin/env python2
import sys


def load_transcript_status(filename, existing_status=None):
	transcript_status={}
	for line in open(filename):
		if line[0]=="#":
			continue
		cut=line.strip().split("\t")
		if len(cut)!=9:
			continue
		status=cut[2]
		info=cut[8].split(";")
		for field in info:
			if field.strip().startswith("transcript_id"):
				transcript = field.split('"')[1]
				break
		transcript_status[transcript] = status
		if existing_status!=None:
			if transcript in existing_status:
				other_status=existing_status[transcript]
				if other_status=="antisense transcript":
					cut[2]="antisense transcript"
				print "\t".join(cut)
	return transcript_status
		


orf_status = load_transcript_status(sys.argv[1])
anno_status = load_transcript_status(sys.argv[2], existing_status=orf_status)










