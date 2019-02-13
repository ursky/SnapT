#!/usr/bin/env python


import sys


f=sys.argv[1]

g=sys.argv[2]


keywords = set()
with open(f) as list_file:
    for line in list_file:
        if line.strip():
            keywords.add(line.strip())

with open(g) as master_file:
	for line in master_file:
		if not set(line.split()[:-1]) & keywords:
			print line.rstrip()

        