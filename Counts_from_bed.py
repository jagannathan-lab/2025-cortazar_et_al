#!/usr/bin/env python3

#This script counts the total number of reads per sequence context in the
#library

import sys

file = sys.argv[1]

#dictionary to store total counts
twist_sequences = dict()


for line in open(file):
    fields = line.strip().split("\t")
    entry = fields[0]
    if entry in twist_sequences:
        twist_sequences[entry] += 1
    else:
        twist_sequences[entry] = int(1)

for key in twist_sequences:
    print key,  twist_sequences[key]





