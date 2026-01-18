#!/usr/bin/env python3

#This script formats counts and prints codon_position, stop identity and
#total count

#Version 2 reads premRNA, spliced RNA and spike-in categories if present
#in the input file

import sys

file = sys.argv[1]

for line in open(file):
    if "premRNA" in line:
        if "_PTC_" in line:
            fields = line.strip().split(" ")
            count = fields[1]
            name = fields[0]
            name_entries = name.split("_")
            codon_position = name_entries[10]
            stop = name_entries[2]
            print codon_position, stop, count, "premRNA"

        elif "_WT_" in line:
            fields = line.strip().split(" ")
            count = fields[1]
            name = fields[0]
            name_entries = name.split("_")
            codon_position = name_entries[10]
            stop = name_entries[2]
            print "WT", stop, count, "premRNA"

        elif "_Syn_" in line:
            fields = line.strip().split(" ")
            count = fields[1]
            name = fields[0]
            name_entries = name.split("_")
            codon_position = name_entries[10]
            stop = name_entries[2]
            print "Syn", stop, count, "premRNA"

        else: continue

    elif "spliced" in line:
        if "_PTC_" in line:
            fields = line.strip().split(" ")
            count = fields[1]
            name = fields[0]
            name_entries = name.split("_")
            codon_position = name_entries[10]
            stop = name_entries[2]
            print codon_position, stop, count, "spliced"

        elif "_WT_" in line:
            fields = line.strip().split(" ")
            count = fields[1]
            name = fields[0]
            name_entries = name.split("_")
            codon_position = name_entries[10]
            stop = name_entries[2]
            print "WT", stop, count, "spliced"

        elif "_Syn_" in line:
            fields = line.strip().split(" ")
            count = fields[1]
            name = fields[0]
            name_entries = name.split("_")
            codon_position = name_entries[10]
            stop = name_entries[2]
            print "Syn", stop, count, "spliced"

        else: continue


    elif "Spike" in line:
        fields = line.strip().split(" ")
        count = fields[1]
        name = fields[0]
        print name, "Stop_NA", count, "RNA_NA"


    else:
        if "_PTC_" in line:
            fields = line.strip().split(" ")
            count = fields[1]
            name = fields[0]
            name_entries = name.split("_")
            codon_position = name_entries[10]
            stop = name_entries[2]
            print codon_position, stop, count

        elif "_WT_" in line:
            fields = line.strip().split(" ")
            count = fields[1]
            name = fields[0]
            name_entries = name.split("_")
            codon_position = name_entries[10]
            stop = name_entries[2]
            print "WT", stop, count

        elif "_Syn_" in line:
            fields = line.strip().split(" ")
            count = fields[1]
            name = fields[0]
            name_entries = name.split("_")
            codon_position = name_entries[10]
            stop = name_entries[2]
            print "Syn", stop, count

        else: continue



