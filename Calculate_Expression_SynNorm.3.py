#!/usr/bin/env python3
import sys

#gDNA counts file
gCounts = sys.argv[1]
#cDNA counts file
rCounts = sys.argv[2]


#Store normalized counts (cDNA/gDNA) in norm_counts_list
norm_counts_list = []

#Store noramlized count value for the synonymous context (cDNAsyn/gDNAsyn)
norm_count_syn = ""

for line in open(gCounts):
    gFields = line.strip().split(" ")
    gCodon = gFields[0]
    gMutation = gFields[1]
    gCount = float(gFields[2])

    for entry in open(rCounts):
        rFields = entry.strip().split(" ")
        rCodon = rFields[0]
        rMutation = rFields[1]
        rCount = float(rFields[2])
        #for the same context in gDNA and cDNA libraries
        #calculate cDNA/gDNA
        if gCodon == rCodon and gMutation == rMutation:
            NormCount = rCount/gCount
            context = " ".join([rCodon, rMutation, str(NormCount)])
            norm_counts_list.append(context)
            
            if rCodon == "Syn":
                norm_count_syn = NormCount
            
            else: continue



        else: continue


#Normalize each (cDNA/gDNA) value to the synonymous context value
#(cDNAsyn/gDNAsyn)

for entry in norm_counts_list:
    fields = entry.split(" ")
    Codon = fields[0]
    Mutation = fields[1]
    NormCount = float(fields[2])
    Rel_expression = NormCount/norm_count_syn
    print " ".join([Codon, Mutation, str(Rel_expression)])














