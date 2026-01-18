import sys

DMSO_Counts = sys.argv[1] #DMSO Counts
SMG1i_Counts = sys.argv[2] #SMG1i Counts

#Variable to store the total count of the syn contex
DMSO_SynCount = 0
SMG1i_SynCount = 0

#find the total count of Syn integrations
for line in open(DMSO_Counts):
    DMSO_Fields = line.strip().split(" ")
    DMSO_Codon = DMSO_Fields[0]
    DMSO_Mutation = DMSO_Fields[1]
    DMSO_Count = DMSO_Fields[2]
    
    if DMSO_Codon == "Syn":
        DMSO_SynCount = float(DMSO_Count)

    else: continue


for line in open(SMG1i_Counts):
    SMG1i_Fields = line.strip().split(" ")
    SMG1i_Codon = SMG1i_Fields[0]
    SMG1i_Mutation = SMG1i_Fields[1]
    SMG1i_Count = SMG1i_Fields[2]

    if SMG1i_Codon == "Syn":
        SMG1i_SynCount = float(SMG1i_Count)

    else: continue


#For each sequence context (j), normalize each count in the library to the Syn count (DMSOj/DMSOsyn)
for line in open(DMSO_Counts):
    DMSO_Fields = line.strip().split(" ")
    DMSO_Codon = DMSO_Fields[0]
    DMSO_Mutation = DMSO_Fields[1]
    DMSO_Count = DMSO_Fields[2]
    DMSO_CountNorm = float(DMSO_Count)/DMSO_SynCount

    for entry in open(SMG1i_Counts):
        SMG1i_Fields = entry.strip().split(" ")
        SMG1i_Codon = SMG1i_Fields[0]
        SMG1i_Mutation = SMG1i_Fields[1]
        SMG1i_Count = SMG1i_Fields[2]
        SMG1i_CountNorm = float(SMG1i_Count)/SMG1i_SynCount
        #For the same context in the DMSO and SMG1i conditions, calculate
        #NMD efficiency
        if DMSO_Codon == SMG1i_Codon and DMSO_Mutation == SMG1i_Mutation:
            DMSOoverSMG1i = DMSO_CountNorm/SMG1i_CountNorm 
            NMDActivity = 1 - DMSOoverSMG1i
            print " ".join([SMG1i_Codon, SMG1i_Mutation, str(NMDActivity)])

        else: continue














