#!/usr/bin/env bash
#BSUB -J Mapping[1-12]
#BSUB -e logs/Mapping.%J.%I.err
#BSUB -o logs/Mapping.%J.%I.out
#BSUB -R "rusage[mem=80] span[hosts=1]"
#BSUB -n 8

#This script calculates mRNA expression levels by 
#diving total counts in the cDNA library by the total counts in the 
#gDNA library for each sequence context in the library.  

#Each value is then further normalized dividing by the normalized value of
#the synonymous context. 

#Relative expression of a sequence context(i) = (cDNAi/gDNAi)/(cDNAsyn/gDNAsyn)
#where cDNA and gDNA are the total counts for the i or syn sequence
#context

mkdir -p logs

#Names of input gDNA files
gDNA_Counts=(gSGE_X10R1_DMSO1.PTC.CountsCondensed
gSGE_X10R1_DMSO1.SNV.CountsCondensed
gSGE_X10R1_DMSO2.PTC.CountsCondensed
gSGE_X10R1_DMSO2.SNV.CountsCondensed
gSGE_X10R1_DMSO3.PTC.CountsCondensed
gSGE_X10R1_DMSO3.SNV.CountsCondensed
gSGE_X10R1_SMG1i1.PTC.CountsCondensed
gSGE_X10R1_SMG1i1.SNV.CountsCondensed
gSGE_X10R1_SMG1i2.PTC.CountsCondensed
gSGE_X10R1_SMG1i2.SNV.CountsCondensed
gSGE_X10R1_SMG1i3.PTC.CountsCondensed
gSGE_X10R1_SMG1i3.SNV.CountsCondensed)

#Names of input cDNA files 
cDNA_Counts=(rSGE_X10R1_DMSO1.PTC.CountsCondensed
rSGE_X10R1_DMSO1.SNV.CountsCondensed
rSGE_X10R1_DMSO2.PTC.CountsCondensed
rSGE_X10R1_DMSO2.SNV.CountsCondensed
rSGE_X10R1_DMSO3.PTC.CountsCondensed
rSGE_X10R1_DMSO3.SNV.CountsCondensed
rSGE_X10R1_SMG1i1.PTC.CountsCondensed
rSGE_X10R1_SMG1i1.SNV.CountsCondensed
rSGE_X10R1_SMG1i2.PTC.CountsCondensed
rSGE_X10R1_SMG1i2.SNV.CountsCondensed
rSGE_X10R1_SMG1i3.PTC.CountsCondensed
rSGE_X10R1_SMG1i3.SNV.CountsCondensed)


gDNA_Count=${gDNA_Counts[$(($LSB_JOBINDEX - 1))]}
cDNA_Count=${cDNA_Counts[$(($LSB_JOBINDEX - 1))]}

gDNA_Files_Location=/beevol/home/osoriom/Jagannathan_Lab/LMNA_SGE_Project/Analysis_of_SGE_Results_LIG4_KO/250804_SGElig4KO_Exons_Exp1_X10/Bowtie
cDNA_Files_Location=/beevol/home/osoriom/Jagannathan_Lab/LMNA_SGE_Project/Analysis_of_SGE_Results_LIG4_KO/250804_SGElig4KO_Exons_Exp1_X10/Bowtie

Programs=/beevol/home/osoriom/Jagannathan_Lab/LMNA_SGE_Project/Analysis_of_SGE_Results/Programs

python $Programs/Calculate_Expression_SynNorm.3.py \
$gDNA_Files_Location/$gDNA_Count \
$cDNA_Files_Location/$cDNA_Count > \
$cDNA_Count.NormExpression 




