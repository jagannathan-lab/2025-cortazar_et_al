#!/usr/bin/env bash
#BSUB -J Mapping[1-6]
#BSUB -e logs/Mapping.%J.%I.err
#BSUB -o logs/Mapping.%J.%I.out
#BSUB -R "rusage[mem=80] span[hosts=1]"
#BSUB -n 8

#This script calculates NMD efficiency from cDNA libraries

#For each sequence context, the total count is divided by the total
#count of the synonymous sequence contex in the same library, in the DMSO
#and in the SMG1i-treated cDNA libraries.

#The normalized value in the DMSO-treated condition is then divided
#by the normalized value in the SMG1i-treated condition. 

#Finally, the product is subtracted from 1. 

#NMD efficiency for a sequence context (j) = 1 - (DMSOj/DMSOsyn)/(SMG1ij/SMG1isyn)
# where DMSO and SMG1i are the total counts in the DMSO or SMG1i treated
# conditions

#Examples
#libraries made from the SMG1i-treated condition
NMDoff_Counts=(rSGE_X10R1_SMG1i1.PTC.CountsCondensed
rSGE_X10R1_SMG1i1.SNV.CountsCondensed
rSGE_X10R1_SMG1i2.PTC.CountsCondensed
rSGE_X10R1_SMG1i2.SNV.CountsCondensed
rSGE_X10R1_SMG1i3.PTC.CountsCondensed
rSGE_X10R1_SMG1i3.SNV.CountsCondensed)

#These are total counts and normalization is also 
#performed within the library first to Syn context
#Then cDNA counts are normalized to the Syn normalized gDNA counts
NMDon_Counts=(rSGE_X10R1_DMSO1.PTC.CountsCondensed
rSGE_X10R1_DMSO1.SNV.CountsCondensed
rSGE_X10R1_DMSO2.PTC.CountsCondensed
rSGE_X10R1_DMSO2.SNV.CountsCondensed
rSGE_X10R1_DMSO3.PTC.CountsCondensed
rSGE_X10R1_DMSO3.SNV.CountsCondensed)

On_Count=${NMDon_Counts[$(($LSB_JOBINDEX - 1))]}
Off_Count=${NMDoff_Counts[$(($LSB_JOBINDEX - 1))]}

Files_Location=/beevol/home/osoriom/Jagannathan_Lab/LMNA_SGE_Project/Analysis_of_SGE_Results_LIG4_KO/250804_SGElig4KO_Exons_Exp1_X10/Bowtie

Programs=/beevol/home/osoriom/Jagannathan_Lab/LMNA_SGE_Project/Analysis_of_SGE_Results/Programs

python $Programs/Calculate_NMDActivity.py \
$Files_Location/$On_Count \
$Files_Location/$Off_Count > \
$Off_Count.NMD 




