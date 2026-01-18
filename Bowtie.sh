#!/usr/bin/env bash
#BSUB -J MapR[1-13]
#BSUB -e logsR/MapR.%J.%I.err
#BSUB -o logsR/MapR.%J.%I.out
#BSUB -R "rusage[mem=80] span[hosts=1]"
#BSUB -n 8

#Mapping of genomic sequences

mkdir -p logsR

#base name of SeqPrep merged files
FILENAMES=(gSGE_X10R1_DMSO1_CKDL250021584-1A_22W75YLT4_L3_
gSGE_X10R1_DMSO2_CKDL250021584-1A_22W75YLT4_L3_
gSGE_X10R1_DMSO3_CKDL250021584-1A_22W75YLT4_L3_
gSGE_X10R1_NCtrl_CKDL250021584-1A_22W75YLT4_L3_
gSGE_X10R1_SMG1i1_CKDL250021584-1A_22W75YLT4_L3_
gSGE_X10R1_SMG1i2_CKDL250021584-1A_22W75YLT4_L3_
gSGE_X10R1_SMG1i3_CKDL250021584-1A_22W75YLT4_L3_
rSGE_X10R1_DMSO1_CKDL250021584-1A_22W75YLT4_L3_
rSGE_X10R1_DMSO2_CKDL250021584-1A_22W75YLT4_L3_
rSGE_X10R1_DMSO3_CKDL250021584-1A_22W75YLT4_L3_
rSGE_X10R1_SMG1i1_CKDL250021584-1A_22W75YLT4_L3_
rSGE_X10R1_SMG1i2_CKDL250021584-1A_22W75YLT4_L3_
rSGE_X10R1_SMG1i3_CKDL250021584-1A_22W75YLT4_L3_)

#New names 
NEW_NAMES=(gSGE_X10R1_DMSO1
gSGE_X10R1_DMSO2
gSGE_X10R1_DMSO3
gSGE_X10R1_NCtrl
gSGE_X10R1_SMG1i1
gSGE_X10R1_SMG1i2
gSGE_X10R1_SMG1i3
rSGE_X10R1_DMSO1
rSGE_X10R1_DMSO2
rSGE_X10R1_DMSO3
rSGE_X10R1_SMG1i1
rSGE_X10R1_SMG1i2
rSGE_X10R1_SMG1i3)

#Reference type "gDNA" or "cDNA"
TYPES=(gDNA
gDNA
gDNA
gDNA
gDNA
gDNA
gDNA
cDNA
cDNA
cDNA
cDNA
cDNA
cDNA)

#Reference name
references=(gX10R1
gX10R1
gX10R1
gX10R1
gX10R1
gX10R1
gX10R1
rX10R1
rX10R1
rX10R1
rX10R1
rX10R1
rX10R1)


Name=${FILENAMES[$(($LSB_JOBINDEX - 1))]}
New_Name=${NEW_NAMES[$(($LSB_JOBINDEX - 1))]}
Type=${TYPES[$(($LSB_JOBINDEX - 1))]}
reference=${references[$(($LSB_JOBINDEX - 1))]}


RefLocation=/beevol/home/osoriom/Jagannathan_Lab/LMNA_SGE_Project/Reference_Files_For_Alignment/Reference/PlusEndogenous/Bowtie/SGE2Refs
MergedFileLocation=/beevol/home/osoriom/Jagannathan_Lab/LMNA_SGE_Project/Analysis_of_SGE_Results_LIG4_KO/250804_SGElig4KO_Exons_Exp1_X10/SeqPrep

Reference=$RefLocation/$Type/$reference/$reference
MergedFile=$MergedFileLocation/$Name.merge

#Bowtie mapping
zcat $MergedFile | \
bowtie \
    -p 8 \
    --sam \
    --un $Name.unmapped \
    -v 0 \
    $Reference \
    - \
    2> $Name.log \
    | samtools view -ShuF4 - \
    | samtools sort - -T $New_Name.temp \
    > $New_Name.bam


bamToBed -i $New_Name.bam \
    | sort -k1,1 -k2,2n \
    > $New_Name.bed


# Count total number of alignmets per mutation in bed file
python Counts_from_bed.py \
$New_Name.bed > $New_Name.Counts

#Format the Counts file into a condensed file
python Format_PTC_Counts_CondensedForm.2.py $New_Name.Counts > \
$New_Name.PTC.CountsCondensed

python Format_SNV_Counts_CondensedForm.2.py $New_Name.Counts > \
$New_Name.SNV.CountsCondensed







