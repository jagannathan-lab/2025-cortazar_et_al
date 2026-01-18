#!/usr/bin/env bash
#BSUB -J SeqPrep[1-13]
#BSUB -e logs/SeqPrep.%J.%I.err
#BSUB -o logs/SeqPrep.%J.%I.out
#BSUB -R "rusage[mem=120] span[hosts=1]"
#BSUB -n 8

mkdir -p logs

#List of fastq base names to process with SeqPrep
NAMES=(gSGE_X10R1_DMSO1_CKDL250021584-1A_22W75YLT4_L3_
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

name=${NAMES[$(($LSB_JOBINDEX - 1))]}

#Location of fastq files
data="/beevol/home/osoriom/Jagannathan_Lab/Sequencing_Results/250804_SGElig4KO_Exons_Exp1_X10/01.RawData"

#Location of SeqPrep
SeqPrep=/beevol/home/osoriom/Jagannathan_Lab/LMNA_SGE_Project/Analysis_of_SGE_Results/SeqPrep

Read1=$data/$name\1.fq.gz
Read2=$data/$name\2.fq.gz

outputR1=$name.R1
outputR2=$name.R2
merge=$name.merge

$SeqPrep/SeqPrep -f $Read1 -r $Read2 -1 $outputR1 -2 $outputR2 \
 -q 0 -L 150 -A AGATCGGAAGAGCACACGTC -B AGATCGGAAGAGCGTCGT \
  -s $merge -o 150 -m 0.001 -n 1


#Report total number of lines in fastq and merge files
echo "Total number of lines in fastq file" > $name.numbers.txt
zcat $Read1 | wc -l >> $name.numbers.txt
echo "Total number of lines in merge file" >> $name.numbers.txt
zcat $merge | wc -l >> $name.numbers.txt

#Remove R1 and R2 files
rm *$name.R1
rm *$name.R2




