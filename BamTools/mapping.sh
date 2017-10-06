#!/bin/bash
# Usage: ./mapping.sh Sample#
# Assumption1: Sample# == Dir name which contains a fwd.fastq, rev.fastq and reference fasta
# Assumption2: Reference Name == scaffold.fasta

set -e
set -u

# Make sure you load the following modules before submitting the scripts
# module load bwa/0.7.9a
# module load samtools/1.0
# module load bedtools/2.22.1

#######################################
##### MAKE PARAMETER CHANGES HERE #####
#######################################
# Only accepted extensions are ".fasta" and ".fastq"

path2scaffolds=$1/scaffold.fasta
path2fwd=$1/fwd.fastq
path2rev=$1/rev.fastq
threads=8

####################################################
##### DO NOT MAKE ANY CHANGES BEYOND THIS LINE #####
#####     unless you know what you're doing    #####
####################################################
version="1.3.3"
alnSam=$(echo $path2scaffolds | sed "s#.fasta#_aligned.sam#")
alnSamLog=${path2scaffolds}.aln.log
sam2bam=$(echo $path2scaffolds | sed "s#.fasta#_fixmate.bam#")
sortBam=$(echo $path2scaffolds | sed "s#.fasta#_sorted.bam#")
readGroup="@RG\\tID:group_${1}\\tSM:Sample_${1}\\tPL:illumina\\tLB:lib_${1}\\tPU:unit_${1}"
bedOut=$(echo $path2scaffolds | sed "s#.fasta#.genomeCovBed.tsv#")

#######################################
######### COMMANDS START HERE #########
#######################################
echo -e "Starting mapping pipeline version: ${version}"
echo -e "[`date`]\tIndexing Database"
bwa index $path2scaffolds
echo -e "[`date`]\tAligning FWD and REV reads using 'bwa mem'"
echo -e "[NOTE]\tUsing Read Group ${readGroup}"
bwa mem -M -R ${readGroup} -t $threads $path2scaffolds $path2fwd $path2rev 1> $alnSam 2> $alnSamLog
echo -e "[`date`]\tFixing alignment artifacts and converting SAM to BAM"
samtools fixmate -O bam $alnSam $sam2bam
echo -e "[`date`]\tSorting BAM"
samtools sort -m 50G -@ ${threads} -O bam -o $sortBam -T $$.samSort.tmp $sam2bam
echo -e "[`date`]\tIndexing Sorted BAM"
samtools index $sortBam
echo "[`date`]\tCalculating Coverage"
echo -e "Contig\tDepth\tbases_with_depth\tcontigLen\tfraction_bases_with_depth" > ${bedOut}
genomeCoverageBed -ibam $sortBam -g $path2scaffolds > $$.bed.tmp
grep -v "^__" $$.bed.tmp >> ${bedOut}
rm $$.bed.tmp
echo -e "[`date`]\tDone!"
