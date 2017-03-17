#!/bin/bash

#velveth assembly_61 61 -fastq -short derep_trimmed_day_fwd.fastq
#velvetg assembly_61 -read_trkg yes
#oases ~/Velvet_Assembly &


module load AMOS/3.1.0-rc1
module load velvet/1.1.07-MAX99-OPENMP
module load oases/0.2.01

args=("$@")

KMER=${args[0]}
FASTQ=${args[@]:1} # everything from element 2 onwards, inc element 2

INTERVAL=120  # Change this number to modify the time interval (in seconds) of the usageStats script.

echo "K-mer Length=		$KMER"
echo "FastQ Files=	$FASTQ"

if [ $# -ne 3 ]; then # if num of arguments not equal to 4
	echo "USAGE: ./$0 <k-mer Length> <Forward Fastq File> <Reverse Fastq File>"
	echo "Yes, it has to be in that exact order!"
	exit
fi

OUTDIR="assembly_paired_$KMER"
# bash check if directory exists
if [ -d $OUTDIR ]; then
	echo "$OUTDIR already exists!"
	exit
fi 

BANK="bank_paired_$KMER"
if [ -d $BANK ]; then
	echo "$BANK already exists!"
	exit
fi 
LOG="$OUTDIR.log"
STATS="usageStats_K$KMER.tsv"

perl /geomicro/data1/COMMON/scripts/usageStats.pl -i $INTERVAL -o $STATS -e &

echo "************************************* VELVETH *******************************************************" > $LOG
echo >> $LOG
velveth $OUTDIR $KMER -fastq -shortPaired $FASTQ >> $LOG
echo >> $LOG
echo "************************************* VELVETG *******************************************************" >> $LOG
echo >> $LOG
velvetg $OUTDIR -read_trkg yes >> $LOG
echo >> $LOG
echo "************************************* OASES **************************************************" >> $LOG
echo >> $LOG
oases $OUTDIR >> $LOG
echo >> $LOG


exit
