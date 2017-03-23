#!/bin/bash

#velveth assembly_61 61 -fastq -short derep_trimmed_day_fwd.fastq
#velvetg assembly_61 -exp_cov auto -ins_length 249
#meta-velvetg assembly_61 -ins_length 249 | tee logfile

module load AMOS/3.1.0
module load velvet/1.1.07-MAX99-OPENMP
module load MetaVelvet/1.0.01

args=("$@")

KMER=${args[0]}
INS=${args[1]}
FASTQ=${args[@]:2} # everything from element 2 onwards, inc element 2

INTERVAL=120  # Change this number to modify the time interval (in seconds) of the usageStats script.

echo "K-mer Length=		$KMER"
echo "Insert Length=		$INS"
echo "FastQ Files=	$FASTQ"

if [ $# -ne 4 ]; then # if num of arguments not equal to 4
	echo "USAGE: ./$0 <k-mer Length> <Insert Length> <Forward Fastq File> <Reverse Fastq File>"
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
if [ -d $OUTDIR ]; then
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
velvetg $OUTDIR -exp_cov auto -ins_length $INS -read_trkg yes >> $LOG
echo >> $LOG
echo "************************************* META-VELVETG **************************************************" >> $LOG
echo >> $LOG
meta-velvetg $OUTDIR -ins_length $INS -amos_file yes -scaffolding yes >> $LOG
echo >> $LOG


exit
: << COMMENT1

REPEATS="$BANK.rep"
PREFIX="paired_scaff_$KMER\_pre"
REDNDCY=5
CONTIG="paired_$KMER.fasta"
SCAFF="paired_scaff_$KMER.fasta"

echo "************************************* Bank-Transact *************************************************" >> $LOG
echo >> $LOG
bank-transact -m $OUTDIR/meta-velvetg.asm.afg -b $BANK -cf >> $LOG
echo >> $LOG
echo "************************************* CLK ***********************************************************" >> $LOG
echo >> $LOG
clk -b $BANK >> $LOG
echo >> $LOG
echo "************************************ Bundler ********************************************************" >> $LOG
echo >> $LOG
Bundler -b $BANK -t M >> $LOG
echo >> $LOG
echo "************************************* MarkRepeats ***************************************************" >> $LOG
echo >> $LOG
MarkRepeats -b $BANK -redundancy $REDNDCY >> $REPEATS
echo >> $LOG
echo "************************************* Archiving *****************************************************" >> $LOG
echo >> $LOG
tar czvf $BANK.tar.gz $BANK
echo >> $LOG
echo "************************************* Orient Contigs ************************************************" >> $LOG
echo >> $LOG
OrientContigs -b $BANK -redundancy $REDNDCY -prefix $PREFIX -repeats $REPEATS -all >> $LOG
echo >> $LOG
echo "************************************* Linearize Scaffolds *******************************************" >> $LOG
echo >> $LOG
untangle -e $PREFIX.evidence.xml -s $PREFIX.out.xml -o $PREFIX.untangle.xml >> $LOG
echo >> $LOG
echo "************************************* Output contig FastA file **************************************" >> $LOG
echo >> $LOG
bank2fasta -d -b $BANK > $CONTIG
echo >> $LOG
echo "************************************* Output scaffold FastA file ************************************" >> $LOG
echo >> $LOG
printScaff -merge -e $PREFIX.evidence.xml -o $PREFIX -s $PREFIX.untangle.xml -l $PREFIX.library -f  >> $LOG
echo >> $LOG
exit

COMMENT1


