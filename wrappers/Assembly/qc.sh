#!/bin/bash
# USAGE: qc.log ./Sample_####
# Module required
# module load Scythe/0.993b (or above)

# Bash error check
set -e
# Bash check for uninitialized variables
# http://www.davidpashley.com/articles/writing-robust-shell-scripts/#id2382181
set -u

if [ ! -h "NexteraPE-PE.fa" ]; then
	ln -s /opt/packages/Trimmomatic/0.32/adapters/NexteraPE-PE.fa .
fi

if [ ! -h "dereplicate.pl" ]; then
	ln -s /geomicro/data1/COMMON/scripts/SeqTools/dereplicate.pl .
fi

if [ ! -h "interleave.pl" ]; then
	ln -s /geomicro/data1/COMMON/scripts/SeqTools/interleave.pl .
fi


sample=$1
echo $sample
SNUM=$(echo $sample | sed "s#./##" |sed "s#Sample_##")

# Decompress Raw Data
echo -e "[`date`]\tDecompressing"
if [ -f $sample/${SNUM}*L007*_R1_004*.gz ]; then
	for fwd in $sample/*R1*.gz; do
		gunzip $fwd
	done
fi

if [ -f $sample/${SNUM}*L007_R2_004*.gz ]; then
	for rev in $sample/*R2*.gz; do
 		gunzip $rev
	done 
fi

# Concatenate Raw Data
echo -e "[`date`]\tConcatenating"
FWD=$sample/${SNUM}_fwd.fastq
REV=$sample/${SNUM}_rev.fastq
cat $sample/*R1*.fastq > $FWD
cat $sample/*R2*.fastq > $REV

# Quality Check 1
echo -e "[`date`]\tQC-1 started in background"
mkdir -p $sample/FASTQC
fastqc -o $sample/FASTQC -f fastq -t 2 $FWD $REV &> $sample/$SNUM.fqc.log &

# Dereplication
echo -e "[`date`]\tDereplicating"
perl dereplicate.pl -fq $FWD -o $sample/derep_fwd_${SNUM}.fastq
perl dereplicate.pl -fq $REV -o $sample/derep_rev_${SNUM}.fastq

# Scythe Adapter Trimming
echo -e "[`date`]\tAdapter Trimming"
scythe -a NexteraPE-PE.fa -q sanger -m $sample/derep_fwd_${SNUM}.matches.fastq -o $sample/derep_scythe_fwd_${SNUM}.fastq $sample/derep_fwd_${SNUM}.fastq
scythe -a NexteraPE-PE.fa -q sanger -m $sample/derep_rev_${SNUM}.matches.fastq -o $sample/derep_scythe_rev_${SNUM}.fastq $sample/derep_rev_${SNUM}.fastq

# Sickle Quality Trimming
echo -e "[`date`]\tQuality Trimming"
sickle se -t sanger -f $sample/derep_scythe_fwd_${SNUM}.fastq -o $sample/derep_scythe_sickle_fwd_${SNUM}.fastq
sickle se -t sanger -f $sample/derep_scythe_rev_${SNUM}.fastq -o $sample/derep_scythe_sickle_rev_${SNUM}.fastq

# Quality Check 2
echo -e "[`date`]\tQC-2 started in background"
fastqc -o $sample/FASTQC -f fastq -t 2 $sample/derep_scythe_sickle_*_${SNUM}.fastq &>> $sample/$SNUM.fqc.log &

# Interleave
echo -e "[`date`]\tInterleaving"
perl interleave.pl -fastq -outfmt fasta -rev $sample/derep_scythe_sickle_rev_${SNUM}.fastq -fwd $sample/derep_scythe_sickle_fwd_${SNUM}.fastq -o $sample/dt &

#done
echo "[`date`]\tDone!"

# Notify user.
my_name="$(basename $0)"
echo -e "[`date`]\t${my_name}:\nOUTPUT:${PWD}" >> email.log
mail -s "QC Finished!" ${USER}@umich.edu < email.log
