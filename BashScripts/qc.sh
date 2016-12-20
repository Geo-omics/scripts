#!/bin/bash
# USAGE: qc.sh ./Sample_####
# Module required
# module load Scythe/0.993b (or above)

# Bash error check
set -e
# Bash check for uninitialized variables
# http://www.davidpashley.com/articles/writing-robust-shell-scripts/#id2382181
set -u

ADAPTERS=/opt/packages/Trimmomatic/0.32/adapters/TruSeq3-PE-2.fa
DEREPLICATE=/geomicro/data1/COMMON/scripts/SeqTools/dereplicate.pl
INTERLEAVE=/geomicro/data1/COMMON/scripts/SeqTools/interleave.pl

sample=$1

if [ ! -d "$sample" ]; then
    echo "Error: $sample is not a directory."
    exit 1
fi

echo $sample
SNUM=$(basename "$sample")
SNUM=${SNUM#Sample_}

if ! ls "$sample/$SNUM"*_R[12]_*.fastq* >/dev/null; then
    echo "Error: $sample does not seem to contain raw fastq[.gz] files with sample id $SNUM."
    exit 1
fi

# Decompress Raw Data
echo -e "[`date`]\tDecompressing"
if ls "$sample/${SNUM}"*.fastq.gz >/dev/null 2>&1; then
    gunzip -v "$sample/${SNUM}"*.fastq.gz
fi

# Concatenate Raw Data
echo -e "[`date`]\tConcatenating"
FWD=$sample/${SNUM}_fwd.fastq
REV=$sample/${SNUM}_rev.fastq
if [ ! -e $FWD ]; then
	cat $sample/*R1*.fastq > $FWD
fi

if [ ! -e $REV ]; then
	cat $sample/*R2*.fastq > $REV
fi

# Quality Check 1
echo -e "[`date`]\tQC-1 started in background"
mkdir -p $sample/FASTQC
fastqc -o $sample/FASTQC -f fastq -t 2 $FWD $REV &> $sample/$SNUM.fqc.log &

# Dereplication
echo -e "[`date`]\tDereplicating"
perl $DEREPLICATE -fq $FWD -o $sample/derep_fwd_${SNUM}.fastq
perl $DEREPLICATE -fq $REV -o $sample/derep_rev_${SNUM}.fastq

# Scythe Adapter Trimming
echo -e "[`date`]\tAdapter Trimming"
scythe -a $ADAPTERS -q sanger -m $sample/derep_fwd_${SNUM}.matches.fastq -o $sample/derep_scythe_fwd_${SNUM}.fastq $sample/derep_fwd_${SNUM}.fastq
scythe -a $ADAPTERS -q sanger -m $sample/derep_rev_${SNUM}.matches.fastq -o $sample/derep_scythe_rev_${SNUM}.fastq $sample/derep_rev_${SNUM}.fastq

# Sickle Quality Trimming
echo -e "[`date`]\tQuality Trimming"
sickle se -t sanger -f $sample/derep_scythe_fwd_${SNUM}.fastq -o $sample/derep_scythe_sickle_fwd_${SNUM}.fastq
sickle se -t sanger -f $sample/derep_scythe_rev_${SNUM}.fastq -o $sample/derep_scythe_sickle_rev_${SNUM}.fastq

# Quality Check 2
echo -e "[`date`]\tQC-2 started in background"
fastqc -o $sample/FASTQC -f fastq -t 2 $sample/derep_scythe_sickle_*_${SNUM}.fastq &>> $sample/$SNUM.fqc.log &

# Interleave
echo -e "[`date`]\tInterleaving"
perl $INTERLEAVE -fastq -outfmt fasta -rev $sample/derep_scythe_sickle_rev_${SNUM}.fastq -fwd $sample/derep_scythe_sickle_fwd_${SNUM}.fastq -o $sample/dt &

#done
echo "[`date`]\tDone!"

# Notify user.
my_name="$(basename $0)"
echo -e "[`date`]\t${my_name}:\nOUTPUT:${PWD}" >> email.log
mail -s "QC Finished!" ${USER}@umich.edu < email.log
