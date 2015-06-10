#/bin/bash

# USAGE: ./assemble.sh
# Make sure you're running this from the same directory
# that you used for "qc.sh".
# This script assumes that the QC was done using "qc.sh".
# If not, please DO NOT use this script.

set -e
set -u

# Load the following modules before you start the script
# module load idba/1.1.1 (or above)
# module load QUAST/2.3 (or above)
# module load blast/2.2.28 (or above)
# module load PhyloSift/1.0.1 (or above)

#######################################
##### MAKE PARAMETER CHANGES HERE #####
#######################################
# Required Databases
path2silva="/omics/PublicDB/silva/release_119/SILVA_119_SSURef_tax_silva.fasta"
path2bact="/omics/PublicDB/NCBI/Bacteria/ncbi_bacteria_10302014.fasta"
path2arch="/omics/PublicDB/NCBI/Archaea/ncbi_archaea_10302014.fasta"
path2markers="~/share/phylosift"

# IDBA Parameters
mink=52
maxk=92
step=8
assembly=k${mink}_${maxk}_s${step}

# Processors Available
threads=15

####################################################
##### DO NOT MAKE ANY CHANGES BEYOND THIS LINE #####
#####     unless you know what you're doing    #####
####################################################
for i in $(find -maxdepth 1 -type d -name "Sample*"); do
    myPath=${i}/Assembly
    mkdir -p $myPath
    cd $myPath
    ln -s ../*int*fasta .

    echo -e "[`date`]\tAssembling ${i} at $PWD"
    idba_ud -o ${assembly} -r *_int.fasta --num_threads ${threads} --mink $mink --maxk $maxk --step $step &> idba_${assembly}.log

    echo -e "[`date`]\tProducing assembly stats for ${i}"
    quast.py -f --meta -T ${threads} -l "Scaffolds, Contigs" ${assembly}/scaffold.fa ${assembly}/contig.fa &> quast.log

    echo -e "[`date`]\tSearch for scaffolds with 16S"
    mkdir -p BLASTN
    blastn -query ${assembly}/scaffold.fa -db $path2silva -outfmt "7 std qlen stitle" -out BLASTN/${i}_vs_silvaSSU119.blastn -num_threads ${threads}
    ln -s /geomicro/data1/COMMON/scripts/BlastTools/top5.pl .
    T1BLAST=BLASTN/${i}_vs_silvaSSU119.topHits.blastn
    perl top5.pl -t 1 -b BLASTN/${i}_vs_silvaSSU119.blastn -o $T1BLAST

    echo -e "[`date`]\tLook up top hits from 16S search for complete genomes in NCBI"
    BBLAST=BLASTN/${i}_subSeq_vs_bactNCBI.blastn
    ABLAST=BLASTN/${i}_subSeq_vs_archaeaNCBI.blastn
    FASTA=${assembly}/scaffold.fa
    SSEQ=BLASTN/silvaSSU119.topHits.fasta

    ln -s /geomicro/data1/COMMON/scripts/SeqTools/extractSubSeq.pl .
    perl extractSubSeq.pl -query -blast $T1BLAST -f $FASTA -o $SSEQ

    blastn -query $SSEQ -db $path2bact -outfmt "7 std qlen qcovs stitle" -out $BBLAST -num_threads ${threads}
    blastn -query $SSEQ -db $path2arch -outfmt "7 std qlen qcovs stitle" -out $ABLAST -num_threads ${threads}

    TBBLAST=$(echo $BBLAST | sed "s#.blastn#.topHits.blastn#")
    TABLAST=$(echo $ABLAST | sed "s#.blastn#.topHits.blastn#")

    perl top5.pl -t 1 -b $BBLAST -o $TBBLAST
    perl top5.pl -t 1 -b $ABLAST -o $TABLAST

    mkdir -p PhyloSift
    if [ -d $path2markers ]; then
	phylosift all --disable_updates --output PhyloSift/Whole_Assembly --threads 10 $FASTA
    else
	phylosift all --output PhyloSift/Whole_Assembly --threads 10 $FASTA
    fi

    echo -e "[`date`]\tFinished with ${i}"
    echo ""
    cd -
done

echo "[`date`]\tDone!"

# Notify user.
my_name="$(basename $0)"
echo -e "[`date`]\t${my_name}:\nOUTPUT:${PWD}" >> email.log
mail -s "Assembly Finished!" ${USER}@umich.edu < email.log
