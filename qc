#!/bin/bash
#
# Copyright (c) 2016 Regents of The University of Michigan.
# All Rights Reserved.
#
# Author: Robert <heinro@med.umich.edu>
#

usage() {
    if [ "$VERBOSITY" -ge 1 ]; then
	echo "
USAGE: $SCRIPT_NAME [-h|--help|--no-color|-v|verbosity n]
"
    fi
}

user_help() {
    echo "$SCRIPT_NAME - quality control for a sample's fastq"
    usage
    echo "Command line parameters:
    --working-dir <path>  Where to retrieve and store files
    --clean-only
-h  --help          Print this help.
    --no-color      Disable colorful output.
-v                  Verbosity: use multiple -v to increase output.
    --verbosity <n> Set level verbosity.
"
}

exception () {
    echo "[${RED}ERROR${ENDCOLOR}] at line $1; error code $2"
}

debug() {
    [ "$VERBOSITY" -ge 3 ] && echo "[$SCRIPT_NAME] [DEBUG] $1"
}

info() {
    [ "$VERBOSITY" -ge 1 ] && echo "[$SCRIPT_NAME] $1"
}

warning() {
    [ "$VERBOSITY" -ge 1 ] && echo -e "[$SCRIPT_NAME] [${ORANGE}WARNING${ENDCOLOR}] $1"
}

error() {
    [ "$VERBOSITY" -ge -1 ] && >&2 echo -e "[$SCRIPT_NAME] [${RED}ERROR${ENDCOLOR}] $1"
}

abort() {
    # write an error message and exit.  If the last positional parameter is 'usage' then also print the usage text.
    >&2 echo -e "[$SCRIPT_NAME] [${RED}ABORT${ENDCOLOR}] $1"
    if [ "$VERBOSITY" -ge 0 ]; then
	if [ $# -ge 0 ]; then
	    [ "${!#}" == "usage" ] && usage
	fi
    fi
    exit 1
}


# get name of this script
SCRIPT_NAME=$(basename "$0")

##########################
# default variable values
##########################

#
ADAPTERS=/opt/packages/Trimmomatic/0.32/adapters/TruSeq3-PE-2.fa
#
CLEAN_ONLY=false
#
DEREPLICATE="perl /geomicro/data1/COMMON/scripts/SeqTools/dereplicate.pl"
#DEREPLICATE=/geomicro/data2/heinro/scripts/SeqTools/dereplicate.pl
#
INTERLEAVE=/geomicro/data1/COMMON/scripts/SeqTools/interleave.pl
# no lane specified by default
WORK_DIR=$(pwd)
# default verbosity level
VERBOSITY=1
# Enable colorful output by default
USE_COLOR=true

###################################3
# command line parameter handling
###################################
GETOPT_SHORT=hv
GETOPT_LONG=clean-only,help,no-color,working-dir:,verbosity:
if which getopt >/dev/null 2>&1; then
    getopt -T
    if [ "$?" == 4 ]; then
        # GNU getopt available 
        _getopt=$(getopt -o "$GETOPT_SHORT" --long "$GETOPT_LONG" --name "$SCRIPT_NAME" -- "$@")
        if [ "$?" != 0 ]; then usage; exit 1; fi
    else
	# for non-GNU getopt (on MacOSX?)
	# try best effort without handling long options
	# TODO: this path needs testing
        _getopt=$(getopt -o "$GETOPT_SHORT" --name "$SCRIPT_NAME" -- "$@")
        if [ "$?" != 0 ]; then usage; exit 1; fi
    fi
    # reset $n parameters
    eval set -- "$_getopt"
fi

# handle options arguments / should work with/without getopt
while [ "$#" -gt 0 ]; do
    case "$1" in
	(--clean-only)
	    CLEAN_ONLY=true
	    ;;
	(--working-dir)
	    WORK_DIR="$2"
            shift
	    ;;
	(-h|--help) user_help; exit 0;;
	(--no-color) USE_COLOR=false;;
	(-v) VERBOSITY=$((VERBOSITY+1));;
	(--verbosity)
	    VERBOSITY="$2"
	    [[ "$VERBOSITY" =~ ^[[:digit:]]+$ ]] || abort "Verbosity level must be a positive integer." usage 
	    shift
	    ;;
	(--) shift; break;;
	(-|-*) abort "invalid command line option: $1" usage;;
	(*) break;; # getopt failure? assume this and following are non-option args
    esac
    shift
done

##########################################################
# color output
#
# for color codes see:
# https://en.wikipedia.org/wiki/ANSI_escape_code#Colors
#
# 
if [ "$USE_COLOR" == true ]; then
    if which tput >/dev/null 2>&1 && [ "$(tput colors)" -ge 8 ]; then
	RED=$(tput setaf 1)
	ORANGE=$(tput bold)$(tput setaf 1)
	ENDCOLOR=$(tput sgr0)
    fi
fi

debug "command line options: $_getopt"
debug "verbosity set to: $VERBOSITY"
# set bash script debugging
[ "$VERBOSITY" -ge 4 ] && set -x
# forward verbosity to other commands
[ "$VERBOSITY" -ge 2 ] && VERBOSE_FLAG=-v
V=$VERBOSE_FLAG

# handle non-option parameters
if [ "$#" -gt 0 ]; then
    abort "Unknown parameters: ${*}"
fi

#########################
# some input sanitation
#########################
[ -d "$WORK_DIR" ] || abort "$WORK_DIR is not accessible."

####################################
# check presence of necessary tools
####################################
for i in fastqc scythe sickle; do
    which $i > /dev/null 2>&1 || abort "$i command is not available."
done

###############################
# enable error trapping
###############################
set -e
trap 'exception $LINENO $?' ERR

#################################
# do stuff
################################
cd "$WORK_DIR"

# clean up
rm $V -rf FASTQC
rm $V -f derep* 
rm $V -f dt_{int,fwd,rev}.fasta
rm $V -f fwd.clust{.list,}
if $CLEAN_ONLY; then
    exit 0
fi


mkdir $V FASTQC
info "pre-QC fastqc run..."
fastqc -o FASTQC -f fastq -t 2 fwd.fastq rev.fastq &> FASTQC/fastqc.log &
info "dereplicating..."
$DEREPLICATE -fq fwd.fastq -o derep_fwd.fastq &
$DEREPLICATE -fq fwd.fastq -o derep_rev.fastq
info "adapter trimming..."
scythe -a $ADAPTERS -q sanger -m derep_fwd.matches.fastq -o derep_scythe_fwd.fastq derep_fwd.fastq &
scythe -a $ADAPTERS -q sanger -m derep_rev.matches.fastq -o derep_scythe_rev.fastq derep_rev.fastq
info "quality score trimming..."
sickle se -t sanger -f derep_scythe_fwd.fastq -o derep_scythe_sickle_fwd.fastq &
sickle se -t sanger -f derep_scythe_rev.fastq -o derep_scythe_sickle_rev.fastq
info "post-QC fastqc run..."
fastqc -o FASTQC -f fastq -t 2 derep_scythe_sickle_*_fastq &>> FASTQC/fastqc.log &
info "interleaving..."
$INTERLEAVE -fastq -outfmt fasta -fwd derep_scythe_sickle_fwd.fastq -rev derep_scythe_sickle_rev.fastq -o dt
info "Done!"


# gmb: input: fwd.fastq
# gmb: input: rev.fastq
# gmb: output: dt_int.fasta
# gmb: cpus: 2
