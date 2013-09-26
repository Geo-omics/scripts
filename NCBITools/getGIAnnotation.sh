#!/bin/bash

# usage: sh getGIAnnotation.sh blastOutput database_type
# example: sh getGIAnnotation.sh test.blastn nucleotide
if [ -z $1 ]; then echo "command failed: Give a file name with list"; exit; fi
if [ ! -s $1 ]; then echo "$1 does not exist"; exit; fi
cut -d '|' -f 2 $1 | sort -u > gi.list

perl /geomicro/data1/COMMON/scripts/NCBITools/getGiInfo.pl -d $2 -o anno.xml -l gi.list

perl /geomicro/data1/COMMON/scripts/NCBITools/GI_info_XMLParser.pl anno.xml gi.desc
