#!/bin/bash

# usage: sh getGIAnnotation.sh blastOutput database_type
# example: sh getGIAnnotation.sh test.blastn nucleotide

cut -d '|' -f 2 $1 | sort -u > gi.list

perl /geomicro/data1/COMMON/scripts/getGiInfo.pl -d $2 -o anno.xml -l gi.list

perl /geomicro/data1/COMMON/scripts/GI_info_XMLParser.pl anno.xml gi.desc
