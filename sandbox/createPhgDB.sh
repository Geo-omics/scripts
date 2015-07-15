#!/bin/bash
set -u
set -e

scripts="/geomicro/data1/COMMON/scripts"

fmtDate=`date +"%m_%Y"`
giFile="phgDB_${fmtDate}.gi.list"
fileBasename=${giFile##*/}
xmlFile=${fileBasename%.*}.xml
rfmtXmlFile=${fileBasename%.*}.refmt.xml
dataFile=${fileBasename%.*}.data

echo "[`date`] Setting up helper scripts..."
for script in ppt_getGI.pl ppt_getXML.pl parseTinySeqXML.xslt derep+alias.pl; do
	if [ ! -h $script ]; then
		ln -s $scripts/PPTTools/$script .
	fi
done

echo "[`date`] Started creation of a new PhgDB from scratch..."

# Use eUtils to fetch Gi numbers for the given search term "Viruses[Organism] AND srcdb_refseq[PROP] NOT cellular organisms[ORGN] NOT AC_000001:AC_999999[pacc] AND gbdiv_phg[properties]"
echo "Searching for phages in NCBI and getting GI numbers matching the search..."
perl ppt_getGI.pl -db protein -o $giFile &> ppt_getGI_${fmtDate}.log 
numGI=$(sort -u $giFile | wc -l)

echo -e "Number of GI's found:\t$numGI"
## 161388 phgDB_07_2015.gi.list

echo "Fetching Protein records in TinySeq XML format for each GI..."
perl ppt_getXML.pl -gi $giFile -db protein -o $xmlFile &> ppt_getXML_${fmtDate}.log

# Using eUtils  messes up the xml output (because it queries your GI numbers 500 records at a time; it concatenates each output; and not in a smart way.) The following is a hack around this; needs to be made more robust.
echo "Reformatting TinySeq XML output..."
head -n 3 $xmlFile > $rfmtXmlFile 
grep -v "<\?xml" $xmlFile | grep -v "<\!DOCTYPE" | grep -v "<TSeqSet>" | grep -v "</TSeqSet>" >> $rfmtXmlFile
tail -n 3 $xmlFile >> $rfmtXmlFile

# Reformatted TinySeq XML to table; for legacy reasons; Format: <GI>\t<NCBI Taxon ID>\t<Prot Sequence>
xsltproc parseTinySeqXML.xslt $rfmtXmlFile > $dataFile

# Sanity checks
echo "Running some numbers..."
## Number of unique GI records extracted
numXmlGI=$(cut -f 1 $dataFile | sort -u | wc -l)
echo -e "\tNumber of unique GI records fetched:\t${numXmlGI}"
## 161388 -- Matches the number in the original list; we're good.

if [ $numGI -ne $numXmlGI ]; then
	echo -e "\t[WARNING] The number of GI numbers searched do not match the number of records downloaded!"
fi

# Number of Unique NCBI Taxon IDs
echo -e "\tNumber of Unique Taxon IDs downloaded:\t$(cut -f 2 $dataFile | sort -u | wc -l)"
# 1633 -- This is up from this time last year(07_2014) == 1429

# Remove duplicate(100%) sequences; assign unique ids and create phgDB in fasta format;
echo "Creating a brand-spankin new phgDB..."
perl derep+alias.pl -l $giFile -d $dataFile -all &> derep+alias_${fmtDate}.log
cat derep+alias.log
# Sample Size:	1633
# Number of Random Taxa:	1633
# Total Sequences:	161388
# Total Unique Sequences:	136329
