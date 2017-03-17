#!/bin/bash
set -e
set -u
#######################################
##### MAKE PARAMETER CHANGES HERE #####
#######################################

# Location of the scripts folder.
scripts="/geomicro/data1/COMMON/scripts"

####################################################
##### DO NOT MAKE ANY CHANGES BEYOND THIS LINE #####
#####     unless you know what you're doing    #####
####################################################
version=1.0.0
fmtDate=`date +"%m_%Y"`
giFile="phgDB_${fmtDate}.gi.list"
fileBasename=${giFile##*/}
xmlFile=${fileBasename%.*}.xml
rfmtXmlFile=${fileBasename%.*}.refmt.xml
dataFile=${fileBasename%.*}.data
phgDB=${fileBasename%.*.*}.db

echo "[`date`] Firing up version: ${version}..."
echo "[`date`] Setting up helper scripts..."
for script in ppt_getGI.pl ppt_getXML.pl parseTinySeqXML.xslt derep+alias.pl; do
	if [ ! -h $script ]; then
		ln -s $scripts/PPTTools/$script .
	fi
done

echo "[`date`] Started creation of a new PhgDB from scratch..."

if [ ! -s $giFile ]; then
	# Use eUtils to fetch Gi numbers for the given search term "Viruses[Organism] AND srcdb_refseq[PROP] NOT cellular organisms[ORGN] NOT AC_000001:AC_999999[pacc] AND gbdiv_phg[properties]"
	echo "Searching for phages in NCBI and getting GI numbers matching the search..."
	perl ppt_getGI.pl -db protein -o $giFile &> ppt_getGI_${fmtDate}.log 
fi

numGI=$(sort -u $giFile | wc -l)
echo -e "Number of GI's found:\t$numGI"
## 161388 phgDB_07_2015.gi.list

if [ ! -s $xmlFile ]; then
	echo "Fetching Protein records in TinySeq XML format for each GI..."
	perl ppt_getXML.pl -gi $giFile -db protein -o $xmlFile &> ppt_getXML_${fmtDate}.log
fi

echo "[`date`] Checking for errors in the XML file..."
# check if there are any errors in the XML.
numErrors=$(grep -c "<ERROR>" $xmlFile) || true
if [ "$numErrors" -gt 0 ]; then
	echo -e "\t[ERROR] Found $numErrors in ${xmlFile}. Exiting..."
	exit
else
	echo -e "\tNo errors found."
fi

# Fixed in Version: 1.0.0; by adding "$params{batch} = -1" in the getXML script.
# Using eUtils  messes up the xml output (because it queries your GI numbers 500 records at a time; it concatenates each output; and not in a smart way.) The following is a hack around this; needs to be made more robust.
# echo "[`date`] Reformatting TinySeq XML output..."
# head -n 3 $xmlFile > $rfmtXmlFile 
# grep -v "<\?xml" $xmlFile | grep -v "<\!DOCTYPE" | grep -v "<TSeqSet>" | grep -v "</TSeqSet>" >> $rfmtXmlFile
# tail -n 3 $xmlFile >> $rfmtXmlFile
$rfmtXmlFile=$xmlFile # since no need to reformat after fix.

# Reformatted TinySeq XML to table; for legacy reasons; Format: <GI>\t<NCBI Taxon ID>\t<Prot Sequence>
xsltproc parseTinySeqXML.xslt $rfmtXmlFile > $dataFile

# Sanity checks
echo "[`date`] Running some numbers..."
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
echo "[`date`] Creating a brand-spankin new phgDB..."
perl derep+alias.pl -l $giFile -d $dataFile -prefix $phgDB -all &> derep+alias_${fmtDate}.log
cat derep+alias_${fmtDate}.log
# Sample Size:	1633
# Number of Random Taxa:	1633
# Total Sequences:	161388
# Total Unique Sequences:	136329
echo "[`date`] All Done!"
