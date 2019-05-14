#DJS 5 September 2018

#Run this shell script in a directory containing bin fasta files to get a summary tab-delimited file to import the bins as a collection in ANVIO.
#It was originally written for VizBin collections, but will work for any group of Bin fastas generated from any binning program

#Make a list of contigs in each bin, and add the bin file in a column next to the contig:
for i in *.fa; do
	grep ">" $i | sed "s/$/	$i/" > ${i}.list;
done

#Concactenate the data into one list file:
cat *.list > cat.list

#remove the file extension from the bin name:
sed 's/.fa//g' cat.list > cat2.list

#Remove the ">" leftover from fasta headers:
sed 's/>//g' cat2.list > VizBin_binning_results.txt

#delete intermediate files:
rm *.list
