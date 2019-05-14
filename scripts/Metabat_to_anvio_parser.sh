#Parse metabat output into format that can be import as a collection into the anvio profile:

#This command will add the Bin ID in a column after the split name
for i in Bin.*; do sed -i "s/$/\t$i/" $i; done
#This command concatenates all the files into one binning results file for anvio:
cat Bin.* > Metabat_binning_results.txt
#Replace the dots with underscores to make anvio happy...
sed -i 's/Bin./Bin_/g' Metabat_binning_results.txt