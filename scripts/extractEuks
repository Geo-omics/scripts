ln -s /geomicro/data1/COMMON/scripts/sandbox/legacy_consolidateJGIdata.pl
perl legacy_consolidateJGIdata.pl -DIR . -OUTDIR consolidated
awk < consolidated/Unclassified.tsv -F'\t' '{ print $6, t, $2 }' > locus_contig.list
awk < *phylodist -F'\t' '{ print $1, "\t", $5 }' | cut -f 1 -d ";" | grep "Eukaryota" | cut -f 1 > eukaryota.list
fgrep -f eukaryota.list locus_contig.list | cut -f 2 | sort -u | sed "s# ##" >  eukaryota_contigs.list
perl /geomicro/data1/COMMON/scripts/SeqTools/extractSeqs.pl -e -l eukaryota_contigs.list -f *.fna -o euksRemoved.fasta
