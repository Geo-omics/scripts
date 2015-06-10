#!/usr/bin/perl

=head1 Usage

	# In multiple blast outputs, count the number of times a query gets a hit above a certain bit score.
	perl tallyWrap.pl -ext blastp -m masterList_output -t combinedTally_output -s 40
	
	# In multiple tab delimited files, for each value in the first column, combine the values of the last column from each file into a combined tabular file.
	perl tallyWrap.pl -ext blastp -m masterList_output -t combinedTally_output -values
	
=head1 Dependencies

	getMasterList.pl -- in TabTools/Tally_Compare/
	tally.pl -- in TabTools/Tally_Compare/
	weave.pl -- in TabTools/Tally_Compare/	

=cut


use strict;
use Getopt::Long;
use File::Spec::Functions;

my $masterList;
my $ext="blastn";
my $combinedTally;
my $bs=0;
my $printValue; # column number (start counting from 1)
my $SCRIPTS="/geomicro/data1/COMMON/scripts/TabTools/Tally_Compare/";
my $fasta; # query fasta
GetOptions(
	'ext:s'=>\$ext,
	'm:s'=>\$masterList,
	't:s'=>\$combinedTally,
	's:f'=>\$bs,
	'value|values:i'=>\$printValue,
	'f|fasta:s'=>\$fasta,
	'scripts:s'=>\$SCRIPTS,
);

#Check for all scripts
my $masterListScript=catfile($SCRIPTS, "getMasterList.pl");
my $tallyScript=catfile($SCRIPTS, "tally.pl");
my $weaveScript=catfile($SCRIPTS, "weave.pl");

die "Could not find dependent scripts. Use the '-scripts' flag to provide a location to the dependencies" if ((! -e $masterListScript)||(! -e $tallyScript)||(! -e $weaveScript));


my @files=glob("*.$ext");
print "Creating MasterList\n";
system("perl $masterListScript -o ".$masterList." -s ".$bs.($ext ? " -e $ext" : ""));
foreach my $f(@files){
	print "\tTally.pl: $f\n";
	my($name, $ext)=split(/\./, $f);
	my $tallyFile=$name.".tally";
	system("perl $tallyScript -m ".$masterList." -i ".$f ." -o ".$tallyFile." -s ".$bs.($printValue ? " -values $printValue" : "").($fasta ? " -fasta $fasta":""));
}
print "Weaving all tally files...\n";
system("perl $weaveScript -o ".$combinedTally);
exit;
