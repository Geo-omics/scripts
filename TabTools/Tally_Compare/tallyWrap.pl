#!/user/bin/perl
use strict;
use Getopt::Long;

my $masterList;
my $combinedTally;
my $bs;
GetOptions(
	'm:s'=>\$masterList,
	't:s'=>\$combinedTally,
	's:f'=>\$bs,
);

my @files=glob("*.out");
print @files." Files Found.\n";
print "Creating MasterList\n";
system("perl getMasterList.pl -o ".$masterList." -s ".$bs);
foreach my $f(@files){
	print "\tTally.pl: $f\n";
	my($name, $ext)=split(/\./, $f);
	my $tallyFile=$name.".tally";
	system("perl tally.pl -m ".$masterList." -i ".$f ." -o ".$tallyFile." -s ".$bs);
}
print "Weaving all tally files...\n";
system("perl weave.pl -o ".$combinedTally);
exit;
