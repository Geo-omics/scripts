#! usr/bin/perl

use strict;

sub checkForCompleteness{
	my $fName=shift;
	chomp($fName);
	open (CONTIGS, $fName) || die "Couldn't open $fName\n";
	$/= ">";
	my %sequences;
	while (my $b = <CONTIGS>) {
		chomp $b;
		next unless $b;
		my ($name, @sequence) = split (/\n/, $b);
		my $seq = join ("", @sequence);
		$sequences{$name} = uc $seq;
	}
	close CONTIGS;

	while (my($n, $s)=each(%sequences)){
		chomp($s);
		print "F:$fName\tSN: $n\n" unless (length($s)>0);
	}
	$/="\n";	
	return ();
}

my $listOfFiles= $ARGV[0];
open (LOF, "$listOfFiles") || die "ERROR: $ARGV[0]\n $!\n";
print "Summary for incomplete Genomes:\n";
while (my $file=<LOF>){
	checkForCompleteness($file);	
}
print "All Done!!\n";
