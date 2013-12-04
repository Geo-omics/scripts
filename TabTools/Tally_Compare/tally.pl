#!/user/bin/perl
use strict;
use Getopt::Long;

my $in;
my $master;
my $out=$$.".out";
my $bs=0;
my $printValues;
GetOptions(
	'i:s'=>\$in,
	'm:s'=>\$master,
	'o:s'=>\$out,
	's:f'=>\$bs,
	'values|value'=>\$printValues,
);

open(IN, $in)|| die "[Error] $in: $!\n";
my %seen;
while (my $line=<IN>){
	next if ($line=~ m/^#/);
	chomp $line;
	next unless $line;
	$line=~ s/\r//g;

	my @lineParts=split(/\t/, $line);
	if(! $printValues){
		$seen{$lineParts[0]}++ if ($lineParts[-1] >= $bs);
	}
	else{
		$seen{$lineParts[0]}=$lineParts[-1];
	}
}
close IN;

open(MASTER, $master) || die "[Error] $master: $!\n";
open(OUT, ">".$out);
while (my $line=<MASTER>){
	next if ($line=~ m/^#/);
	chomp $line;
	next unless $line;
	$line=~ s/\r//g;

	if ($seen{$line}){
		print OUT $line."\t".$seen{$line}."\n";
	}
	else{
		print OUT $line."\t0\n";
	}
}
close MASTER;
close OUT;
exit;
