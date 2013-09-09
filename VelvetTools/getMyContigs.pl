#! /usr/bin/perl

# USAGE: perl getMyContigs.pl <CoverageOutputFromReadCoverage Script> <Your Contig Shortlist file> <OUTPUT.list>

use strict;

my $readCov=$ARGV[0];
my $list=$ARGV[1];
my $OUT=$ARGV[2];

die "Incorrect number of files input\nUSAGE: perl getMyContigs.pl <CoverageOutputFromReadCoverage Script> <Your Contig Shortlist file> <OUTPUT.list>" if (scalar(@ARGV) != 3);

open(LIST, $list)|| $!;
my %LIST;
while(my $line=<LIST>){
	chomp;
	next unless $line;
	next if $line=~ /^#/;
#	NODE_14_length_2679_cov_8.406121
	my @headerParts=split(/\_/, $line);
	$LIST{$headerParts[1]}++;
}
close LIST;

my %READS;
open(READ, $readCov)|| $!;
while(my $line=<READ>){
	next if $line=~ /^#/;
	chomp $line;
	next unless $line;

	my ($contigName, $size, @reads)=split(/\t/, $line);

	next unless $LIST{$contigName};
	foreach my $r(@reads){
		$READS{$r}++;
	}
}
close READ;
undef %LIST;

print "Total # Reads Mapped to this bin:".keys(%READS)."\n";
open(OUT, ">".$OUT)|| die $!;
foreach my $r(keys %READS){
	print OUT $r."\n";
}
close OUT;
