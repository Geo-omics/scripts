#! /usr/bin/perl

use strict;

my $fasta=$ARGV[0];

## Read Fasta File and compute N50, L50, N95 and L95 ##
my $length;
my $totalLength;
my $totalContigs;
my @allLen;
open(FASTA, $fasta)|| die $!;
$/=">";
while(my $line=<FASTA>){
	chomp $line;
	next unless $line;

	my ($header, @sequence)=split(/\n/, $line);
	my $length=length(join("", @sequence));

	push (@allLen, $length);
	$totalLength += $length;
	$totalContigs++;
}
$/="\n";
close(FASTA);

my @sortedLen = sort {$b <=> $a} @allLen;
my $cumLen;
my $numContig;
print "Total_Contigs:\t$totalContigs\n";
foreach my $len(@sortedLen){
	$cumLen+=$len;
	$numContig++;
	if ($cumLen >= $totalLength * 0.95) {
		print "N95:\t$len\n";
		print "L95:\t$numContig\n";
	}
	if($cumLen >= $totalLength * 0.50){
		print "N50:\t$len\n";
		print "L50:\t$numContig\n";
		last;
	}
}
