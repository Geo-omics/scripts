#! /usr/bin/perl

use strict;

my $fasta=$ARGV[0];

## Read Fasta File and compute length ###
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
}
$/="\n";
close(FASTA);

my @sort = sort {$b <=> $a} @allLen;
my $n50; 
my $L50;
foreach my $val(@sort){
	$n50+=$val;
	$L50++;
	if($n50 >= $totalLength/2){
		print "N50:\t$val\n";
		print "L50:\t$L50\n";
		last; 
	}
}
