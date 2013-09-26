#!/usr/bin/perl

=head1 Description

	renameHeaders.pl : Changes the names of the sequence headers by inserting a user defined prefix.

=head2 Usage

	perl renameHeaders.pl -f fasta.file -prefix "seq header prefix" -d " "
	OR
	perl renameHeaders.pl -f fasta.file -suffix "seq header prefix" -d " "

=head3 Options
	
	-suffix	add after the header
	-prefix add before the header
	-d delimiter to follow the prefix; default= "_"
	-o output file; default= "processID.fasta"

=head1 Comments/Feedback/Beer

	Sunit Jain, sunitj [AT] umich [DOT] edu
	May 2012

=cut

use strict;
use Getopt::Long;

my $fasta;
my $prefix;
my $suffix;
my $out=$$.".fasta";
my $delim="_";

GetOptions(
	'f:s'=>\$fasta,
	'prefix:s'=>\$prefix,
	'suffix:s'=>\$suffix,
	'o:s'=>\$out,
	'd:s'=>\$delim,
	'h'=>sub{system('perldoc', $0); exit;},
);

my $insert = ($prefix ? $prefix : $suffix);
if (! $fasta || ! $insert){	system('perldoc', $0); exit; }

open (IN, $fasta)||die "[ERROR] $fasta: $!\n";
open (OUT, ">".$out);
my %sampleName;
my $h=0;
while(my $line=<IN>){
	next if ($line=~ m/^#/);
	chomp $line;
	$line=~ s/\r//;
	next unless $line;

	if ($line=~ m/^>/){
		$line=~ s/^>//;
		$sampleName{$line}++;
		my $newName;
		if ($prefix){
			$newName= $insert.$delim.$line;
		}
		elsif ($suffix){
			$newName= $line.$delim.$insert;
		}
		print OUT ">".$newName."\n";
		$h++;
	}
	else{
		print OUT $line."\n";
	}
}
close IN;
close OUT;

#print "Name of the Sample was:\t";
foreach my $k(keys %sampleName){
	print $k."\t".$sampleName{$k}."\n" if ($sampleName{$k} > 1);
}
print "Total Sequences: ".$h."\n";
