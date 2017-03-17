#!/usr/bin/perl

=head1 DESCRIPTION

	countInstances.pl: count the number of times a value is seen in a column

=head1 USAGE

	perl countInstances.pl -in tab-delimited-input.txt -out output-filename.txt -col column-number

=head2 Options

	-in	[character]	input file name
	-out	[character]	output file name
	-col	[integer] column to subtotal; start counting from 1.
	-sum	[integer] sum the column # instead of incrementing by 1.

=head1 Author

	Sunit Jain, (Thu Jul 18 20:37:35 EDT 2013)
	sunitj [AT] umich [DOT] edu

=cut

use strict;
use Getopt::Long;

my $version="0.0.1b";
my $col=2;
my ($in, $out, $sum);
GetOptions(
	'in=s'=>\$in,
	'o|out=s'=>\$out,
	'col:i'=>\$col,
	'sum:i'=>\$sum,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system('perldoc', $0); exit;},
);

$col--;
#$sum--;

my %counts;
open(IN, "<".$in)|| die $!;
while(my $line=<IN>){
	chomp $line;
	next if $line=~ /^#/;
	
	my @cols=split(/\t/, $line);
	if($sum){
		$sum--;
		$counts{$cols[$col]}+=$cols[$sum];
	}
	else{
		$counts{$cols[$col]}++;
	}
}
close IN;

open(OUT, ">".$out)||die $!;
foreach my $keys(keys %counts){
	print OUT $keys."\t".$counts{$keys}."\n";
}
close OUT;
