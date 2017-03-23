#!/usr/bin/perl

use strict;

my $in=$ARGV[0];
my $out=$ARGV[1];

open(IN, $in)|| die $!."\n";
open(OUT, ">".$out);

while(my $line=<IN>){
	next if $line=~ /^#/;
	chomp $line;
	$line=~ s/\r//;
	next unless $line;

	print OUT $line."\n";
}
close IN;
close OUT;
