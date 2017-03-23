#!/usr/bin/perl

use strict;

my $in= $ARGV[0];
my $out= $ARGV[1];

open( SUMM, $in)|| die "$!\n";
open (OUT, ">".$out);

while(my $line=<SUMM>){
	my ($id);
	if($line=~ m/<ID>(\d*)<\/ID>/i){
		print OUT $1."\t";	
	}
	elsif($line=~ m/<Item Name\=\"Title\" Type\=\"String\">([\w\W]*)<\/Item>/i){
		print OUT $1;
	}
	elsif($line=~ m/<\/DocSum>/i){
		print OUT "\n";
	}
	else{
		next;
	}
}
