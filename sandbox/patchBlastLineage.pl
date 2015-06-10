#!/usr/bin/perl

use strict;

# $ARGV[0]; Blast output.
# $ARGV[1]; name of hit (y/n).

my $fName= $ARGV[0];
my $lFile= "l_".$fName.".txt";

my %index;
open (LF, $lFile) || die "[err] $lFile not found\n".$!."\n";
while (my $desc=<LF>){
	my($gi, $taxa, $rank)=split(/\t/, $desc);
	chomp($gi);
	chomp($taxa);
	$index{$gi}=$taxa;
}
close LF;

open (OUT, ">taxaBlast_".$ARGV[0]);
open (BO, $fName) || die "[err] $fName not found\n".$!."\n";
while(my $line=<BO>){
	next if ($line=~ m/^\#/);
	my @blast=split(/\t/, $line);
	chomp(@blast);
	my($giTag, $gi, $id, $name)=split(/\|/, $blast[1]);
	chomp($gi);
	$blast[1]=$index{$gi}."\|".$id;
	$blast[1].="\|".$name if (lc($ARGV[1]) eq 'y');
	my $bo=join("\t", @blast);
	print OUT $bo."\n";
}
close BO;
close OUT;
