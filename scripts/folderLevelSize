#!/usr/bin/perl

# This script will only look at folders and sub-folders of the present working directory. You'll HAVE TO paste this script to the other folder if you want it's stats. Also make sure you have read permissions for the folders and sub folders before you run this script.

use strict;

#my $path= `pwd`;

my $level=$ARGV[0];
my $tmp=$$.".tmp";

`du -h > $tmp`;


if(! $ARGV[1]){ $level =1;}
$level++;

open (IN, $tmp);
while (my $line=<IN>){
	next if $line=~ /^#/;
	$line=~ s/\r//;
	chomp $line;
	next unless $line;

	my($size, $path)=split(/\t/, $line);
	my @levels=split(/\//, $path);
	
	print $size."\t".$path."\n" if (scalar(@levels) == $level);
}
unlink $tmp;
close IN;
