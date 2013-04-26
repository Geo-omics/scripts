#!/usr/bin/perl
=head1 Name

	getRandomData - version 0.0.2

=head2 Description

	Gets a 'random' subset of the data provided.

=head2 Usage

	perl getRandomLines.pl -i <File Name>

=head3 Optional

	-o: output file name; default: random<processID>.out
	-p: What % of the data set do you want as your subset?; default: 10
	-d: dataset type; default: "tab"

=head3 File Types

	This script maybe used for the following file types:

	"tab" tab-delimited data.
	"space" space-delimited data.
	"fastq"	fastq file w/ seq and quality scores.
	"qual"	fastq file with just the quality scores and each new header begins with a '+'
	"fasta" fasta files.
	"newline" each entry in your dataset is in a 'newline'.
	"comma" comma-seperated file
	
	Use the above mentioned key-words in the "-d" flag to describe the filetype to the script.

=head2 Author

	Sunit Jain, July 2011
	sunitj [AT] umich [DOT] edu

=cut

use strict;
use Getopt::Long;


my $in;
my $out="random".$$.".out";
my $delimiter="tab";
my $per= 10;
my $version="0.0.2";
GetOptions(
	'i:s'=> \$in,
	'p:f'=>\$per,
	'o:s'=>\$out,
	'd:s'=>\$delimiter,
	'v'=> sub{ print $0.":\tVersion ".$version."\n"; exit; },
	'h'=> sub{ system("perldoc", $0); exit; }
);

my($delim, $fileSize)=&delimiter;
print "Dataset Size: ". $fileSize."\n";

my $totalLines= $fileSize * ($per/100);
my %rNum;
while (scalar(keys %rNum) < $totalLines){
	my $r=int(rand($totalLines));
	$rNum{$r}++;	
}
print "Sample Size: ".int($totalLines)."\n";

open (IN, $in) || die "[ERROR] $in $!\n";
open (OUT, ">".$out);

$/=$delim;
my $lNum=0;
while(my $line=<IN>){
	next if ($line=~ m/^#/);
	chomp $line;
	next unless $line;
	$lNum++;
	if ($rNum{$lNum}){
		print OUT $delim unless ($delim=~ m/\n/o);
		print OUT $line."\n";
	}
}
close IN;
close OUT;

sub delimiter{
	if ((lc($delimiter) eq "tab")|| (lc($delimiter) eq "newline") || (lc($delimiter) eq "comma") || (lc($delimiter) eq "space")){
		my $fs=`wc -l $in`;
		return ("\n", $fs);
	}
	elsif (lc($delimiter) eq "fastq"){
		my $fs = `grep -c \'\@\' $in`;
		return ("\@", $fs);
	}
	elsif (lc($delimiter) eq "qual"){
		my $fs = `grep -c \'\+\' $in`;
		return ("\+", $fs);
	}
	elsif (lc($delimiter) eq "fasta"){
		my $fs = `grep -c \'\>\' $in`;
		return ("\>", $fs);
	}
	else{
		print "Invalid Delimiter! $delimiter \n. See \"".$0." -h \" for more information\n";
	}
}
