#!/usr/bin/perl

use warnings;
use strict;
use File::Basename;

my $inFasta = $ARGV[0];
my $baseName = basename($inFasta, qw/.fasta .fna/);
my $inQual = $baseName . ".qual";
my $outFastq = $baseName . ".fastq";

my %seqs;

$/ = ">";

open (FASTA, "<$inFasta");
my $junk = (<FASTA>);

while (my $frecord = <FASTA>) {
	chomp $frecord;
	my ($fdef, @seqLines) = split /\n/, $frecord;
	my $seq = join '', @seqLines;
	$seqs{$fdef} = $seq;
}

close FASTA;

open (QUAL, "<$inQual");
$junk = <QUAL>;
open (FASTQ, ">$outFastq");

while (my $qrecord = <QUAL>) {
	chomp $qrecord;
	my ($qdef, @qualLines) = split /\n/, $qrecord;
	my $qualString = join ' ', @qualLines;
	my @quals = split / /, $qualString;
	print FASTQ "@","$qdef\n";
	print FASTQ "$seqs{$qdef}\n";
	print FASTQ "+\n";
	foreach my $qual (@quals) {
		print FASTQ chr($qual + 33);
	}
	print FASTQ "\n";
}

close QUAL;
close FASTQ;
