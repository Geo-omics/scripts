#! /usr/bin/perl

use strict;
use Getopt::Long;

=head2 NAME
	
	U2T : Converts  U -> T and removes gaps.
	
=head2 Description

	 perl U2T.pl -in <input fasta file> -out <output fasta file>

=cut

my $seqFile;
my $out;

GetOptions(
	"in:s" => \$seqFile,
	"out:s"	=>	\$out,
	"h|help"	=>	sub {system('perldoc', $0); exit;},
);

$/=">";
open (SEQ, $seqFile) || die "Couldn't open $seqFile\n";
open (OUT, ">".$out);
while (my $line = <SEQ>) {
	next if $line=~ m/^#/;
    chomp $line;
	$line=~ s/\r//;
    next unless $line;

	my($seqDesc, @sequence)=split(/\n/, $line);
	my $seq=join("", @sequence);
	
	$seq=~ tr/ACGTU/ACGTT/;
	$seq=~ s/[\.\-\s]//g;
	print OUT ">". $seqDesc."\n".$seq."\n";
}
close SEQ;
close OUT;
