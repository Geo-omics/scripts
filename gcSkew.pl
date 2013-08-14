#!/usr/bin/perl
=head1 DESCRIPTION

	gcSkew.pl v0.0.2	-	takes a fasta formatted nucleotide file and calculates GC skew [(C-G)/(C+G)]*100 content for the whole sequence in sliding windows. You can specify the window length, and amount of overlap between successive windows.

=head1 USAGE

	perl gcSkew.pl -f file.fasta

=head2 OPTIONS

	-f or -fasta:	Fasta file	(REQUIRED)
	-o or -out:	Output File	(DEFAULT: process_id.gcskew)
	-w or -window:	Window Size (DEFAULT: 4)
	-overlap:	Overlap between two windows (DEFAULT: 0)

	-h or -help:	This Documentation
	-v or -version:	Script version

=head1 AUTHOR

	Sunit Jain, May 2013
	sunitj [AT] umich [DOT] edu

=cut

use strict;
use Getopt::Long;

my $help;
my $version=0.0.3;
my $fasta;
my $out=$$.".gcskew";
my $window=4;
my $overlap=0;

GetOptions(
	'f|fasta=s'=>\$fasta,
	'w|window:i'=>\$window,
	'o|out:s'=>\$out,
	'overlap:i'=>\$overlap,
	'h|help'=>\sub{system('perldoc', $0); exit;},
);

die "Overlap has to be less than the Window size\n" if ($overlap >= $window);
my $step = $window - $overlap;

print "Input File:\t".$fasta."\n";
print "Output File:\t".$out."\n";
print "Window:\t".$window."\n";
print "Step:\t".$step."\n";
print "Overlap:\t".$overlap."\n";

open(FASTA, $fasta)|| die $!;
open(OUT, ">".$out)|| die $!;
print OUT "# Header\tStartPos\tStopPos\tGC-Skew(\%)\n";
$/=">";
my $totalWindows;
my $totalSequences;
while(my $line=<FASTA>){
	chomp $line;
	next unless $line;

	my($header,@seq)=split(/\n/, $line);
	my $sequence=join("", @seq);

	$totalWindows+= gcSkew($header, $sequence);
	$totalSequences++;
}
$/="\n";
close FASTA;
close OUT;

print "Total Windows:\t".$totalWindows."\n";
print "Total Sequences:\t".$totalSequences."\n";


exit;

sub gcSkew{
	my ($header, $seq)=@_;
	my $seqLen=length($seq);
	my $windows=0;
	foreach(my $start=0; $start<=$seqLen; $start+=$step){
		my ($A, $T, $G, $C)=(0,0,0,0);
		my $windowSeq=substr($seq, $start, $window);
		if (length($windowSeq) == $window){
		    while ( $windowSeq =~ /G/ig ) { $G++ }
		    while ( $windowSeq =~ /C/ig ) { $C++ }

			my $skew=0;
			if (($C + $G) > 0){
				$skew= (($C - $G)/($C + $G))* 100;
			}
			my $GCSkew = sprintf( "%.4f", $skew );
			print OUT $header."\t".$start."\t".($start+$step)."\t".$GCSkew."\n";
		}
		else{
			print OUT $header."\t".$start."\t".$seqLen."\t0.0000\n";
		}
		$windows++;
	}
	return $windows;
}
