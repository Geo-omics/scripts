#!/usr/bin/perl

=head1 USAGE

	perl sangerSeqParser.pl -p Folder_Path

=head2 Options

	-o	output file name; default= processID.fasta
	-e	file extension to look for in folder; default= fasta
	-h	this page.

=head1 Suggestions/Corrections/Feedback/Beer

	Sunit Jain, sunitj@umich.edu

=cut

use strict;
use Getopt::Long;
use File::Basename;

my $path; # Folder path
my $ext="fasta";
my $out;

GetOptions(
	'p:s'=>\$path,
	'o:s'=>\$out,
	'e:s'=>\$ext,
	'h|help'=>sub{system('perldoc', $0); exit;},
);

$path= `pwd` if !$path;
chomp $path;
$out=$$.".".$ext if !$out;

my @files=<$path/*.$ext>;

open(FASTA, ">".$out) || die $!;

foreach my $f(@files){
	my $fhIN;
	
	open($fhIN, $f) || die $!;
	my @sequence;
	while(my $line=<$fhIN>){
		chomp $line;
		$line=~ s/\r//;
		next unless $line;
		push(@sequence, $line);
	}
	close($fhIN);
	my $seq= join("", @sequence);
	
	my $nuHead=fileparse($f);

	print FASTA ">".$nuHead."\n".$seq."\n";
}
close (FASTA);
