#! /usr/bin/perl

=head1 DESCRIPTION

	Use this script to create an annotation file and a concatenated fasta for ESOM binning.

=head1 USAGE

	perl makeAnnotationFile -path Folder_Path -ext extension_of_files

=head2 Options

	-p or path	path to folder; use "." (dot, without the quotes) for current folder.
	-e or ext	file extension to look for in folder; default= fasta
	-a or ann	name of the output annotation file; default= esom.ann
	-c or cat	name of the output concatenated fasta file; default= esom.fasta
	-O or dir	name of the output directory; default= ESOM
	-h	this page.

=head1 Suggestions/Corrections/Feedback/Beer

	Sunit Jain, sunitj@umich.edu
	January 2013

=cut

use strict;
use Getopt::Long;
use File::Basename;
use File::Spec;

my $path; # Folder path
my $ext="fasta";
my $annotationFile="esom.ann";
my $concatenatedFasta="esom.fasta";
my $outDir="ESOM";

GetOptions(
	'p|path:s'=>\$path,
	'e|ext:s'=>\$ext,
	'a|ann:s'=>\$annotationFile,
	'c|cat:s'=>\$concatenatedFasta,
	'O|dir:s'=>\$outDir,
	'h|help'=>sub{system('perldoc', $0); exit;},
);

die "[ERROR: $0] Folder Path Required! See $0 -h for help on the usage" if !$path;


unless (-e $outDir){mkdir($outDir, 0755)};
my $ann=File::Spec->catfile( $outDir, $annotationFile);
my $catFasta=File::Spec->catfile( $outDir, $concatenatedFasta);


my @files=<$path/*.$ext>;

$|++;

open(FASTA, ">".$catFasta) || die $!;
open(ANN, ">".$ann) || die $!;

my $class=0;
foreach my $fileName(@files){
	my $countSeqs=	&parseFasta($fileName);
	$class++;
	print $fileName."\t".$countSeqs."\n";
}
close(IN);
close(FASTA);


sub parseFasta{
	my $fileName=shift;

	open(IN, $fileName) || die $!;

	my ($prevHeader, $flag);
	$/=">";
	my $countSeqs=0;
	while(my $line=<IN>){
		chomp $line;
		$line=~ s/\r//;
		next unless $line;

		my($header, @sequence)=split(/\n/, $line);
		my $seq= join("", @sequence);
		if (length($seq)==0){
			$prevHeader=$header;
			$flag=1;
			next;
		}
		elsif($flag==1){
			$header=$prevHeader."_".$header;
			$flag==0;
			$prevHeader="";
		}

# Beautify
		$header=~ s/\s+/\_/g;
		$header=~ s/\W+/\_/g;
		$header=~ s/\_+/\_/g;
		$header=~ s/\_+$//;
		$header=~ s/^\_+//;

		$countSeqs++;
		print FASTA ">".$header."\n".$seq."\n";
		print ANN $header."\t".$header."\t".$class."\n";
	}
	$/="\n";
	return $countSeqs;
}

