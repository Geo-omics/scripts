#!/usr/bin/perl

=head1 DESCRIPTION

	Use this script to append the name of the file to all the sequences in your fasta/fastq format files. Output will be a new file of the same name, created in a directory called "renamed" in the current working directory. Remember this was written to function as a batch job, so it doesn't require a file name just the folder and file extension.

=head1 USAGE

	perl addFileName2header.pl -path Folder_Path -ext extension_of_files

=head2 Options

	-prefix	prefix to the filename; this will appear at the end of the current header BEFORE the filename.
	-suffix	suffix to the filename; this will appear at the end of the current header AFTER the filename.
	-after	add the filename AFTER the existing header. (default: BEFORE the existing header.)
	-sep	separator between existing header and filename. (default: "_")
	-replace	replace something in the file name.
	-p	path to folder; use "." (dot, without the quotes) for current folder.
	-e	file extension to look for in folder; default= fasta
	-h	this page.

=head1 Suggestions/Corrections/Feedback/Beer

	Sunit Jain, sunitj@umich.edu

=cut

use strict;
use Getopt::Long;
use File::Basename;

my $help;
my $path; # Folder path
my $ext="fasta";
my $description=$$.".desc";
my $prefix="";
my $suffix="";
my $after;
my $sep="_";
my $replace=" ";
my $version=fileparse($0)."\tv2.1.3";
GetOptions(
	'prefix:s'=>\$prefix,
	'suffix:s'=>\$suffix,
	'after'=>\$after,
	'sep:s'=>\$sep,
	'replace:s'=>\$replace,
	'p|path:s'=>\$path,
	'e|ext:s'=>\$ext,
	'h|help'=>sub{system('perldoc', $0); exit;},
	'v|version'=>sub{print $version."\n"; exit;},
);
print "\# $version\n";
die "[ERROR: $0] Folder Path Required! See $0 -h for help on the usage" if !$path;

unless (-e "renamed"){mkdir("renamed", 0755)};

my @files=<$path/*.$ext>;

foreach my $fileName(@files){
	my $countSeqs=$ext eq "fastq" ? &parseFastq($fileName) : &parseFasta($fileName);
	print $fileName."\t".$countSeqs."\n";
}

sub parseFasta{
	my $fileName=shift;
	my $f = fileparse($fileName, ".fasta");
	my $fOut="./renamed/".$f.".fasta";
	$f=~ s/$replace//oe;
	my $addition=$prefix.$f.$suffix;
	open(IN, $fileName) || die $!;
	open(FASTA, ">".$fOut) || die $!;
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

		my $nuHead=$after ? $header.$sep.$addition : $addition.$sep.$header;
		$countSeqs++;
		print FASTA ">".$nuHead."\n".$seq."\n";
	}
	$/="\n";
	close(IN);
	close(FASTA);
	return $countSeqs;
}

sub parseFastq{
	my $fileName=shift;
	my $f = fileparse($fileName, ".fastq");
	my $fOut="./renamed/".$f.".fastq";
	$f=~ s/$replace//oe;
	my $addition=$prefix.$f.$suffix;
	my $countSeqs=0;
	open(FILE, $fileName) || die $!;
	open(FASTQ, ">".$fOut) || die $!;
	while (my $line=<FILE>){
		$line=&trim($line);
		if ($line=~ /^@/){
			$line=~s/\@//;
			
			my $nuHead=$after ? $line.$sep.$addition : $addition.$sep.$line;
			print FASTQ "\@".$nuHead."\n"; # Sequence Header
			$countSeqs++;
		
			$line=<FILE>; # Sequence
			$line=&trim($line);
			print FASTQ $line."\n";

			$line=<FILE>; # Quality Header
			print FASTQ "+\n";

			$line=<FILE>; # Quality
			$line=&trim($line);
			print FASTQ $line."\n";
		}
		else{ die "ERROR: Script Borked! Get Sunit (sunitj [ AT ] umich [ DOT ] edu)\n"; }
	}
	close FILE;
	close FASTQ;
	return $countSeqs;
}

sub trim{
	my $line=shift;
	chomp($line);
	$line=~ s/\r//;
	$line=~ s/^\s+//;
	$line=~ s/\s+$//;
	return $line;
}


