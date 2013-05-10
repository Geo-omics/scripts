#!/usr/bin/perl

=head1 DESCRIPTION

	Use this script to append the name of the file to all the sequences in your fasta/fastq format files. Output will be a new file of the same name, created in a directory called "renamed" in the current working directory. Remember this was written to function as a batch job, so it doesn't require a file name just the folder and file extension.

=head1 USAGE

	perl addFileName2header.pl -path Folder_Path -ext extension_of_files

=head2 Options

	-p	path to folder; use "." (dot, without the quotes) for current folder.
	-e	file extension to look for in folder; default= fasta
	-h	this page.

=head1 Suggestions/Corrections/Feedback/Beer

	Sunit Jain, sunitj@umich.edu

=cut

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;

my $help;
my $path; # Folder path
my $ext="fasta";
my $description=$$.".desc";

GetOptions(
	'p|path:s'=>\$path,
	'e|ext:s'=>\$ext,
#	'h|help'=>sub{system('perldoc', $0); exit;},
	'h|help'=>\$help,
);

pod2usage(1) if $help;

die "[ERROR: $0] Folder Path Required! See $0 -h for help on the usage" if !$path;

unless (-e "renamed"){mkdir("renamed", 0755)};

my @files=<$path/*.$ext>;

$|++;

foreach my $fileName(@files){
	my $countSeqs=$ext eq "fastq" ? &parseFastq($fileName) : &parseFasta($fileName);
	print $fileName."\t".$countSeqs."\n";
}

sub parseFasta{
	my $fileName=shift;
	my $f = basename($fileName, ".fasta");
	my $fOut="./renamed/".$f.".fasta";

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

		my $nuHead=$f."_".$header;
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
	my $f = basename($fileName, ".fastq");
	my $fOut="./renamed/".$f.".fastq";
	my $countSeqs=0;
	open(FILE, $fileName) || die $!;
	open(FASTQ, ">".$fOut) || die $!;
	while (my $line=<FILE>){
		$line=&trim($line);
		if ($line=~ /^@/){
			$line=~s/\@//;
			print FASTQ "\@".$f."_".$line."\n"; # Sequence Header
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


