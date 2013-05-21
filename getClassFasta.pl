#! /usr/bin/perl

=head1 Description

	This program takes a class file, a names file, a fasta file and extracts the seqs for the contigs presents in the desired class.

=head2 Usage

	perl getClassFasta.pl -cls <CLASS File> -names <NAMES File> -fasta <Split FASTA File> -num <CLASS NUMBER>
	
=head3 Options

	-cls	:	the .cls file
	-names	:	the .names file
	-fasta	:	The Concatenated fasta file that was used with the tetramer script.
	-num	:	the class number you're interested in.
    
=head1 Questions/Comments/Suggestions/Accolades/Beer

	Sunit Jain, sunitj [AT] umich [DOT] edu

=cut


use strict;
use Getopt::Long;
use Pod::Usage;

my ($help, $class, $names, $fasta, $classNum);
my $version="0.0.02";

GetOptions(
	"cls=s"=>\$class,
	"names=s"=>\$names,
	"fasta=s"=>\$fasta,
	"num=i"=>\$classNum,
	'h|help'=>\$help,
);

pod2usage(1) if $help;

print "getClassFasta v$version\n";

die if (! $class or ! $names or ! $fasta or ! $classNum);

# Parse *.cls file to get SeqID for all Seqs in the desired class
my %clsHash;
open ( CLS, $class) || die "ERROR: $class.\n".$!;
	while (my $line=<CLS>){
		chomp($line);
		unless ($line=~ /^%/){
			my ($seqNum, $cls)=split(/\t/,$line);
			if ($cls==$classNum){
				$clsHash{$seqNum}=$cls;	# %clsHash {Sequence Number  => Class Number}
			}
		}
	}
close CLS;

# Parse the *.names file to id the seq names in the fasta file using the classHash from above.
my %seqNames;
open (NAMES, $names) || die "ERROR: $names.\n".$!;
	while (my $line=<NAMES>){
		chomp($line);
		unless ($line =~ /^%/){
			my ($seqNum, $seqName, $seqFastaName)=split(/\t/, $line);
			if ($clsHash{$seqNum}){
				$seqNames{$seqFastaName}=$seqName; # %seqNames {Name of Seq part of the contig => Name of the whole contig}
			}
		}
	}
close NAMES;
undef %clsHash;

# Parse the fasta file to get the desired sequences, using the seqNames hash from above.
my $outFile=$classNum.".fasta";
open (OUT,">".$outFile )|| die "ERROR: $outFile.\n".$!;
open (FASTA, $fasta) || die "ERROR: $fasta.\n".$!;
	$/= ">";
	while (my $line = <FASTA>) {
		chomp $line;
		next unless $line;
		my ($name, @sequence) = split (/\n/, $line);
		my $seq = uc(join ("", @sequence));
		if ($seqNames{$name}){
			print OUT ">$name\n$seq\n";
		}
	}
	$/= "\n";
close FASTA;
close OUT;

