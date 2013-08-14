#! /usr/bin/perl
=head1 Description

	meanCov - Calculates the mean coverage for a list of 'features' from a Sorted Bam file.

=head2 Usage

	1) 	load the samtools module.
	2)	perl meanCov.pl -ref reference.fasta -l listOfFeatures.txt -b sortedBamFile.bam

=head3 Optional

	-ref	Reference Fasta (From NCBI)
	-l	List of all the features. (col1 = Accession number of contig; col2 = start position; col3 = stop position).
	-b	Sorted bam file

	IMPORTANT NOTE: Currently the script assumes that the refernce file has NCBI like headers. If you know perl feel free to tinker with your own copy of the script. Else contact Sunit before you use it.

=head2 For Suggestions/Comments/Feedback/Beer Contact:

	Sunit Jain, Aug 2011
	sunitj [AT] umich [DOT] edu

=cut

use strict;
use List::Util qw (sum);
use Getopt::Long;

my ($refFasta,$list,$sortedBam);
GetOptions(
	'ref=s'=>\$refFasta,
	'l|list=s'=>\$list,
	'b|bam=s'=>\$sortedBam,
	'h'=>sub{system('perldoc', $0); exit;},
);

my @cov;

open (LIST , $list);
my %positions;
while(my $line=<LIST>){
	next if $line=~ /^#/;
	chomp $line;
	next unless $line;

	my ($name, $p1, $p2)=split(/\t/, $line);
	my($start, $stop)= sort{$a <=> $b} ($p1, $p2);
	my $pos=$start."-".$stop;
	push(@{$positions{$name}},$pos);
}

open(FASTA, $refFasta);
my %contigs;
$/=">";
while(my $line= <FASTA>){
	chomp $line;
	next unless $line;

	my($desc, @seq)=split(/\n/, $line);
#	my $desc=~ s/^\>//;

	my @parts=split(/\|/, $desc);
	if ($positions{$parts[3]}){
		my $longName= join("\\\|", (@parts[0..3]));
		$contigs{$parts[3]}= $longName."\\\|";
	}
}
$/="\n";

my @cov;

print "#Accession no.	Start-Stop	MeanCoverage\n";
foreach my $c(keys %contigs){
# samtools mpileup -C50 -f ob3b_ref.fasta -r gi\|296253650\|gb\|ADVE01000127.1\|:2539-3309 ob3b_trimAln.sorted.bam | cut -f 4
	foreach my $pos(@{$positions{$c}}){
		my $R=$contigs{$c}.":".$pos;
		@cov= `samtools mpileup -f ob3b_ref.fasta -r $R $sortedBam | cut -f 4`;
		chomp(@cov);
		if (scalar(@cov) > 0){
			my $sum= sum(0, @cov);
			my $coverage=$sum/scalar(@cov);
			print $c."\t".$pos."\t".$coverage."\n";
		}
		else{
			print $c."\t".$pos."\t0\n";
		}
	}
}
