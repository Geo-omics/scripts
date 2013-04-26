#! /usr/bin/perl
=head1 Description

	bamTools: Yup, until I think of a better name that's what i'm calling it!

	give me a sorted bam(-b) file, a list of contigs(-l). I'll tell you which reads mapped(-reads) onto those list of contigs and the mean coverage for each contig.

=head2 Preparation

	[1] Make sure you have an indexed sorted bam file in the same folder; and
	[2] The samtools module is loaded.

=head2 Usage

	perl bamTools.pl -f reference.fasta -b sorted.bam -l contigNames.list [-mean output.mean.txt -reads output.reads.txt]

=head2 Options

	-b or bam:	sorted bam file;
	-l or list:	list of contig used for the mapping.
	-mean :	Mean Coverage for each contig.
	-r or reads : Number of reads for each contig.
	-f or fasta : Reference fasta file.

=head2 Questions/Comments/Feedback/Beer

	October 2012, Sunit Jain (sunitj [AT] umich [DOT] edu)

=cut


use strict;
use List::Util qw (sum);
use Getopt::Long;

my ($refFasta,$list, $sortedBam, $mean, $reads);

GetOptions(
	"b|bam=s"=>\$sortedBam,
	"l|list=s"=>\$list,
	"f|fasta=s"=>\$refFasta,
	"mean:s"=>\$mean,
	"r|reads:s"=>\$reads,
	"h|help"=>sub{system('perldoc', $0); exit;},
);

if (!$reads && !$mean){
	$mean=$$.".mean.txt";
	$reads=$$.".numReads.txt";
}

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
while(my $line=<FASTA>){
	chomp $line;
	next unless $line;

	my($desc, @seq)=split(/\n/, $line);
#	my $desc=~ s/^\>//;
############################################################
	my @parts=split(/\|/, $desc);
	if ($positions{$parts[3]}){
		my $longName= join("\\\|", (@parts[0..3]));
		$contigs{$parts[3]}= $longName."\\\|";
	}
############################################################
}
$/="\n";

open(MEAN, ">".$mean) if ($mean);
open(READS, ">".$reads) if ($reads);
print MEAN "#Accession no.\tStart-Stop\tMeanCoverage\n";
print READS "#Accession no.\tStart-Stop\tNumReads\n";
foreach my $c(keys %contigs){
# samtools mpileup -C50 -f ob3b_ref.fasta -r gi\|296253650\|gb\|ADVE01000127.1\|:2539-3309 ob3b_trimAln.sorted.bam | cut -f 4
	foreach my $pos(@{$positions{$c}}){
		my $R=$contigs{$c}.":".$pos;
		if ($mean){
			my @cov= `samtools mpileup -f ob3b_ref.fasta -r $R $sortedBam | cut -f 4`;
			chomp(@cov);
			if (scalar(@cov) > 0){
				my $sum= sum(0, @cov);
				my $coverage=$sum/scalar(@cov);
				print MEAN $c."\t".$pos."\t".$coverage."\n";
			}
			else{
				print MEAN $c."\t".$pos."\t0\n";
			}
		}
		if ($reads){
			my @cov= `samtools view -F0x4 $sortedBam $R | cut -f 1 `;
			my $numReads=scalar(@cov);
			if ($numReads > 0){
				print READS $c."\t".$pos."\t".$numReads."\n";
			}
			else{
				print READS $c."\t".$pos."\t0\n";
			}
		}
	}
}

