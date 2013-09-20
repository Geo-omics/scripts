#! /usr/bin/perl
=head1 Description

	bamTools: Yup, until I think of a better name that's what i'm calling it!

	give me a sorted bam(-b) file, a list of contigs(-l). I'll tell you which reads mapped(-reads) onto those list of contigs and the mean coverage for each contig.

=head2 Dependencies

	Samtools 1.1.17 or above

=head2 Preparation

	[1] Make sure you have an indexed sorted bam file in the same folder; and
	[2] The samtools module is loaded.
	[3] The contigs list file may contain 2 additional columns in the form of start and stop position on the contig.

=head2 Usage

	perl bamTools.pl -f reference.fasta -b sorted.bam -l contigNames.list [-out processID.bamstats]

=head2 Options

	-b or bam:	sorted bam file;
	-l or list:	list of contig used for the mapping.
	-no_cov :	Don't calculate 'Mean Coverage' for each contig/region.
	-all	:	Prints results for everything, even the reference regions that didn't recruit any reads.
	-o or out : Number of reads for each contig/region.
	-f or fasta : Reference fasta file.

=head2 Questions/Comments/Feedback/Beer

	October 2012, Sunit Jain (sunitj [AT] umich [DOT] edu)

=cut


use strict;
use List::Util qw (sum);
use Getopt::Long;

my ($refFasta,$list, $sortedBam, $noCov, $all);
my $out=$$.".bamstats";

GetOptions(
	"b|bam=s"=>\$sortedBam,
	"l|list=s"=>\$list,
	"f|fasta=s"=>\$refFasta,
	"no_cov"=>\$noCov,
	"o|out:s"=>\$out,
	"all"=>\$all,
	"h|help"=>sub{system('perldoc', $0); exit;},
);

open (LIST , $list)|| die "ERROR: Something wrong with LIST file: $list\n: $!\n";
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

open(FASTA, $refFasta)|| die "ERROR: Something wrong with REF FASTA file: $refFasta\n: $!\n";
my %contigs;
$/=">";
while(my $line=<FASTA>){
	chomp $line;
	next unless $line;

	my($desc, @seq)=split(/\n/, $line);
	my @description=split(/\s+/, $desc);
	$desc=shift @description;
	$desc=~ s/^\>//;
	next unless $positions{$desc};
	
	my $seqLen=length(join("", @seq));
	$contigs{$desc}= $seqLen;
# This bit was modified for a different set of input files.
############################################################
#	my @parts=split(/\|/, $desc);
#	if ($positions{$parts[3]}){
#		my $longName= join("\\\|", (@parts[0..3]));
#		$contigs{$parts[3]}= $longName."\\\|";
#	}
############################################################
}
$/="\n";

my $log="log";
open(READS, ">".$out)|| die $!;
print READS "#Ref_Sequence\tStart-Stop\tNumMappedReads\t#Bases_Covered\tMeanCoverage\t#Reads_Normalized_by_Gene/Contig_Length\n";
foreach my $c(keys %contigs){
# samtools mpileup -C50 -f ob3b_ref.fasta -r gi\|296253650\|gb\|ADVE01000127.1\|:2539-3309 ob3b_trimAln.sorted.bam | cut -f 4
	foreach my $pos(@{$positions{$c}}){
		my $R=$c.($pos eq "-" ? "" : ":$pos");
		my ($coverage, $numReads, $numBasesCovered)=(0,0,0);
		unless($noCov){
#			print "samtools mpileup -d 8000 -f $refFasta -r $R $sortedBam  | cut -f 4\n";
			my @cov=`samtools mpileup -d 8000 -f $refFasta -r $R $sortedBam  | cut -f 4 2>$log`;
			chomp(@cov);
			$numBasesCovered=scalar(@cov);
			unless($all){ next if $numBasesCovered==0; }
			my $sum= sum(@cov);
			if($numBasesCovered>0){
				$coverage=$sum/$numBasesCovered;
			}
			else{
				$coverage=0;
			}
		}
		$numReads= `samtools view -c -F0x4 $sortedBam $R`;
		chomp $numReads;
		my $numReads_by_geneLen;
		if($pos eq "-"){
			$pos=$contigs{$c};
			$numReads_by_geneLen= $numReads/$contigs{$c}; # (number of reads / seqLength)
		}
		else{
			my ($a, $b)=split(/\-/, $pos);
			my $geneLen=$b-$a;
			my $numReads_by_geneLen=$numReads/$geneLen;
		}
		print READS $c."\t".$pos."\t".$numReads."\t".$numBasesCovered."\t".$coverage."\t".$numReads_by_geneLen."\n";
	}
}

