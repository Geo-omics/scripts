#!user/bin/perl

=head1 Description

	parseFastq - parses illumina Fastq files and splits them into fasta and qual files.

=head2 Usage

	perl parseFastq.pl [-i Illumina FastQ file]

=head3 Optional

	-f	Fasta Output File.
	-q	Qual Output File.
	-p	Phred offset; default 33;
	-summary	provides a sequence length distribution for the fastq file.

=head2 For Suggestions/Comments/Feedback/Beer Contact:

	Sunit Jain, Aug 2011
	sunitj [AT] umich [DOT] edu

=cut

use strict;
use Getopt::Long;
use List::Util qw(sum);

my $in;
my $fasta; #output
my $qual; #output
my $offset=33;
my $summary;
my $qLim=1;
my $version="parseFastq.pl v0.6.1";

GetOptions(
	'i|in:s'=>\$in,
	'f|fasta:s'=>\$fasta,
	'q|qual:s'=>\$qual,
	'p|offset:i'=>\$offset,
	'summary'=>\$summary,
	'min_qual:i'=>\$qLim,
	'v|version'=>\sub{print $version."\n"; exit;},
	'h|help'=>sub{system('perldoc', $0); exit;},
);

if (! $offset){$offset=33;}
my (@fileName)=split(/\./, $in);
pop(@fileName);
my $fName=join("\.", @fileName);

if (! $in){
	sub{system('perldoc', $0); exit;}
}
if (! $fasta){
	$fasta=$fName.".fasta";
	print "#Fasta File: $fasta\n";
}

open (ILL, $in)||die "[ERROR: $0] $in : $! \n";
open (FASTA, ">".$fasta)|| die "[ERROR: $0] $fasta : $! \n";
open (QUAL, ">".$qual)|| die "[ERROR: $0] $qual : $! \n" if ($qual);

my $lNum=0;
my $seqName;
my $numSeqs=0;
my (%stats, %qualSummary);
while(my $line=<ILL>){
	$line= &trim($line);
	next unless $line;

	if ($line=~ /^@(\w+)/){
		#Sequence Header
		my $header=$line;
		$header=~ s/^@/>/;
		$numSeqs++;

		#Sequence
		$line=<ILL>;
		my $seq= &trim($line);
		$stats{length($seq)}++ if ($summary);

		#Quality Header
		$line=<ILL>;

		#Quality Scores
		$line=<ILL>;
		$line= &trim($line);
		my @score=split (//, $line);
		
		# Error in script if length of the sequence != length of the quality score.
		die "[ERROR: line $.] $header \nScript Borked! Unexpected Sequence Format.\nGet Sunit (sunitj [ AT ] umich [ DOT ] edu)\n" if (length($seq)!=length($line));

		# score conversion.
		@score=map{(ord($_)-$offset)} @score;
		my $scores=join(",", @score);

		my $avgScore=sum(@score)/scalar(@score);
		$header.="_avgScore_".$avgScore;

		if(($avgScore > 41) || ($avgScore < 0)){
			print STDERR "[ERROR: $0] $header has an Average Quality Score of $avgScore\n";
			die "Incorrect offset: $offset\nChange to the correct offset by using the \"-offset\" flag\n";
		}
		unless($summary){		
			if(($avgScore > 30) && ($avgScore <=41)){
				$qualSummary{"Q1"}++;
			}
			elsif(($avgScore > 20) && ($avgScore <=30)){
				$qualSummary{"Q2"}++;
			}
			elsif(($avgScore > 10) && ($avgScore <=20)){
				$qualSummary{"Q3"}++;
			}
			elsif(($avgScore > 0) && ($avgScore <=10)){
				$qualSummary{"Q4"}++;
			}
			elsif($avgScore == 0){
				$qualSummary{"Zero"}++;
			}
		}
		print2Files($header, $seq, $scores) unless ($avgScore < $qLim);
	}
	else{ die "[ERROR: $0] Script Borked! Unexpected Sequence Format.\nGet Sunit (sunitj [ AT ] umich [ DOT ] edu)\n"; }
}
close FASTA;
close QUAL if ($qual);
close ILL;

if ($summary){
	my ($lenSum, $total);
	print "\>".$in."\n";
	while(my($len, $multiple)=each(%stats)){
		print "Length:".$len."\tTotal:".$multiple."\n";
		$lenSum+=($len*$multiple);
		$total+=$multiple;
	}
	my $avgLen=$lenSum/$total;

	print "Avg Len:\t".$avgLen."\n\n";

	my @bins=qw(Q1 Q2 Q3 Q4 Zero);
	foreach my $b(@bins){
		if (!$qualSummary{$b}){$qualSummary{$b}=0;}
		my $perc=sprintf("%.2f", (($qualSummary{$b}/$numSeqs)*100));
		print $b."\t".$perc."\n";
	}
}


sub print2Files{
	my @data=@_;

	print FASTA $data[0]."\n";
	print FASTA $data[1]."\n";
	
	if ($qual){
		print QUAL $data[0]."\n";
		print QUAL $data[2]."\n";
	}
	return;
}

sub trim{
	my $line=shift;
	chomp $line;
	$line=~ s/\r//;
	$line=~ s/^\s+//;
	$line=~ s/\s+$//;
	return $line;
}

exit;
