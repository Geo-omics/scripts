#!/usr/bin/perl

=head1 Description

Extract blast outputs based on any of the following: evalue, bitscore, percent ID, alignment length, % query coverage and & inferred global coverage (requires slightly modified blast output). Read `help` for more.

=head2 Acceptable BLAST Output Format

In order to filter your blast outputs by query coverage or inferred global coverage, you'll need to run your command line blast with:

	'-outfmt "6 std qcovs qlen slen"'

This will output your standard tabular blast output (6 std) and then will add the columns for query coverage (qcovs), query length (qlen) and subject length (slen).

Once your blast is finished, you may use the flags:

	"-cov 80 -covCol 13"
	to accept hits with query coverage greater than or equal to 80%, where Column 13 is your Query Coverage column; OR

	"-globalCov 80 -qLenCol 14 -sLenCol 15"
	to accept hits with inferred global query coverage greater than or equal to 80%, where Column 14 and 15 are your QueryLength and SubjectLength columns respectively.

=head1 USAGE

perl postBlast.pl -b <blast output>

=head2 Options

=head3 Thresholds

	-e	maximum e-value; default=1
	-s	minimum bit score; default=0
	-p	minimum % id; default=0
	-l	minimum alignment length; default=0
	-c	-localCoverage	minimum Query Coverage in %;
		requires a number between 0-100;
		requires '-covCol'
	-covCol	Column number for query coverage in blast output; (start counting at '1')
	-globalCov	minimum Inferred Global Coverage, calculated as [((Alignment_Length - Gaps - Mismatches)/ min(Query_Length, Subject_Length))*100]
				requires a number between 0-100;
				requires '-qLenCol' and '-sLenCol'
	-qLenCol	Query Length column number (start counting at '1')
	-sLenCol	Subject Length column number (start counting at '1')
	-list	Just get me the queries/subj in the list File;
		if you wish to search in subj column for matches; mention with flag "-subj";
		default= will look in query;

=head3 Dependencies

	-subj	look for items in list in the subj col; required if using '-list' AND wish to only check subjects for matches to the list.

=head3 Output

	-o	output file name; default= processID.pb.out
	-r	inverse of output file; default= not printed at all.

=head1 Author

	Sunit Jain, January 2010 : sunitj [AT] umich [DOT] edu
	- updated May 2016

=cut

use strict;
use Getopt::Long;
use List::Util qw(min);

my $setEval=1;
my $setPer=0;
my $setLen=0;
my $setScore=0;
my $out=$$.".pb.out";
my $version="postBlast.pl\tv1.0.12";

my(
$blast,
$list,
$subj,
$queryFasta,
$covCol,
$setGlobalCov,
$setCoverage,
$getCoverage,
$qLenCol,
$sLenCol,
$remainder,
);

GetOptions(
	'b|blast:s'=>\$blast,
	'e|evalue:f'=>\$setEval,
	's|bit_score:f'=>\$setScore,
	'p|per:i'=>\$setPer,
	'l|len:i'=>\$setLen,
	'o|out:s'=>\$out,
	'r|remainder:s'=>\$remainder,
	'list:s'=>\$list,
	'c|cov:i'=>\$setCoverage,
	'globalCov:i'=>\$setGlobalCov,
	'covCol:i'=>\$covCol,
	'qLenCol:i'=>\$qLenCol,
	'sLenCol:i'=>\$sLenCol,
	'f|q|fasta|query:s'=>\$queryFasta,
	'subj'=>\$subj,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system('perldoc', $0); exit;},
);
print "\# $version\n";
die "[ERROR] Blast output file required. See $0 -h for help.\n" if (!$blast);

if (($list) && (!$subj)){ print "# Searching Queries\n";}
elsif (($list) && ($subj)){ print "# Searching Subjects\n";}

## create index file from list;
my %index;
&readListFile if ($list);

## set coverage value
my (%lengths, $printCov);
if ($setCoverage){
	die "Coverage column required. Use -covCol to specify the column number\n" if (! $covCol);
	$printCov=$setCoverage;
	$setCoverage=$setCoverage/100;
}
else{
	$setCoverage=0;
	$printCov=0;
}

## set global coverage value
my ($printGlobalCov);
if ($setGlobalCov){
	die "Query and Subject length column numbers required. Use '-qLenCol' and '-sLenCol' to specify the column numbers\n" if ((! $qLenCol)||(! $sLenCol));
	$printGlobalCov=$setGlobalCov;
	$setGlobalCov=$setGlobalCov/100
}
else{
	$setGlobalCov=0;
	$printGlobalCov=0;
}

open (BLAST, $blast) || die "[ERROR]: $!\n";
open (OUT, ">".$out) || die "[ERROR]: $!\n";
if($remainder){open (REM, ">".$remainder) || die "[ERROR]: $!\n";}
print OUT "#Thresholds: Blast Output=$blast; \%ID=$setPer; aln_len=$setLen; e-value=$setEval; bitScore=$setScore; QueryCoverage=$printCov\%; GlobalCoverage=$printGlobalCov\%\n";
my $header="#Columns: Query\tSubject\tPercent_ID\tAln_Len\tMismatch\tGap\tQuery_Start\tQuery_Stop\tSubject_Start\tSubject_Stop\tE-Value\tBit_Score";
$header.= ($getCoverage ? "\tCoverage" : "");
$header.= ($setGlobalCov ? "\tQuery_Length\tSubject_Length\tGlobal_Coverage" : "");
$header.="\n";
print OUT $header;
while(my $line=<BLAST>){
	next if $line=~ /^#/;
	$line=~ s/\r//;
	chomp($line);
	next unless $line;

	my ($q, $s, $per, $len, $mm, $gap, $q_start, $q_stop, $s_start, $s_stop, $eval, $score, @extra)=split(/\t/, $line);
	my $cov=$setCoverage > 0 ? $extra[$covCol-13] : 0;
	my $qLen=$extra[$qLenCol-13];
	my $sLen=$extra[$sLenCol-13];
	my $gCov = $setGlobalCov > 0 ? (($len-$mm-$gap)/min($qLen,$sLen)) : 0;
	my $item=$subj ? $s : $q;
	if (($per >= $setPer) && ($len >= $setLen) && ($eval <= $setEval) && ($score >= $setScore) && ($cov >= $setCoverage) && ($gCov >= $setGlobalCov)){
		my $printLine="";
		if (($list) && ($index{$item})){
			$printLine = $line
		}
		elsif(($list) && (! $index{$item})){
			next;
		}
		else{
			$printLine = $line;
		}

		if($getCoverage){
			$cov=sprintf("%.3f", ($cov*100));
			$printLine.="\t".$cov;
		}

		if($setGlobalCov > 0){
			$gCov=sprintf("%.3f", ($gCov*100));
			$printLine.="\t".$gCov;
		}

		print OUT $printLine."\n";
	}
	elsif($remainder){
		print REM $line."\n";
	}
}

## Sub-routines ##
sub getLengths{
	open(FASTA, $queryFasta)|| die "[ERROR]: $!\n";
	$/=">";
	while(my $line=<FASTA>){
		chomp $line;
		my ($head, @sequence)= split(/\n/, $line);
		my ($query, @stuff)=split(/\s+/, $head);

		if ($list){
			next unless $index{$query};
		}

		my $seqLen=length(join("", @sequence));
		$lengths{$query}=$seqLen;
	}
	$/="\n";
	close FASTA;
	return;
}

sub readListFile{
	open (LIST, $list)|| die "[ERROR]: $!\n";
	while(my $line=<LIST>){
		next if $line=~ /^#/;
		$line=~ s/\r//;
		chomp($line);
		next unless $line;

		$index{$line}++;
	}
	close LIST;
	return;
}
