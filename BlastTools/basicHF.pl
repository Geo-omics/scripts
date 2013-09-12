#!/usr/bin/perl

=head1 NAME
	
	basicHF - Basic no-frills Homolog Finder.
	
=head2 Dependencies

	Blast 2.2.20 or below (old)

=head2 USAGE

	perl basicHF.pl [-i: list of all files you wish to compare, the seed should be the first one on the list]
	[-p: Blast percent identity] [-a: alignment length ratio (0-1)]	[-e: E-Value] [-t: file type (protein or nucleotide)]
	[-h: help (this page)]
	
=head2 DEFAULTS
	
	If you dont mention any of the following, this is what will be assumed.
	-p	blast % id= 50%
	-a	alignment length ratio= 0.75
	-e	E-Value= -4;
	-t	File type= 'p' (protein).
	
=head2 AUTHOR
	
	Sunit Jain, December 2010
	sunitj_at_umich_dot_edu

=cut

use strict;
use Getopt::Long;
use Benchmark;


my $start=new Benchmark;

## Default Values ##
my $list;
my $per=50;
my $aLen=0.75;
my $eVal=-4;
my $fType='p';

GetOptions(
	'i|input:s'=>\$list,
	'p|percent:f'=>\$per,
	'a|len:f'=>\$aLen,
	'e|eval:f'=>\$eVal,
	't|fileType:s'=>\$fType,
	'h|help'=> sub{system('perldoc', $0); exit;},
);

if (!$list){system('perldoc', $0); exit;}

## MAIN ##

uc($fType);
open (LIST, $list) || die "[ERROR] $list:$!\n";
my @fileList=<LIST>;
close LIST;

chomp(@fileList);

my($dbFormat,$blastType);
if (lc($fType) eq 'p'){
	$dbFormat='T';
	$blastType='blastp';
}
elsif(lc($fType) eq 'n'){
	$dbFormat='F';
	$blastType='blastn';
}
	
my $num=0;
my @edList;
foreach my $f(@fileList){
	my $ed=prepFasta($f, $num, $dbFormat);
	push(@edList, $ed);
	$num++;
}
my $seed=@edList[0];

my $totalGenomes=@fileList;

my %qLens=getLengths($seed);

open(NAME, ">".$seed."P".$per."L".$aLen."E".$eVal."B121.names");
open(NUM, ">".$seed."P".$per."L".$aLen."E".$eVal."B121.tsv");

## PARAMETERS ##
print NAME "\#E-Value:\t".$eVal."\n";
print NAME "\#\% ID:\t".$per."\n";
print NAME "\#Alignment Length Ratio:\t".$aLen."\n";

print NUM "\#E-Value:\t".$eVal."\n";
print NUM "\#\% ID:\t".$per."\n";
print NUM "\#Alignment Length Ratio:\t".$aLen."\n";

print NAME "\#Ref\t";
print NUM "\#Ref\t";

my %homologs;
for (my $i=0; $i<$totalGenomes; $i++){
	my %qsHits;
	my $blastOut="Qvs".$i.".bout";
	print NAME $fileList[$i]."\t";
	print NUM $fileList[$i]."\t";
	print "Blast:\t".$seed."\tvs\t".$fileList[$i]."\n";

	my %sLens=getLengths($edList[$i]);	
	my $subjSize=keys(%sLens);
	if ($subjSize>=1){
		system ("blastall -p ".$blastType." -e 1e".$eVal." -m 9 -a 8 -i ".$seed." -d ". $edList[$i]." -o ".$blastOut);

		my %qsHits=parseBlastOutput($blastOut, \%qLens, \%sLens, $per, $aLen);
	
		my %seen;
		while(my ($k, $value)=each(%qsHits)){
			chomp(@{$value});
			foreach my $v(@{$value}){
				push(@{$homologs{$k}}, $v) unless ($seen{$v} eq $k);
				$seen{$v}=$k;
			}
		}
	}
}

my %numbers;
print NAME "\n";
while (my($k, $value)=each(%homologs)){
	my($n, $num)=split(/\@/, $k);
	print NAME $n."\t";
	foreach my $v(@{$value}){
		my($name, $pos)=split(/\@/, $v);
		$numbers{$n}{$pos}++;
		print NAME $name."\t";
	}
	print NAME "\n";
}
close NAME;

print NUM "Total\n";
foreach my $k(keys %numbers){
	print NUM $k."\t";
	my $sum=0;
	for (my $i=0; $i<$totalGenomes; $i++){
		if ($numbers{$k}{$i}){
			print NUM $numbers{$k}{$i}."\t";
		}
		else{
			print NUM "0\t";
		}
		$sum +=$numbers{$k}{$i};
	}
	print NUM "$sum\n";
}
close NUM;

#print "Removing temporary files...\n";
system ("rm *.phr");
system ("rm *.pin");
system ("rm *.psq");
system ("rm *_ed.faa");
system ("rm *.bout");


# Timing
my $end=new Benchmark;
my $time=timediff($end,$start);
print "For $totalGenomes Genomes the script took: ". timestr($time, 'all')."\n";


## Edit the headers of all fasta files to make-up for BLAST's stupid way of dealing with subject names ##
sub prepFasta{
	
	my $in=shift;
	my $gNum=shift;
	my $format=shift;

	my %wg;

	chomp($in);
	open (FAA, $in) || die "ERROR: $in\n".$!;
	$/= ">";
	while (my $a = <FAA>) {
		chomp $a;
		next unless $a;
		my ($query, @seqs) = split (/\n/, $a);
		my ($numbers, $taxa);
		if ($query=~ m/^gi\|/){
			(my $giTag,my  $gi, $taxa,my $name)=split(/\|/, $query);
			$numbers=$gi."\|".$taxa;
		}
		else{
			my ($num,$name)=split(/\s/, $query);
			$numbers=$num;
		}
		my $nuQuery=$numbers."@".$gNum;
		my $seq = join ("", @seqs);
		chomp($nuQuery, $seq);
		$wg{$nuQuery}=$seq;
	}
	close FAA;
	$/= "\n";
	
	## print out edited version.
	my($fName, @ext)=split(/\./, $in);
	my $out=$fName."_ed.faa";
	open (ED, ">$out");
	while (my($k, $v)=each(%wg)){
		chomp($k, $v);
		print ED ">".$k."\n";
		print ED $v."\n";
	}
	
#	print "\tMaking Blast DB for $out\n";
	system ("formatdb -i ".$out." -p ".$format);
#	print "\tDONE!!\n\t\tTotal Contigs:\t".keys(%wg)."\n\n";

	return $out;
}

## Parse blast output, keep the useful info and discard the rest ##
sub parseBlastOutput {
## Sub-Routine Input i.e. BLAST Output Format: m 8 ##
#	0		1		2	3		4		5		6		7		8		9		10		11
#	query	sbjct	%id	a_len	m_cnt	g_cnt	q_start	q_end	s_start	s_end	e_val	b_score

##	Sub-routine Output hash structure ##
#	key	0	1	2		
#	query	sbjct	%id	a_len

	my $bOut = shift;
	my $qL=shift;
	my $sL=shift;
	my %qLens=%{$qL};
	my %sLens=%{$sL};
	my $setPer= shift;
	my $setLen= shift;
	
	my %hits;

	open (BOUT, "$bOut") || die "ERROR: $bOut \n:".$!;
	while(my $line=<BOUT>){
		next if ($line=~ m/^\#/);
		chomp($line);
		my @blastOut=split(/\t/, $line);
		chomp(@blastOut);
		my @everything_else=@blastOut[1..3];
		my $query=$blastOut[0];
		my $subj=$blastOut[1];
		my $per=$blastOut[2];
		my $aLen=$blastOut[3];
		chomp($query,$subj,$per,$aLen);
# Get all hits	
		my $tempLen;
		if ($qLens{$query}){
			# print "\|";
			my $subjLen=$sLens{$subj};
			my $queryLen=$qLens{$query};
			# print $subjLen."\t".$queryLen."\n";
			if ($subjLen==$queryLen || $subjLen < $queryLen){
				$tempLen=$aLen/$subjLen;
			}
			else{
				$tempLen=$aLen/$queryLen;
			}
		}
		if (($tempLen > $setLen) && ($per > $setPer)){	
			chomp($tempLen);
			my ($sName, $gNum)=split(/\@/, $subj);
			my $subjName=$sName."\(".$per."\%\)\@".$gNum;
			#print $subjName."\t";
			push(@{$hits{$query}},$subjName);
		}
	}
	close BOUT;
	return %hits;
}

## Calculate length of all sequences, for alignment Length ratio ##
sub getLengths{
	my $in=shift;
	chomp($in);
	open (FAA, $in) || die "ERROR: $in\n".$!;
	$/= ">";
	my %lens;
	while (my $a = <FAA>) {
		chomp $a;
		next unless $a;
		my ($query, @seqs) = split (/\n/, $a);
		my $seq = join ("", @seqs);
		my $len = length($seq);
		$lens{$query}=$len;
		@seqs=();
	}
	close FAA;
	$/= "\n";
	return %lens;
}
