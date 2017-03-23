#!/usr/bin/perl

=head1 NAME
	
	uClustHomology - Does a combined UClust on the list of files that you provide.
	Dependencies: USearch 5.0
	
=head2 USAGE
	
	perl uClustHomology.pl [-i: list of all files you wish to compare, the seed should be the first one on the list]
		[-p: uClust percent identity]
		[-t: type of files provided, 'p' for proteins and 'n' for nucleotides]
		[-o: name of the output file]
		[-h: help; brings this screen]
						   
=head2 OPTIONS

	-i	REQUIRED
	-p	OPTIONAL		Default: 95
	-t	OPTIONAL (p/n)	Default: p (protein)
	-o	OPTIONAL		Default: queryFileName.tsv
	-h	This screen.
	
=head2 AUTHOR

	Sunit Jain, Dec 2010
	sunitj_at_umich_dot_edu
	
=cut

use strict;
use Getopt::Long;
use Benchmark;

sub scriptHelp {
	system('perldoc', $0);
	exit;
}

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
	return $out;
}

sub uClust{
	my $genome=shift;
	my $iden=shift;
	my $uType=shift;
	my $id=$iden/100;
	my ($fName, $ext)=split (/\./, $genome);
	my $sorted=$fName.".sorted";
	my $out=$fName."_".$iden.".count";
	system("uclust --mergesort ".$genome." --output ".$sorted." --".$uType." --quiet");
	system("uclust --input ".$sorted." --uc ".$out." --id ".$id." --".$uType." --quiet");
	my (%c)=clustCount($out);
	return (%c);
}

sub clustCount{
	my $genome=shift;
	open (CLUST, "$genome") || die "[ERROR]: $genome\n$!";
	my %clusters;
	while (my $line=<CLUST>){
		unless ($line=~/^\#|^C/){
			chomp $line;
			my @stuff=split(/\t/, $line);
			my ($number, @etc)=split (/\s/, $stuff[-2]); #$number=name of the seq;
			if ($line=~/^S/){
				#$stuff[1] =clusterNumber;
				push (@{$clusters{$stuff[1]}}, $number);
			}	
			elsif($line=~/^H/){
				push (@{$clusters{$stuff[1]}}, $number);
			}	
		}
	}
	
	my($fName, $ext)=split (/\./, $genome);
	my $out=$fName.".clust";
	open (OUT, ">$out");
	print OUT "\#CNum\tTotal\tClusteredSeqs\n";

	my $totalGenes=0;
	foreach my $key(keys %clusters){
		my $arraySize= scalar (@{$clusters{$key}});
		my $everything= join("\t", @{$clusters{$key}});
		print OUT $key."\t".$arraySize."\t".$everything."\n";
		$totalGenes=$totalGenes+$arraySize;
	}
	print OUT "\#TotalGenes:\t".$totalGenes."\n";
	return (%clusters);	
}

sub makeItReadable{
	my $clusterHash=shift;
	my $totalGenomes=shift;
	my %clusters=%{$clusterHash};

	my %numbers;
	my %queryGenome;

	while(my($k, $value)=each(%clusters)){
		foreach my $v(@{$value}){
			my ($name, $gNum)=split(/\@/, $v);
			$numbers{$k}{$gNum}++;
			if ($gNum==0){
				$queryGenome{$k}=$name if (!$queryGenome{$k});
			}
		}
	}
	return(\%numbers, \%queryGenome);
}

## MAIN ##

## Start Timer
my $start=new Benchmark;

## Default Values
my $listFile;
my $pid=95;
my $fileType='p';
my $out;

## Get Options
GetOptions(
	'i|list:s'=>\$listFile,
	'p|pid:f'=>\$pid,
	't|type:s'=>\$fileType,
	'o|out:s'=>\$out,
	'h|help'=> sub{system('perldoc', $0);exit;},
);

## Parse Options
&scriptHelp if (!$listFile);

my $fType;
if (lc($fileType) eq 'p'){
	$fType='amino';
}
elsif(lc($fileType) eq 'n'){
	$fType='nucl';
}
else{
	&scriptHelp
}

open (LIST, $listFile)||die "[ERROR]$listFile: $! \n";
my @list=<LIST>;
chomp(@list);
close LIST;

if (!$out){
	my($fName, @etc)=split(/\./, $list[0]);
	$out=$fName.".tsv";
}

open (NUM, ">".$out)|| die "[ERROR]$out: $! \n";
print NUM "\# \% ID\t".$pid."\n";
print NUM "\# File Type\t".$fType."\n";
print NUM "\# Query Ref\t".$list[0]."\n";
print NUM "\# CNum\t$list[0]\t";

my $count=0;
my @edList;
foreach my $file(@list){
	next if ($file=~ m/^\#|^\s/);
	my $edFile=prepFasta($file, $count);
	print NUM $file."\t";
	$count++;
	push(@edList, $edFile);
}

my $concatenate=join(" ", @edList);
my $combined=$out.".cfasta";
print "Concatenating ".@edList." Files...\n";
system("cat ".$concatenate." > ".$combined);

my %uClustOut=uClust($combined, $pid, $fType);

my($h1, $h2)=makeItReadable(\%uClustOut, $count);
my %numberd=%{$h1};
my %qGenome=%{$h2};

undef %uClustOut;
print NUM "Total\n";
while (my($k, $value)=each(%numberd)){
	print NUM "$k\t";
	if ($qGenome{$k}){
		print NUM $qGenome{$k}."\t";
	}
	else{
		print NUM "0\t";
	}
	my $sum=0;
	for (my $i=0; $i<$count; $i++){
		if ($numberd{$k}{$i}){
			print NUM $numberd{$k}{$i}."\t";
			$sum+=$numberd{$k}{$i};
		}
		else{
			print NUM "0\t";
		}
	}
	print NUM "$sum\n";
}

system("rm *.count");
system("rm *.sorted");
system ("rm *_ed.faa");

# stop timer
my $stop=new Benchmark;
my $time=timediff($stop,$start);
print "For $count Genomes the Process took ". timestr($time, 'all')."\n";
