#!/usr/bin/perl

=head1 DESCRIPTION

	img_Bin_Classifier.pl -- Use the IMG taxonomic classification of contigs/scaffolds to get the taxonomic makeup of each bin.

=head1 USAGE

	perl img_Bin_Classifier.pl -img consolidated.tsv -conf binning_confidence.conf -out bin.taxa

=head2 Options

	-img	<CHAR>	Consolidated JGI-IMG data from consolidateJGIdata.pl script version:0.1.3 and up; previous versions may have issues with column orientation.
	-conf    <CHAR>	Concatenated file of Confidence values for all bins obtained from getClassFasta.pl script
	-out    <CHAR>	Output with taxonomic distribution for each bin
	-top	<INT>	Top X taxonomies per bin; [default=5]
	-virus    <INT>	Value between [0-100]	Minimum percent id for a viral taxa; [default=30]
	-bact    <INT>	Value between [0-100]	Minimum percent id for a bacterial taxa; [default=50]
	-other    <INT>	Value between [0-100]	Minimum percent id for any other taxa; [default=50]
	-min_conf    <FLOAT>	Value between [0-1]	Minimum confidence value for a contig/scaffold to be considered a part of the bin; [default=0.6]

	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Authors

	Sunit Jain, (Tue Nov 26 14:32:26 EST 2013)
	sunitj [AT] umich [DOT] edu

	Daniel Marcus
	dnmarc [AT] umich [DOT]	edu

=head1 License

	This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

	This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;

my ($IMG, $names, $cls, $map, $taxaSumm, $conf);
my $top = 5;
my $taxaSumm=$$.".bin.taxa";
my $help;
my $version="img_Bin_Classifier.pl\tv0.0.6b";
my $vpid=30;
my $bpid=50;
my $other_pid=50;
my $setMinConf=0.6;
GetOptions(
	'img|jgi:s'=>\$IMG,
	'conf:s'=>\$conf,
	'virus:i'=>\$vpid,
	'bact:i'=>\$bpid,
	'other:i'=>\$other_pid,
	'min_conf:i'=>\$setMinConf,
	'top:i'=>\$top,
	'o|out:s'=>\$taxaSumm,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

my $rRNA=$taxaSumm.".rrna";
my %binIndex;
if ($conf){
	open(CONF, ">".$conf)||die $!;
	while(my $line=<CONF>){
		chomp $line;
		$line=strip($line);
		next if ($line=~ /^#/);
		next unless $line;
	
		my($bin, $contig, $confVal)=split(/\t/, lc($line));
		$binIndex{$contig}{$bin}=$confVal;
	}
	close CONF;
}

my (%taxaData, %binSize, %unfit, %seen);
open(CONS, "<".$IMG)|| die $!;
open(RNA, ">".$RNA)|| die $!;
while(my $line=<CONS>){
	chomp $line;
	next unless $line;
	next if $line=~ /^#/;
	
	my (@data)=split(/\t/, $line);
	my $rRNA=$data[7];
	my $contig=strip(lc($data[2]));
	my $bin;
	if ($conf){
		$bin=$seen{$contig} ? $seen{$contig} : getBin($contig);
	}
	else{
		$bin=strip(lc($data[0]));
	}
	my $geneID=strip(lc($data[6]));
	my $taxa_pid=strip(lc($data[13]));
	my $taxa=lc($data[14]);
	
	if (! $taxa){
		$taxaData{$bin}{"Domain"}{"Unknown"}++;
		$taxaData{$bin}{"Phyla"}{"Unknown"}++;
		$taxaData{$bin}{"Class"}{"Unknown"}++;
		$taxaData{$bin}{"Order"}{"Unknown"}++;
		$taxaData{$bin}{"Family"}{"Unknown"}++;
		$taxaData{$bin}{"Genus"}{"Unknown"}++;
		$taxaData{$bin}{"Species"}{"Unknown"}++;
		$binSize{$bin}++;
		next;
	}
	# Domain\tPhyla\tClass\tOrder\tFamily\tGenus\tSpecies
	my($domain,$phyla,$class,$order,$family,$genus,$species)=split(/\;/, $taxa);
	
	$domain=strip($domain);
	if($domain=~ /virus/i){
		next if ($taxa_pid < $vpid);
	}
	elsif($domain=~ /bacteria/i){
		next if ($taxa_pid < $bpid);
	}
	else{
		next if ($taxa_pid < $other_pid);
	}
	
	$phyla=strip($phyla);
	$class=strip($class);
	$order=strip($order);
	$family=strip($family);
	$genus=strip($genus);
	$species=strip($species);
	
	$taxaData{$bin}{"Domain"}{$domain}++;
	$taxaData{$bin}{"Phyla"}{$phyla}++;
	$taxaData{$bin}{"Class"}{$class}++;
	$taxaData{$bin}{"Order"}{$order}++;
	$taxaData{$bin}{"Family"}{$family}++;
	$taxaData{$bin}{"Genus"}{$genus}++;
	$taxaData{$bin}{"Species"}{$species}++;
	if ($geneType eq "rRNA"){
		print RNA $bin."\t".$contig."\t".$taxa."\n";
	}
	$binSize{$bin}++;
}
close CONS;
close RNA;

my @colOrder=qw(Domain Phyla Class Order Family Genus Species);
my %printMatrix;
foreach my $bin(keys %taxaData){
	foreach my $col(@colOrder){
		my @sortedKeys= sort{ $taxaData{$bin}{$col}{$b} <=> $taxaData{$bin}{$col}{$a} } keys %{$taxaData{$bin}{$col}};
		foreach my $k(@sortedKeys){
			push(@{$printMatrix{$bin}{$col}},$k);
		}
	}
}

open(SUMM, ">".$taxaSumm)||die $!;
print SUMM "# Bin\tTotal_Genes";
print SUMM "\t".$_."\t".$_." \%" foreach(@colOrder);
print SUMM "\n";
foreach my $bin(keys %printMatrix){
	for(my $row=0; $row < $top; $row++){
		print SUMM $bin."\t".$binSize{$bin};
		foreach my $t(@colOrder){
			if($printMatrix{$bin}{$t}[$row]){
				my $value=ucfirst($printMatrix{$bin}{$t}[$row]);
				my $pid=sprintf("%.3f",(($taxaData{$bin}{$t}{$value}/$binSize{$bin}) * 100));
				print SUMM "\t".$value."\t".$pid."\%";
			}
			else{
				print SUMM "\t\t";
			}
		}
		print SUMM "\n";
	}
}

sub strip{
	my $data=shift;
	chomp $data;
	$data=~ m/^\s+/;
	$data=~ m/\s+$/;
	return $data;
}

sub getBin{
	my $contig=shift;
	return "Unclassified" if (! $binIndex{$contig});	
	my $Contig=ucfirst($contig);
	# sort by value
	my @sortedKeys= sort{ $binIndex{$contig}{$b} <=> $binIndex{$contig}{$a} } keys %{$binIndex{$contig}};
	if((scalar @sortedKeys) > 1 ){
		print "[WARNING]\t".$Contig." was put in multiple bins.\n";
		foreach my $k(@sortedKeys){
			print "[WARNING]\t".$k."\t".$binIndex{$contig}{$k}."\n";
		}
	
		if($binIndex{$contig}{$sortedKeys[0]} >= $setMinConf){
			my $bin=strip($sortedKeys[0]);
			$seen{$contig}=$bin;
			return $bin;
		}
		else{
			$seen{$contig}="Unclassified";
			print "[WARNING]\tNothing met the minimum confidence criterion(".$setMinConf.")\n";
			print "[WARNING]\t".$Contig." Excluded from analysis\n";
			return "Unclassified";
		}
	}
	elsif($binIndex{$contig}{$sortedKeys[0]} >= $setMinConf){
		my $bin=strip($sortedKeys[0]);
		$seen{$contig}=$bin;
		return $bin;
	}
	else{
		$seen{$contig}="Unclassified";
		my $k=strip($sortedKeys[0]);
		print "[WARNING]\t".$Contig." failed to meet minimum confidence criterion(".$setMinConf.")\n";
		print "[WARNING]\t".$k."\t".$binIndex{$contig}{$k}."\n";
		print "[WARNING]\t".$Contig." Excluded from analysis\n";
		return "Unclassified";
	}
}
