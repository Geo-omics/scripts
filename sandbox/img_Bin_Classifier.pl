#!/usr/bin/perl

=head1 DESCRIPTION

	img_Bin_Classifier.pl -- Use the IMG taxonomic classification of contigs/scaffolds to get the taxonomic makeup of each bin.

=head1 USAGE

	perl img_Bin_Classifier.pl -img consolidated.tsv -conf binning_confidence.conf -prefix outputName

=head2 Options

	-img	<CHAR>	Consolidated JGI-IMG data from consolidateJGIdata.pl script version:0.1.3 and up; previous versions may have issues with column orientation.
	-conf    <CHAR>	Concatenated file of Confidence values for all bins obtained from getClassFasta.pl script
	-prefix    <CHAR>	Prefix for output files with taxonomic distribution for each bin
	-top	<INT>	Top X taxonomies per bin; [default=5]
	-virus    <INT>	Value between [0-100]	Minimum percent id for a viral taxa; [default=30]
	-bact    <INT>	Value between [0-100]	Minimum percent id for a bacterial taxa; [default=50]
	-other    <INT>	Value between [0-100]	Minimum percent id for any other taxa; [default=50]
	-min_conf    <FLOAT>	Value between [0-1]	Minimum confidence value for a contig/scaffold to be considered a part of the bin; [default=0.6]
	-img_names	<BOOLEAN>	Use IMG contig names. If your conf file contains IMG assigned contig names
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

my ($IMG, $names, $cls, $conf, $useIMGName);
my $top = 5;
my $prefix=$$;
my $help;
my $version="img_Bin_Classifier.pl\tv0.0.11b";
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
	'p|prefix:s'=>\$prefix,
	'img_names'=>\$useIMGName,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

my $taxaSumm=$prefix.".taxa";
my $RNA=$prefix.".rrna";
my $LOG=$prefix.".multibin";
my $BINs=$prefix.".bins.list";
my %binIndex;
if ($conf){
	open(CONF, "<".$conf)||die $!;
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

my (%taxaData, %binSize, %binData, %seen, %rnaBin, %no_taxa);
open(CONS, "<".$IMG)|| die $!;
open(BINS, ">".$BINs)|| die $!;
print BINS "#Bins\tContigs/Scaffolds\n";
my $TMP=$$.".tmp";
open(TMP, ">".$TMP)|| die $!;
open(LOG, ">".$LOG)|| die $!;
while(my $line=<CONS>){
	chomp $line;
	next unless $line;
	next if $line=~ /^#/;
	
	my (@data)=split(/\t/, $line);
	my $contig=$useIMGName ? strip(lc($data[1])) : strip(lc($data[2]));
	my $contigLen=strip($data[4]);
	my $bin;
	if ($conf){
		$bin=$seen{$contig} ? $seen{$contig} : getBin($contig);
	}
	else{
		$bin=strip(lc($data[0]));
	}

	my $geneType=strip(lc($data[7]));
	if ($geneType eq "rrna"){
		my $geneDesc=strip($data[15]);
		my $geneLen=strip($data[10]);
		print TMP $bin."\t".$geneType."\t".$contig."\t".$geneDesc."\t".$geneLen."\n";
		$rnaBin{$bin}{"Numbers"}++;
		$rnaBin{$bin}{"16s"}++ if ($geneDesc=~ /16S rRNA/);
		$binSize{$bin}++;
		next;
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
		$binData{$bin}{"Len"}+=$contigLen unless $seen{$contig};
		$binSize{$bin}++;
		$no_taxa{$geneType}++;
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
	$binData{$bin}{"Len"}+=$contigLen unless $seen{$contig};
	$binSize{$bin}++;
}
close CONS;
close TMP;
close LOG;
close BINS;

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
print SUMM "# Bin\tTotal_Genes\trRNA\tSize";
print SUMM "\t".$_."\t".$_." \%" foreach(@colOrder);
print SUMM "\n";
my %bins_with_16s;
foreach my $bin(keys %printMatrix){
	for(my $row=0; $row < $top; $row++){
		print SUMM $bin."\t".$binSize{$bin}."\t".$rnaBin{$bin}{"Numbers"}."\t".$binData{$bin}{"Len"};
		if($rnaBin{$bin}{"16s"}){ print SUMM "*"; $bins_with_16s{$bin}++;}
		foreach my $t(@colOrder){
			if($printMatrix{$bin}{$t}[$row]){
				my $value=$printMatrix{$bin}{$t}[$row];
				my $pid=sprintf("%.3f",(($taxaData{$bin}{$t}{$value}/$binSize{$bin}) * 100));
				print SUMM "\t".ucfirst($value)."\t".$pid."\%";
			}
			else{
				print SUMM "\t\t";
			}
		}
		print SUMM "\n";
	}
}

print "# ".scalar(keys %bins_with_16s)." out of ".scalar(keys %binSize)." Bins had rRNAs\n";
print "# These Gene types had no taxonomic Information and were clubbed into the 'Unknown' category:\n";
print $_."\t".$no_taxa{$_}."\n" foreach(keys %no_taxa);

system("echo '# BIN\tGene_Type\tContig/Scaffold\tDescription\tLength' > $RNA");
system("sort -V -k 3,1 $TMP >> $RNA");

unlink $TMP;

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
		print LOG "[WARNING]\t".$Contig." was put in multiple bins.\n";
		foreach my $k(@sortedKeys){
			print LOG "[WARNING]\t".$k."\t".$binIndex{$contig}{$k}."\n";
		}
	
		if($binIndex{$contig}{$sortedKeys[0]} >= $setMinConf){
			my $bin=strip($sortedKeys[0]);
			$seen{$contig}=$bin;
			print BINS $bin."\t".ucfirst($contig)."\n";
			return $bin;
		}
		else{
			$seen{$contig}="Unclassified";
			print LOG "[WARNING]\tNothing met the minimum confidence criterion(".$setMinConf.")\n";
			print LOG "[WARNING]\t".$Contig." Excluded from analysis\n";
			return "Unclassified";
		}
	}
	elsif($binIndex{$contig}{$sortedKeys[0]} >= $setMinConf){
		my $bin=strip($sortedKeys[0]);
		$seen{$contig}=$bin;
		print BINS $bin."\t".ucfirst($contig)."\n";
		return $bin;
	}
	else{
		$seen{$contig}="Unclassified";
		my $k=strip($sortedKeys[0]);
		print LOG "[WARNING]\t".$Contig." failed to meet minimum confidence criterion(".$setMinConf.")\n";
		print LOG "[WARNING]\t".$k."\t".$binIndex{$contig}{$k}."\n";
		print LOG "[WARNING]\t".$Contig." Excluded from analysis\n";
		return "Unclassified";
	}
}
