#! /usr/bin/perl

use strict;
use Getopt::Long;

my $absolute=50;
my $percent;
my $all;
my $phgGI_list="../phg_gi.list";  # original GI list; obtain this from NCBI.
my $phgGI_data="../phg_gi.data";  # original data file; obtain this from 'blastdbcmd'

my $tmpGI_data="tmp_phg.data";
my $tmpGI_taxa="tmp_taxa.list";

my $testGI_data="test_phg.data";
my $testGI_fasta="test_phg.fasta";
my $testGI_meta="test_phg.meta";
my $testGI_taxa="test_phg.taxa";
my $testGI_distr="test_phg.distr";

GetOptions(
	'p|per:i'=>\$percent,
	'abs|absolute:i'=>\$absolute,
	'all'=>\$all,
	'l|list:s'=>\$phgGI_list,
	'd|data:s'=>\$phgGI_data,
);

#system("blastdbcmd -entry_batch phg_gi.list -db nr -outfmt "%g %T %s" -out phg_gi.data");

# Read the original GI list
# Create an index of all GI numbers.
open(LIST, $phgGI_list); 
my %index;
while(my $line=<LIST>){
	chomp $line;
	next unless $line;

	$index{$line}++;
}	
close LIST;

# Remove the extra GI numbers that 'blastdbcmd' produces based on the original list.
# Create a new temporary data file.
# Create an index of all taxa.
open(DATA, $phgGI_data);
open(TDATA, ">".$tmpGI_data);
my %taxaIndex;
while(my $line=<DATA>){
	chomp $line;
	next unless $line;
	next if $line=~ /^#/;

	my($gi, $taxa, $seq)=split(/\s+/, $line);

	next if $taxa==0;	
	next unless $index{$gi};
	
	print TDATA $line."\n";
	$taxaIndex{$taxa}++;
}
close TDATA;
close DATA;

# Calculate size of subset
my $size;
if ($percent){
	my $totalLines= `wc -l $tmpGI_data`;
	$size= int(($totalLines * ($percent/100))+0.5);
}
elsif($all){
	$size=keys(%taxaIndex);
}
else{
	$size=$absolute;
}

# Randomize
my %rNum;
while (scalar(keys %rNum) < $size){
	my $r=int(rand(keys(%taxaIndex))+1);
	$rNum{$r}++;
}
print "Sample Size:\t".(keys %rNum)."\n";
my ($i, %randomTaxa);
foreach my $t(keys %taxaIndex){
	$i++;
	$randomTaxa{$t}++ if ($rNum{$i});
}
print "Number of Random Taxa:\t".(keys %randomTaxa)."\n";

# Create a data file with the random taxa subset 
open(WDATA, ">".$testGI_data);
open(TDATA, $tmpGI_data);
my %seen;
while(my $line=<TDATA>){
	chomp $line;
	next unless $line;
	next if $line=~ /^#/;

	my($gi, $taxa, $seq)=split(/\s+/, $line);

	next unless $randomTaxa{$taxa};
	next if $seen{$gi};

	print WDATA $line."\n";
	$seen{$gi}++;
}
close WDATA;
close TDATA;
unlink $tmpGI_data;

# Read a data file with the random taxa subset
# Dereplicate
# Data file is of format gi, ncbi_taxa, sequence
open(RDATA, $testGI_data);
my $i=0;
my (%unique, %gi_per_taxa);
while(my $line=<RDATA>){
	chomp $line;
	next unless $line;

	$i++;
	my($gi, $taxa, $seq)=split(/\s+/, $line);

	my $gi_taxa=join("_", $gi, $taxa);
	push(@{$unique{$seq}}, $gi_taxa);

	$gi_per_taxa{$taxa}++;
}
close RDATA;

my $size=keys %unique;

print "Total Sequences:\t".$i."\n";
print "Total Unique Sequences:\t".$size."\n";

# Calculate general Stats per taxa
open(DIST, ">".$testGI_distr);
print DIST "# Taxa\t# of GIs\n";
foreach my $t(keys %gi_per_taxa){
	print DIST $t."\t".$gi_per_taxa{$t}."\n";
}
close DIST;

# Generate alias
# Generate metadata files
open(FASTA, ">".$testGI_fasta);
open(META, ">".$testGI_meta);
open(TAXA, ">".$testGI_taxa);
my %alias_taxa;
my $c=0;
foreach my $u(keys %unique){
	$c++;
	my $alias=sprintf("%010d", $c);
	print FASTA ">$alias\n$u\n";

	print META "$alias\t";
	my %taxon;
	my $i=0;
	foreach my $gt(@{$unique{$u}}){
		print META "$gt\t";
		my ($gi, $taxa)=split(/\_/, $gt);
		$i++;
		$taxon{$taxa}++;
	}
	print META "\n";

	print TAXA "$alias\t";
	foreach my $t(keys %taxon){
#		$alias_taxa{$alias}{$t}= $taxon{$t}/$i;
		
		my $taxaScore=$taxon{$t}/$i;
		print TAXA $t."_".$taxaScore."\t";
	}
	print TAXA "\n";
}
close FASTA;
close META;
close TAXA;
