#!/usr/bin/perl

use strict;

my $list=$ARGV[0];
my $pptDBmeta="/geomicro/data1/COMMON/publicDB/pptDB/pptDB.meta";
my $pptDBgi="/geomicro/data1/COMMON/publicDB/pptDB/pptDB_gi.desc";
my $pptDBtaxa="/geomicro/data1/COMMON/publicDB/pptDB/pptDB_taxa.desc";

open(LIST, $list) || die $!;
my %index;
while(my $line=<LIST>){
	next if $line=~ /^\D/;
	chomp $line;
	next unless $line;
	
	my $alias=sprintf("%010d", $line);
	$index{$alias}++;
}
close LIST;

open(META, $pptDBmeta) || die $!;
open(GI, ">ppt_gi.list") || die $!;
my $counter=0;
my (%gi_alias, %taxa);
while(my $line=<META>){
	next if $line=~ /^#/;
	chomp $line;
	next unless $line;

	my($alias, @gi_taxa)=split(/\t/, $line);

	next unless $index{$alias};
	foreach my $gt(@gi_taxa){
		my($gi, $taxon)=split(/\_/, $gt);
		print GI $gi."\n" unless $gi_alias{$gi};
		$counter++  unless $gi_alias{$gi};
		$gi_alias{$gi}=$alias;
		$taxa{$taxon}++;
	}
}
close GI;

open(TAXA, $pptDBtaxa) || die $!;
open(TOUT, ">$list.taxaReport") || die $!;
print TOUT "#Taxonomy\tNCBI TaxID\tNumber Of Hits\t\% of given subset\tLineage\n";
while (my $line=<TAXA>){
	chomp $line;
	next unless $line;

	my($taxon, $sciName, $div, $lineage)=split(/\t/, $line);
	
	next unless $taxa{$taxon};
	my $perc=($taxa{$taxon}/(keys %taxa)) * 100;
	print TOUT $sciName."\t".$taxon."\t".$taxa{$taxon}."\t".$perc."\t".$lineage."\n";
}
close TAXA;
close TOUT;

open(DESC, $pptDBgi) || die $!;
open(OUT, ">".$list."_anno.desc")|| die $!;
print OUT "#SeqID Number\tDescription\n";
while(my $line=<DESC>){
	next if $line=~ /^#/;
	chomp $line;
	next unless $line;

	my($gi, $desc)=split(/\t/, $line);

	next unless $gi_alias{$gi};
	print OUT $gi_alias{$gi}."\t".$desc."\n";
}
close OUT;
close DESC;


