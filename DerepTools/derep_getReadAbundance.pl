#! /usr/bin/perl

use strict;
use Getopt::Long;

my ($list,$bam, $stats, $fwdClust, $revClust);
my $out=$$.".abundance";

GetOptions(
	'b|bam=s'=>\$list,
	'fwd=s'=>\$fwdClust,
	'rev=s'=>\$revClust,
	'stats=s'=>\$stats,
	'o|out:s'=>\$out,
	'h'=>sub{system('perldoc', $0); exit;},
);
my %clust;

open (FCLUST, $fwdClust) || die "[ERROR] $fwdClust:$!\n";
while(my $line=<FCLUST>){
	next if ($line=~ /^#/);
	chomp($line);
	$line=~ s/\r//;
	next unless $line;

	my($cNum, $size, $rep, @seqNames)=split(/\t/, $line);
	my ($name, $strand)=split(/\s/, $rep);
	$name=~ s/^@//;

	$clust{$name}=$size;
}
close FCLUST;

open (RCLUST, $revClust) || die "[ERROR] $revClust:$!\n";
while(my $line=<RCLUST>){
	next if ($line=~ /^#/);
	chomp($line);
	$line=~ s/\r//;
	next unless $line;

	my($cNum, $size, $rep, @seqNames)=split(/\t/, $line);
	my ($name, $strand)=split(/\s/, $rep);
	$name=~ s/^@//;

	if ($clust{$name}){
		$clust{$name}+=$size;
	}
	else{
		$clust{$name}=$size;
	}
}
close RCLUST;

print "All Clusters read into memory...\n";

open (STATS, $stats);
my %index;
while(my $line=<STATS>){
	chomp $line;
	next unless $line;

	my ($name, $len, $mapped, $unmapped)=split(/\t/, $line);
	$name=~ s/\|/\\\|/;

	$index{$name}=$mapped;
}
close STATS;

print "All Stats read into memory...\n";

open(OUT, ">".$out);
print OUT "Gene\tDerepMapped\tTotalMapped\n";
my $allReads=0;
foreach my $i(keys %index){
	my @list=`samtools view -F0x4 $bam $i | cut -f 1`;
	chomp @list;

	my $totalMapped=0;
	foreach my $l(@list){
		$totalMapped+= $clust{$l};
		$allReads+=$clust{$l};
	}
	my $derepMapped= $index{$i};
	print OUT $i."\t".$derepMapped."\t".$totalMapped."\n";
}

print "Total Reads mapped: $allReads\n";
