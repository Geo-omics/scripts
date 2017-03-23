#!/usr/bin/perl

use strict;
use Getopt::Long;

###################################
## Parameters
###################################

my $seqFile;
my $kmerFile=$$.".kmer.out";
my $kps;
my $isFastq;
my $k=4;

GetOptions(
	'i|in:s'=>\$seqFile,
	'k|kmer:i'=>\$k,
	'o|out:s'=>\$kmerFile,
	's|each_seq:s'=>\$kps,
	'fq|fastq'=>\$isFastq,
);

###################################
## Main
###################################

$/= $isFastq ? "@" : ">";
my %kmers;
my %kmersPerSeq;

open (SEQ, $seqFile)|| die "$! : $seqFile\n";
while(my $line=<SEQ>){
	next if ($line=~ /^#/);
	chomp $line;
	$line=~ s/ //;
	next unless $line;

	$isFastq ? &parseFastq($line) : &parseFasta($line);
}
close SEQ;

open (SUMM, ">".$kmerFile);
my %seen;
my @kmerArray;
while(my($kmer, $count)=each(%kmers)){
	next if $seen{$kmer};
	
	my $rc_kmer=&rev_comp($kmer);
	
	# if rev_comp exists, combine counts and mark with '*'.
	my $totalCount= $kmers{$rc_kmer} ? ($kmers{$kmer} + $kmers{$rc_kmer})."\t\*" : $kmers{$kmer} ;

	# Print to Output file
	print SUMM $kmer."\t".$totalCount."\n";
	push(@kmerArray, $kmer);
	$seen{$kmer}++;
	$seen{$rc_kmer}++;
}
close SUMM;
unlink %kmers;
unlink %seen;

if($kps){
	open (OUT, ">".$kps);
	print OUT "#SeqNames.\t";
	foreach my $km(@kmerArray){ print OUT $km."\t"; }
	print OUT "\n";

	foreach my $desc(keys %kmersPerSeq){
		print OUT $desc."\t";
		foreach my $km(@kmerArray){
			my $rc_km=&rev_comp($km);
			my $totalCount= $kmersPerSeq{$desc}{$km} + $kmersPerSeq{$desc}{$rc_km};
			$totalCount = $totalCount ? $totalCount : 0;
			print OUT $totalCount."\t";
		}
		print OUT "\n";
	}
	close OUT;
}
unlink %kmersPerSeq;
exit;

###################################
## Sub-Routines
###################################

sub parseFastq{
	my 	$line=shift;
	my($seqDesc,$seq,$qualDesc,$qual)=split(/\n/, $line);
	$seq=~ s/ //g;
	&getKmers($seq, $seqDesc);
	return;
}

sub parseFasta{
	my 	$line=shift;
	my($seqDesc,@sequence)=split(/\n/, $line);
	my $seq=join("",@sequence);
	$seq=~ s/ //g;
	&getKmers($seq, $seqDesc);
	return;
}

sub getKmers{
	my $seq = shift;
	my $desc=shift;
	$seq=uc($seq);
	my $windows=length($seq) - $k;
	
	for (my $pos=0; $pos <= $windows; $pos++){
		my $kmer=substr $seq, $pos, $k;
		$kmers{$kmer}++;
		$kmersPerSeq{$desc}{$kmer}++;
	}
	return;
}

sub rev_comp{
	my $seq=shift;
	$seq=uc($seq);
	$seq=reverse($seq);
	$seq=~ tr/ATGCN/TACGN/;
	return $seq;
}
