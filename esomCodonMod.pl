#!/usr/bin/perl
=head1 USAGE

	perl esomCodonMod.pl -lrn file.lrn -o outputFile.lrn

=cut

use strict;
use Getopt::Long;

my ($lrn, $tri);
my $out=$$."_".$lrn;
my $version=0.2.0

GetOptions(
	'lrn=s'=>\$lrn,
	'o|out:s'=>\$out,
	'tri'=>\$tri,
	'v|version'=>sub{print STDERR $0."\tversion:".$version."\n"},
	'h|help'=>sub{system('perldoc', $0); exit;},
);

&help  if (! $lrn);
sub help{system('perldoc', $0); exit;}

my @removeCodons=qw (ATG TAG TAA TGA);
my @nucl=qw(A T G C);

my %removeTetra;
foreach my $c(@removeCodons){
	foreach my $n(@nucl){
		$removeTetra{$n.$c}++;
		$removeTetra{$c.$n}++;		
	}
}
if($tri){	
	foreach (@removeCodons){	$removeTetra{$_}++; }
}

print "Tetramers Removed:\t".keys(%removeTetra)."\n";
#foreach my $k(keys %removeTetra){
#	print $k."\n";
#}
#print "CodonMod says:".$lrn."\n";
open(LRN, $lrn) || die $!;
open(OUT, ">".$out) || die $!;
my @codonOrder;
while(my $line=<LRN>){
	chomp $line;
	next unless $line;
	
	my $pos=0;
	if ($line=~ /^\% Key/){
		@codonOrder=split(/\t/, $line);
		foreach (@codonOrder){	print OUT $_."\t" unless $removeTetra{$_}};
	}
	elsif($line=~ /^\d/){
		my $thisLine;
		my @frequencies=split(/\t/, $line);
		foreach my $freq(@frequencies){
			next if ($removeTetra{$codonOrder[$pos]});
			$thisLine.=$freq."\t";
			$pos++;
		}
		$thisLine=~ s/\t$/\n/;
		print OUT $thisLine;
	}
	else{
		print OUT $line."\n";
	}
}

exit 0;
