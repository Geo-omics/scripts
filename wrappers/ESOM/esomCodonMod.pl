#!/usr/bin/perl
=head1 USAGE

	perl esomCodonMod.pl -lrn file.lrn -o outputFile.lrn

=cut

use strict;
use Getopt::Long;

my ($lrn, $tri);
my $out=$$."_".$lrn;
my $version=0.2.0;

GetOptions(
	'lrn=s'=>\$lrn,
	'o|out:s'=>\$out,
	'tri'=>\$tri,
	'v|version'=>sub{print STDERR $0."\tversion:".$version."\n";},
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

#print "Possible Tetramers that can be Removed:\t".keys(%removeTetra)."\n";

open(LRN, $lrn) || die $!;
my (@codonOrder);
my ($cols, $secondPart, $firstLine, $removed);
while(my $line=<LRN>){
	chomp $line;
	next unless $line;
	if ($line=~ /^\% Key/){
		@codonOrder=split(/\t/, $line);
		my $thisLine;
		foreach (@codonOrder){
			if ($removeTetra{$_}){
				$removed++;
				next;
			}
			$thisLine.=$_."\t";
			$cols++;
		}
		$thisLine=~ s/\t$/\n/;
		$secondPart.=$thisLine;
	}
	elsif($line=~ /^\d/){
		my $thisLine;
		my @frequencies=split(/\t/, $line);
		my $pos=-1;
		foreach my $freq(@frequencies){
			$pos++;
			next if ($removeTetra{$codonOrder[$pos]});
			$thisLine.=$freq."\t";
		}
		$thisLine=~ s/\t$/\n/;
		$secondPart.=$thisLine;
	}
	elsif($.==1){
		$firstLine=$line;
	}
}
close LRN;

print "Tetramers Removed:\t".$removed."\n";

open(OUT, ">".$out) || die $!;
print OUT $firstLine."\n";
print OUT "% ".$cols."\n";
my $thisLine.="% 9\t";
for (my $i=$cols; $i > 1; $i--){
	$thisLine.="1\t";
}
$thisLine=~ s/\t$/\n/;

print OUT $thisLine;
print OUT $secondPart;
close OUT;

exit 0;
