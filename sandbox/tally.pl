#!/user/bin/perl
use strict;
use Getopt::Long;

my $in;
my $master;
my $out=$$.".out";
my $bs=0;
my $printValues; # which column?
my $fasta;
GetOptions(
	'i:s'=>\$in,
	'm:s'=>\$master,
	'o:s'=>\$out,
	's:f'=>\$bs,
	'f|fasta:s'=>\$fasta,
	'v|values|value:i'=>\$printValues,
);

my %seqLen;
if($fasta){
	$/=">";
	open(FASTA, "<".$fasta)|| die $!;
	while(my $line=<FASTA>){
		chomp $line;
		next unless $line;
		my($header, @s)=split(/\n/, $line);
		my $seq=join("",@s);
		$header=~ s/^>//;
		$seqLen{$header}=length($seq);
	}
	close(FASTA);
	$/="\n";
}

$printValues--;
open(IN, $in)|| die "[Error] $in: $!\n";
my %seen;
while (my $line=<IN>){
	next if ($line=~ m/^#/);
	chomp $line;
	next unless $line;
	$line=~ s/\r//g;

	my @lineParts=split(/\t/, $line);

	if(($fasta) && ($seqLen{$lineParts[0]})){
		# Alternative to previous condition; don't need it if query 100% identical with same start and stop positions.
		my($subjStart,$subjStop)=sort{$a <=> $b} ($lineParts[8],$lineParts[9]);
		next if(($lineParts[6]==$lineParts[8]) && ($lineParts[7]==$lineParts[9]) && ($lineParts[2]==100) && ($subjStop==$seqLen{$lineParts[0]}));
	}
	else{
		next if ($lineParts[0] eq $lineParts[1]); # Don't want query and subj to be the same
	}
	
	next if($seen{$lineParts[0]}); # Only need the top hit.
	next if ($lineParts[-1] < $bs); # Don't need anything with a bitscore less than user provided BS.

	if(! $printValues){
		$seen{$lineParts[0]}++;
	}
	else{
		$seen{$lineParts[0]}=$lineParts[$printValues];
	}
}
close IN;

open(MASTER, $master) || die "[Error] $master: $!\n";
open(OUT, ">".$out);
while (my $line=<MASTER>){
	next if ($line=~ m/^#/);
	chomp $line;
	next unless $line;
	$line=~ s/\r//g;

	if ($seen{$line}){
		print OUT $line."\t".$seen{$line}."\n";
	}
	else{
		print OUT $line."\t0\n";
	}
}
close MASTER;
close OUT;
exit;
