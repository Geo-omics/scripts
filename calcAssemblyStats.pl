#! /usr/bin/perl

use strict;

my $fasta=$ARGV[0];

## Read Fasta File and compute length ###
my $length;
my $totalLength; 
my $totalContigs;
my @arr;
open(FASTA, $fasta)|| die $!;
$/=">";
while(my $line=<FASTA>){
	chomp $line; 
	next unless $line;

	my ($header, @sequence)=split(/\n/, $line);
	my $length=length(join("", @sequence));

	push (@arr, $length);
	$totalLength += $length; 
}
$/="\n";
close(FASTA);

my @sort = sort {$b <=> $a} @arr;
my $n50; 
my $L50;
foreach my $val(@sort){
	$n50+=$val;
	$L50++;
	if($n50 >= $totalLength/2){
		print "N50:\t$val\n";
		last; 
	}
}

__END__

my $in=$ARGV[0];
my (@path)=split("\/", $in);
my $file=$path[-1];
my %lenCount;
my $count=0;

open(IN, $in) || die $!."\n";
open(OUT, ">len_".$file) || die $!."\n";

$/=">";
## Read Fasta File and compute length ###
my $totalLength;
my $totalContigs;
my $lenMax=0; 
while(my $line=<IN>){
	chomp $line; 
	next unless $line;

	my ($header, @sequence)=split(/\n/, $line);
	my $length=length(join("", @sequence));

	$lenCount{$length}++;
	$lenMax= $lenMax < $length ? $length : $lenMax;
	$totalLength+=$length;
	$totalContigs++;
}
close(IN);
$/="\n";

my @sort = sort {$a <=> $b} keys %lenCount;
my @Lr;
foreach my $k(@sort){
	my $v=$lenCount{$k};
	for (my $i==0; $i< $v; $i++){
		push(@Lr, $k);
#		print $k." ";
	}
#	print "\n";
}

my $median= int((scalar(@Lr)/2) + 0.5); #l - 1;

print "N50:\t$Lr[$median]\n";
print "Number of Bases: $totalLength\n";
print "Number of Contigs: $totalContigs\n";
print "Longest Contig:$lenMax\n";
