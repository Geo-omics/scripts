#!/usr/bin/perl -w

=head1 MOTIVATION

	In a dereplicated file genererated by my 'dereplicate.pl' script, each sequences may be a representative of a cluster of sequences.
	If you've used the dereplicated file for your Blasts, you may wish to find out the actual number of sequences that 'would have' gotten
	the same result had the file not been dereplicated. By using this script you can get to that number	by simply multiplying this
	result to your blast output.

	Why use a dereplicated file at all, you say?
	One of the many reasons is that it cuts down on the blast times by a HUGE margin!
	for more on dereplication and its affects, goto: 

=head2 USAGE

	derepClustMap.pl -clust <name of the .clust file> -l <log file> -m <original mapped file>

=head3 OPTIONAL

	-o Output File Name

=head2 AUTHOR

	Sunit Jain, November, 2011
	sunitj [AT] umich [DOT] edu

=cut


use strict;
use Getopt::Long;

my $log;
my $mapped;
my $clustF;
my $out=$$.".derep.map";
my $version= "derep_ClusterMap.pl v 0.1.9";
GetOptions(
	'l|log=s'=>\$log,
	'm|mapped=s'=>\$mapped,
	'clust=s'=>\$clustF,
	'o|out:s'=>\$out,
	'h'=>sub{system('perldoc', $0); exit;},
	'v|version'=>sub{print "# $version\n"; exit;}
);
print "# $version\n";

my @clustFiles=split(/\,/, $clustF);

open (LOG, $log) || die "[ERROR] $log:$!\n";
my (%index, %map);
while(my $line=<LOG>){
	next if ($line=~ /^#/);
	chomp($line);
	$line=~ s/\r//;
	next unless $line;

	my($query, @hits)=split(/\t/, $line);
	foreach my $h(@hits){
		$index{$h}++;
		push (@{$map{$query}},$h);
	}
}
close LOG;

my $numNames=scalar(keys %index);
print "[".$0."] Looking for [ ".$numNames." ] items.\n";


my $totalSeqs=0;
my %multiplier;
my $expandedNumSeqs=0;
foreach my $clust(@clustFiles){
	open (CLUST, $clust) || die "[ERROR] $clust:$!\n";
	print "Parsing: $clust\n";
	while(my $line=<CLUST>){
		next if ($line=~ /^#/);
		chomp($line);
		$line=~ s/\r//;
		next unless $line;

		my($cNum, $size, $rep, @seqNames)=split(/\t/, $line);
		my ($name, $strand)=split(/\s/, $rep);
		$name=~ s/^@//;
		if ($index{$name}){
			$multiplier{$name}=$size;
			$expandedNumSeqs+=$size;
			delete $index{$name};
		}
		$totalSeqs+=$size;
	}
	close CLUST;
}

my $notFound=scalar(keys %index);

print "# [".$0."] ".$numNames." dereplicated sequences correspond to ".$expandedNumSeqs." real sequences.\n";
print "# [".$0."] ".$notFound." not found!\n" unless $notFound==0;

open (OUT, ">".$out) || die "$mapped : $!";
open(MAPPED, $mapped) || die "$mapped : $!";
while(my $line=<MAPPED>){
	chomp $line;
	next unless $line;


	if ($line=~ /^#/){
		if($.==4){
			print OUT $line."\t\t\tRead Counts:\t".$expandedNumSeqs."\n";
		}
		else{
			print OUT $line."\n";
		}
 		next;
	}
	
	my($query, @data)=split(/\t/, $line);
	my $newTotal=0;
	if ($map{$query}){
		foreach my $h(@{$map{$query}}){
			$newTotal+=$multiplier{$h} if ($multiplier{$h});
		}
		$data[-3]=$newTotal;
#		$data[-1]=($newTotal/$totalSeqs) * 100;
		$line=join("\t", $query, @data);
	}
	else{
		print STDERR "[".$0."]".$query." not found!\n";
	}
	print OUT $line."\n";
}
close MAPPED;
