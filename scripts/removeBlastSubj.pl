#!/user/bin/perl

use strict;

my $in=$ARGV[0];
my $blastOut= $ARGV[1];
my $out= $$.".QueriesFromListRemoved.out";
my $list= $$.".QueriesFromList.out";

my %exclude;
open(LIST, $in)|| die $!;
while (my $line=<LIST>){
	next if ($line=~ m/^#/);
	chomp ($line);
	next unless ($line);
	$line=~ s/ //g;
	$line=~ s/\r//g;
	$line=lc($line);
	$exclude{$line}++;
}
close LIST;

print keys(%exclude)."\n";

open(BOUT, $blastOut) || die $!;
open(OUT, ">".$out);
open(OUT2, ">".$list);
my $count=0;
while (my $line= <BOUT>){
	next if ($line=~ m/^#/);
	chomp ($line);
	next unless ($line);

	my ($query, $subj, @etc)=split(/\t/, $line);
	chomp($query, $subj);
	$subj=~ s/ //g;

	$subj=lc($subj);
	if ($exclude{$subj}){
		$count++;
		print OUT2 $line."\n";
	}
	else{
		print OUT $line."\n";
	}
}
print "Matches Found:".$count."\n";
close BOUT;
close OUT;
