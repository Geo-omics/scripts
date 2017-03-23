#!/user/bin/perl

=head1 Motivation

	Get a subset of query list from the mapper script

=head1 USAGE
	
	perl mapper_getQueryList.pl -log <log file> -list <list of Queries of interest> -out <output>
	
=head1 Author

	Sunit Jain, July 2013
	
=cut

use strict; 
use Getopt::Long;

my ($logFile,$list);
my $version="0.3.1";
my $compatible="0.3.0 +";
my $out=$$.".list";

GetOptions(
	'log=s'=>\$logFile,
	'list=s'=>\$list,
	'o|out:s'=>\$out,
	'v|version'=>sub{print $version."\n"."Compatible with mapper script version $compatible"; exit;},
	'h|help'=>sub{system('perldoc', $0); exit;},
);

my %index;
open(LIST, "<".$list)|| die $!;
while(my $line=<LIST>){
	next if ($line=~ /^#/);
	chomp($line);
	$line=~ s/\r//;
	next unless $line;
	
	$index{$line}++;
}
close LIST;

open (LOG, "<".$logFile)|| die $!;
open (OUT, ">".$out)|| die $!;
while(my $line=<LOG>){
	next if ($line=~ /^#/);
	chomp($line);
	$line=~ s/\r//;
	next unless $line;
	
	my ($subj, @queries)=split(/\t/, $line);
	next unless $index{$subj};

	print OUT "#".$subj."\n";
	foreach my $q(@queries){
		print OUT $q."\n";
	}
}
close LOG;
close OUT;
