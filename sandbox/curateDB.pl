#!/usr/bin/perl
use strict;
use Getopt::Long;

my $in;
my $out = $$.".fasta";
my $delTmp;
my $hasGItags;
my $printDuplicates;

GetOptions(
	'i|in:s'=> \$in,
	'o|out:s'=> \$out,
	't|del'=> \$delTmp,
	'g|gi'=>\$hasGItags,
	'd|duplicates'=>\$printDuplicates,
	'h|help'=> sub{	system('perldoc', $0); exit; },
);

if (! $in){	system('perldoc', $0); exit; }

open (IN, $in)||die "[ERROR] $in: $!\n";
my $tmpFile= $$.".tmp";
open (TMPW, ">".$tmpFile);
my $u=0;
my $s=0;
my $t=0;
my %seen;

while(my $line=<IN>){
	chomp $line;
	next unless $line;
	if ($line=~ m/^>/){
		$t++;
		$line=~ s/\>/\+/g;
		$line=~ s/^\+/>/;
		print TMPW $line."\n";
	}
	else{
		print TMPW $line."\n";
	}
}
close IN;
close TMPW;
print "#Number of Headers before Curation: $t\n";

open (TMPR, $tmpFile) || die "[ERROR] $tmpFile: $!\n";
open (OUT, ">".$out);
my $printSeq= 0;
$/= ">";
while (my $line=<TMPR>){
	chomp $line;
	next unless $line;

	if ($hasGItags){
		my ($giTag, $giNum, @etc)=split(/\|/, $line);
		print OUT ">".$line unless ($seen{$giNum});
		$seen{$giNum}++;
	}
	else{
		my($name, @seqs)=split(/\n/,$line);
		print OUT ">".$line."\n" unless $seen{$name};
		$seen{$name}++;
	}
}
$/= "\n";
close TMPR;
close OUT;
print "#Number of Sequences after Curation: ".keys(%seen)."\n";
if ($printDuplicates){
	
	foreach my $k(keys %seen){
		#print ".";
		print $k."\t".$seen{$k}."\n" if ($seen{$k} >1);
	}
}
if (! $delTmp){
	print "#Deleting temporary file...\n";
	unlink $tmpFile;
}
exit;
