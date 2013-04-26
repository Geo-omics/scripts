#AA:
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=17230130&rettype=fasta

#!/usr/bin/perl -w

 # USAGE: perl generic_ncbi_dwnld.pl <File containing List of GI Numbers> <Database> <Return File Type>
 #  <Database> : protein, nucleotide, etc 
 #  <Return File Type>: fasta, gb, etc
 # refer eutils website for other database and filetype options.

# or you could use Batch Enterez.

my %g2gi;
open(G2GI, "ARGV[0]")||die "ERROR: $!";
while (my $line2=<G2GI>){
	my($geneName, $giNumber)= split(/\s/, $line2);
	chomp($geneName);
	chomp($giNimber);
	$g2gi{$geneName}=$giNumber;
}
close G2GI;

my %geneLinks;
foreach my $gName(keys %geneKey){
	my $url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$ARGV[1]&id='.$g2gi{$gName}.'&rettype=$ARGV[2]' if defined ($g2gi{$gName});
	$geneLinks{$gName}=$url if defined ($g2gi{$gName});
}
my $time=int(keys(%geneLinks))/3
print "The script will take atleast ".$time." seconds\n";

use LWP::Simple;
open (LEFT, ">notfound.txt");
my $tmpCogs="allCogs.tmp";
open (COGS, ">>".$tmpCogs);
my $count;
while (my($gene, $link)=each(%geneLinks)){
	print ".";
	my $content =get($link);
	die "Couldn't get $gene\:$link" unless defined $content;
	
	if ($content){
		print COGS $content."\n";
	}
	else{
		print LEFT "$gene\n";
	}
	sleep 0.35;
}
print "\n";
close LEFT;
close COGS;

open (COGED, "$tmpCogs") || die "ERROR: $!";
open (OUT, ">allCogs.faa");
while (my $line3=<COGED>){
	if ($line3=~ m/^\>/){
		chomp($line3);
		$line3=~ s/\>//g;
		print OUT ">".$line3."\n";
	}
	else{
		print OUT $line3;
	}
}
close OUT;
close COGED;
unlink ($tmpCogs);
