#! usr/bin/perl -w
chomp($ARGV[0]);
$f=$ARGV[0];
@out=split(/\./, $f);
$o=$out[0].".jn.fasta";
open (OUT, ">".$o);
open(IN, "$f") || die "ERROR:".$!;

while(<IN>){
	$l = $_;
	chomp($l);
	if (/^>(.+?)\s/) {
		print OUT ">$1\n";
		}
		else {
		print OUT "$l\n";
		
	
		}
		}

