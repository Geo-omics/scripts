#! usr/bin/perl
chomp($ARGV[0]);
$f=$ARGV[0];
@out=split(/\./, $f);
$o=$out[0]."_justnames.txt";
open (OUT, ">".$o);
open(IN, "$f") || die "ERROR:".$!;

my $oldQuery = "";
my $oldRef = "";
my $oldScore = -1;

while(<IN>){
	chomp($_);
	my @row = split /\s/, $_;

	 if (/^>(.+)/) {
		print OUT "@row[0]\t@row[1]\n";
			
	}
	next;
	
}
close(OUT);
close(IN);