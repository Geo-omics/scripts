#!/user/bin/perl
use strict;
use Getopt::Long;

my $ext="tally";
my $out;
GetOptions(
	'e:s'=>\$ext,
	'o:s'=>\$out,
);

my @files=glob("*.".$ext);
open(OUT, ">".$out);
print OUT "#Transcripts\t";
my %master;
print @files." Files will be tallied...!\n";
foreach my $f(@files){
	my ($dbName, $etc)=split(/\_/, $f);
	print OUT $dbName."\t";
	my $fh;
	open($fh, $f) || die "[error] $f: $! \n";
	while (my $line=<$fh>){
		next if ($line=~ m/^#/);
		chomp $line;
		$line=~ s/\r//g;
		next unless $line;

		my @cols=split(/\t/, $line);
		push(@{$master{$cols[0]}}, $cols[1]);
	}
	close $fh;
}
print OUT "DB-presence\n";

foreach my $key(keys %master){
	print OUT $key."\t";
	my $total=0;
	foreach my $v(@{$master{$key}}){
		print OUT $v."\t";
		$total++ if ($v>0);
	}
	print OUT $total."\n";
}

