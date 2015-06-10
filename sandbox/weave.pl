#!/user/bin/perl
use strict;
use Getopt::Long;
use File::Basename;

my $ext="tally";
my $out;
GetOptions(
	'e:s'=>\$ext,
	'o:s'=>\$out,
);

my @DBs;
my @files=glob("*.".$ext);
open(OUT, ">".$out);
print OUT "#Transcripts\t";
my %master;
print @files." Files will be tallied...!\n";
foreach my $f(@files){
	my $dbName=basename($f,"\.$ext"); #split(/\_/, $f);
	push(@DBs, $dbName);
	print OUT $dbName."\t";
	my $fh;
	open($fh, $f) || die "[error] $f: $! \n";
	while (my $line=<$fh>){
		next if ($line=~ m/^#/);
		chomp $line;
		$line=~ s/\r//g;
		next unless $line;

		my @cols=split(/\t/, $line);
		$master{$cols[0]}{$dbName}=$cols[1];
	}
	close $fh;
}
print OUT "DB-presence\n";

foreach my $key(keys %master){
	print OUT $key."\t";
	my $total=0;
	foreach my $db(@DBs){
		my $v;
		if($master{$key}{$db}){$v = $master{$key}{$db}}
		else{$v=0}
		print OUT $v."\t";
		$total++ if($v != 0);
	}
	print OUT $total."\n";
}
