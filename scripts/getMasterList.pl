#!/user/bin/perl
use strict;
use Getopt::Long;

my $ext="out";
my $out=$$.".list";
my $col=1;
my $bs=0;
GetOptions(
	'e:s'=>\$ext,
	'o:s'=>\$out,
	'c:i'=>\$col,
	's:f'=>\$bs,
);

my @listOfFiles=glob("*.".$ext);
print @listOfFiles." Filenames provided\n";

my $c= $col-1;
my %masterList;
open (OUT, ">".$out);
foreach my $f(@listOfFiles){
	my $fh;	
	open($fh, $f) || die "[error] $f: $!\n";
	while (my $line=<$fh>){
		next if ($line=~ m/^#/);
		chomp $line;
		$line=~ s/\r//g;
		next unless $line;

		my @cols=split(/\t/, $line);
		print OUT $cols[$c]."\n" unless ($masterList{$cols[$c]});
		$masterList{$cols[$c]}++;
	}
	close $fh;
}
close OUT;
exit;
