#!/usr/bin/perl

=head2 USAGE

	perl topX.pl -b <blast output> -t <top X hits; default: 5> -o <output_File_Name>
	
=cut

use strict;
use Getopt::Long;

my $in;
my $x=5;
my $out=$$.".top5";
GetOptions(
	'b:s'=>\$in,
	't:i'=>\$x,
	'o|out:s'=>\$out,
	'h'=> sub{system("perldoc $0 \| cat"); exit;},
);

open (TB, $in) || die "[err] $in: $! \n";
my %file;
while (my $line=<TB>){
	next if ($line=~ m/^#/);
	chomp $line;
	$line=~ s/\r//;
	next unless $line;
	
	my($query, @etc)=split(/\t/, $line);
	push(@{$file{$query}}, $line);
}
close TB;

my($nothing, $fileName)=split(/\_/, $in);
open (OUT, ">".$out)|| die "[err] $out: $! \n";;
while(my($k,$val)= each(%file)){
	my $count=0;
	foreach my $v(@{$val}){
		if ($count < $x){
			print OUT $v."\n";
		}
		$count++;
	}
}
close OUT;
