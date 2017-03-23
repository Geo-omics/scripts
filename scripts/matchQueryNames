#!/usr/bin/perl

=head1 DESCRIPTION

matchQueryNames.pl -- Do this.

=head1 USAGE

perl matchQueryNames.pl

=head2 Options


	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

Sunit Jain, (Fri Jul 18 14:46:20 EDT 2014)
sunitj [AT] umich [DOT] edu

=head1 License

This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;

my ($meta1, $meta2, $out);
my $help;
my $version="matchQueryNames.pl\tv0.0.1b";
GetOptions(
	'1|meta1:s'=>\$meta1,
	'2|meta2:s'=>\$meta2,
	'o|out:s'=>\$out,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";
my %metaDataIndex;
print $meta1."\n";
open(META1, "<".$meta1)|| die $!;
while(my $line=<META1>){
	next if ($line=~/^#/);
	chomp $line;
	next unless $line;
	match("1",$line);
}
close META1;

open(META2, "<".$meta2)|| die $!;
while(my $line=<META2>){
	next if ($line=~/^#/);
	chomp $line;
	next unless $line;
	match("2",$line);
}
close META2;

open(OUT, ">".$out)||die $!;
print OUT $meta1."\t".$meta2."\n";
my $gt_2=0;
my %seen;
foreach my $meta(keys %metaDataIndex){
	if (@{$metaDataIndex{$meta}}>2){
		$gt_2++;
	}
	elsif(@{$metaDataIndex{$meta}}==2){
		my $line;
		foreach my $name(@{$metaDataIndex{$meta}}){
			$line.=$name."\t";
		}
		$line=~ s/\t$/\n/;
		next if $seen{$line};
		print OUT $line;
		$seen{$line}++;
	}
}
close OUT;
print "gt_2 = ".$gt_2."\n";

sub match{
	my $prefix=shift;
	my $line=shift;

	my($alias, @metadata)=split(/\t/,$line);
	foreach my $m(@metadata){
		push(@{$metaDataIndex{$m}},$prefix."_".$alias);
	}
}

