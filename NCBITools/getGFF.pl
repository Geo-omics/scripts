#!/usr/bin/perl

=head1 DESCRIPTION

	getGFF.pl -- Given a list of contig names extract GFF data.

=head1 USAGE

	perl getGFF.pl -list contig_names.list -gff annotated_metagenome.gff

=head2 Options

	-list 	<CHAR>	list of contigs
	-gff	<CHAR>	metagenome GFF file.
	-col	<INT>	Column number that contains the contig names; start count from 1. [ default = 1]

	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

	Sunit Jain, (Thu Jan  2 12:41:53 EST 2014)
	sunitj [AT] umich [DOT] edu

=head1 License

	This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

	This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;
use File::Basename;

my $help;
my $version="getGFF.pl\tv0.1.0";
my ($list, $gff);
my $col=1;
GetOptions(
	'l|list:s'=>\$list,
	'g|gff:s'=>\$gff,
	'c|col:i'=>\$col,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

if ($col==0){ warn "Column 0, does not compute! Start your counts from 1\nAssuming you meant Column 1\n"; $col=1}

print $list;
open(LIST, "<".$list)|| die $!;
my %index;
while(my $line=<LIST>){
	chomp $line;
	$index{uc($line)}++;
}
close LIST;

print "\t..";

$col--;
my $out=fileparse($list, ".list");
open(GFF, "<".$gff)|| die $!;
open(OUT, ">".$out.".gff")|| die $!;
while(my $line=<GFF>){
	chomp $line;
	my (@data)=split(/\t/, $line);
	my $contig=$data[$col];
	print OUT $line."\n" if ($index{uc($contig)});
}
close GFF;
close OUT;
print ".Done.\n"

