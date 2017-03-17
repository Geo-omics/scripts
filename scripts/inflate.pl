#!/usr/local/bin/perl

=head1 DESCRIPTION

inflate.pl -- Use this script if you've used dereplicate.pl to remove duplicate reads and wish to know what the original number of reads would have been had you used the original dataset.

=head1 USAGE

perl inflate.pl -b file.blast -c file.clust -o revized.blast

=head2 Options

	-blast	-b	<CHAR>	blast output
	-clust	-c	<CHAR>	dereplicate output clust file.
	-out	-o	<CHAR>	output file; will contain an extra column called "multiplier" denoting the size of the cluster of identical reads that the query represents.
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

Sunit Jain, (Fri Oct  9 15:00:04 EDT 2015)
sunitj [AT] umich [DOT] edu

=head1 License

This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;
use FileHandle;
use File::Basename;

my $blastFile= "DOES_NOT_EXIST.blast";
my $clustFile= "DOES_NOT_EXIST.clust";
my $outFile= $$.".blast.output";
my $help;
my $version=fileparse($0)."\tv0.0.1b";
GetOptions(
	'b|blast:s'=>\$blastFile,
	'c|clust:s'=>\$clustFile,
	'o|out:s'=>\$outFile,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

die "[FATAL] Problem with required option '-b' or '-blast'. See '-h' for more.\n" if(! -e $blastFile);
die "[FATAL] Problem with required option '-c' or '-clust'. See '-h' for more.\n" if(! -e $clustFile);

# ClusterNumber  Size    Representative  SeqHeaders
# c1      1       HWI-D00709:31:C6UHRANXX:8:1202:18201:19009 1:N:0:GATCAGCG       
# c2      2       HWI-D00709:32:C6RYLANXX:4:1313:11459:92771 1:N:0:GATCAGCG       HWI-D00709:32:C6RYLANXX:1:2207:3019:66349 1:N:0:GATCAGCG 

my %clust;
my $CLUST=FileHandle->new();
open( $CLUST, "<", $clustFile) || die $!;
while(my $line=<$CLUST>){
	chomp $line;
	next unless $line;
	next if ($line=~ /^#/);

	my ($name, $size, $rep, @members)=split(/\t/, $line);
	my ($head, $strand)=split(/ /, $rep);	
	$clust{$head}=$size;
}
close $CLUST;

my $BLAST=FileHandle->new();
my $OUT=FileHandle->new();
open( $BLAST, "<", $blastFile) || die $!;
open( $OUT, ">", $outFile) || die $!;
while(my $line=<$BLAST>){
	chomp $line;
	next unless $line;
	next if ($line=~ /^#/);

	my($query, @etc)=split(/\t/, $line);
	print $OUT $line."\t".$clust{$query}."\n";
}
close $BLAST;
close $OUT;

exit 0;
