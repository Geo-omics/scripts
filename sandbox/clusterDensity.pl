#!/usr/local/bin/perl

=head1 DESCRIPTION

clusterDensity.pl -- Do this.

=head1 USAGE

perl clusterDensity.pl

=head2 Options


	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

Sunit Jain, (Mon Mar  2 16:41:57 EST 2015)
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

my $help;
my $version=fileparse($0)."\tv0.0.1b";
my $clustFile="results.txt";
GetOptions(
        'c|clusters:s'=>\$clustFile,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

my $FILE=FileHandle->new();
open( $FILE, "<", $clustFile) || die $!;
while(my $line=<$FILE>){
    chomp $line;
    next unless $line;
    
    my($thresh, $fMeasure, $combinedClust)=split(/\t/, $line);
    my @clustComma=split(/\,/,$combinedClust);
    my @clustSemiColon=split(/\;/, $combinedClust);
    my $totalClusters=scalar(@clustSemiColon);
    my $totalNodes=scalar(@clustComma);
    print $thresh."\t".$totalNodes."\t".$totalClusters."\t".($totalNodes/$totalClusters)."\n";
}
close $FILE;


