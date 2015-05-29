#!/usr/bin/perl

=head1 DESCRIPTION

coveragePerScaffold.pl -- Using the GenomeCoverageBed default output calculate the coverage per scaffold and the whole genome.

=head1 USAGE

perl coveragePerScaffold.pl -bed genomeCovBed.txt

=head2 Options

        -bed    -b      <CHAR>      Default output format from the genomeCoverageBed command in bedtools.
        -min            <INT>       Only count bases with coverage more than or equal to this; default=0;
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

Sunit Jain, (Thu Feb 26 15:24:40 EST 2015)
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
my $version=fileparse($0)."\tv0.0.2";
my $bedFile;
my $minCov=0;
GetOptions(
        'b|bed:s'=>\$bedFile,
        'min:i'=>\$minCov,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

my %cov;
my $FILE=FileHandle->new();
open( $FILE, "<", $bedFile) || die $!;
while(my $line=<$FILE>){
    chomp $line;
    next unless $line;
    my($scaffold,$depth,$numBases, $scafLen,$fraction)=split(/\t/, $line);
    next if ($depth!~ /\d+/);
    next if ($depth < $minCov);

    $cov{$scaffold}{"LEN"}=$scafLen;
    $cov{$scaffold}{"CMCOV"}+=($depth*$numBases);
}
close $FILE;

my $genomeCov=0;
foreach my $scaf(keys %cov){
    my $coverage=$cov{$scaf}{"CMCOV"}/$cov{$scaf}{"LEN"};
    if ($scaf eq "genome") {
        $genomeCov=$coverage;
        next;
    }
    
    print $scaf."\t".$coverage."\t".$cov{$scaf}{"LEN"}."\n";
}
print "Genome\t".$genomeCov."\t".$cov{"genome"}{"LEN"}."\n";

exit 0;

