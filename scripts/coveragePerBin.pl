#!/usr/bin/perl

=head1 DESCRIPTION

coveragePerBin.pl -- Do this.

=head1 USAGE

perl coveragePerBin.pl

=head2 Options

    -conf   <FILE>      Concatenated confidence file. # make sure you've used the loyalty scores in getClassFasta.pl.
                        # This script won't handle scaffold mapping to multiple bins.
    -cov    <FILE>      Sample Scaffold file.
    -prefix <TEXT>      Prefix that should be used for scaffold names in the `-cov` file.
    -out    <TEXT>      Output file name.
    -min    <INT>       Only count scaffolds with cov >= min; default=0
    
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

Sunit Jain, (Fri Jul 17 13:19:32 EDT 2015)
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

my $minCov=0;
my ($confFile,$covFile,$prefix,$outFile);
my $help;
my $version=fileparse($0)."\tv0.0.1";
GetOptions(
    'conf:s'=>\$confFile,
    'cov:s'=>\$covFile,
    'prefix:s'=>\$prefix,
    'out:s'=>\$outFile,
    'min:i'=>\$minCov,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

my %binHash;
my %errorTrigger;
my $CONF=FileHandle->new();
open( $CONF, "<", $confFile) || die $!;
while(my $line=<$CONF>){
    next if ($line=~ /^#/);
    chomp $line;
    next unless $line;
    
    my($bin,$scaffold,$conf)=split(/\t/, $line);
    
    if (($binHash{$scaffold}) && ($binHash{$scaffold} != $bin)){
        $errorTrigger{$scaffold}{$binHash{$scaffold}}++;
        $errorTrigger{$scaffold}{$bin}++
    }
    
    $binHash{$scaffold}=$bin;
}
close $CONF;

if (scalar(keys %errorTrigger) > 0) {
    my $numOfChar=80;
    my $error=("-" x $numOfChar)."\n";
    $error.="[ERROR] Scaffold assigned to multiple bins!\n";
    $error.="Can't decide which bin to add it to.\n";
    $error.="Use the `getClassFasta.pl' script with a loyalty value >= 51 (recommend: 60)\n";
    $error.="Here are the offending scaffolds and their bins:\n";
    foreach my $scaffold(keys %errorTrigger){
        $error.="$scaffold";
        foreach my $bin(keys %{$errorTrigger{$scaffold}}){
            $error.="\t$bin";
        }
        $error.="\n";
    }
    $error.=("-" x $numOfChar)."\n";
    die $error;
}



my %binCov;
my $COV=FileHandle->new();
open( $COV, "<", $covFile) || die $!;
while(my $line=<$COV>){
    next if ($line=~ /^#/);
    chomp $line;
    next unless $line;
    
    my($scaffoldName,$avgCov,$len)=split(/\t/, $line);
    next if ($avgCov < $minCov);
    
    my $scaffold=$prefix.$scaffoldName;
    next unless $binHash{$scaffold};
    
    my $bin= $binHash{$scaffold};
    $binCov{$bin}{"COV"}+=$avgCov;
    $binCov{$bin}{"COUNT"}++;
}
close $COV;

my $OUT=FileHandle->new();
open( $OUT, ">", $outFile) || die $!;
print $OUT "#Bin\tCoverage\n";
foreach my $bin(keys %binCov){
    print "No scaffolds found for Bin: $bin\n" if ($binCov{$bin}{"COUNT"} == 0);
    print $OUT $bin."\t".($binCov{$bin}{"COV"}/$binCov{$bin}{"COUNT"})."\n";
}
close $OUT;

exit 0;