#!/usr/bin/perl

=head1 DESCRIPTION

map_project_names.pl -- Map IMG project names to your own. Creates Symbolic links with your project names to extracted IMG tar balls.

=head1 USAGE

perl map_project_names.pl

=head2 Options


	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

Sunit Jain, (Mon Feb 23 12:55:40 EST 2015)
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
my $mapFile="names_map.txt";
my $path2col1="./";
my $outPath="./";

GetOptions(
        'm|map:s'=>\$mapFile,
        'p1|path2col1:s'=>\$path2col1,
        'o|outdir:s'=>\$outPath,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

my %projects;
my $MAP=FileHandle->new();
open( $MAP, "<", $mapFile) || die $!;
while(my $line=<$MAP>){
    chomp $line;
    next unless ($line=~/^\d+/);
    next unless $line;
    my($imgName, $sampleNum, $desc)=split(/\t/, $line);
    my $oldPath=File::Spec->catdir($path2col1,$imgName);
    my $newPath=File::Spec->catdir($outPath,$sampleNum);
    symlink($oldPath, $newPath);
}
close $MAP;


