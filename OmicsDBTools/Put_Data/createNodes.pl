#!/usr/bin/perl

=head1 DESCRIPTION

createNodes.pl -- Do this.

=head1 USAGE

perl createNodes.pl

=head2 Options

	-config	-c	<CHAR>		Config file explicitly declaring which columns to make nodes out of and which the properties; Recommended but Optional; default = guess.
	-nodes	-n	<CHAR>		Read this nodes file created by one of the parsers and create nodes in the database; Required.
	-port	-p	<INT>		Port which the database is listening to; default: 7474.
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

Sunit Jain, (Tue Jun 30 08:14:43 EDT 2015)
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
my ($configFile, $nodesFile);
my $port = 7474;
GetOptions(
	'c|config:s'=>\$configFile,
	'n|nodes:n'=>\$nodesFile,
	'p|port:i'=>\$port,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

my $NODES=FileHandle->new();
open( $NODES, "<", $nodesFile) || die $!;
while(my $line=<$NODES>){

}
close $NODES;


