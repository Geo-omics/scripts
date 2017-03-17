#!/usr/bin/perl

=head1 DESCRIPTION

	nameClassFiles.pl -- Rename class files according to a tab-separated list of old/original names to new/more sensible names.

=head1 USAGE

	perl nameClassFiles.pl -tsv old_and_new_filenames.tsv -ext fasta

=head2 Options

	-tsv	<CHAR>	col1=old name; <TAB> col2=new name
	-out	<CHAR>	Output Folder [default: "Renamed"]
	-ext	<CHAR>	extensions for the old and new names if the tsv doesn't already have them.
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

	Sunit Jain, (Mon Feb 10 14:26:12 EST 2014)
	sunitj [AT] umich [DOT] edu

=head1 License

	This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

	This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;
use File::Spec;
use File::Copy "cp";

my ($tsv);
my $out="Renamed";
my $ext="fasta";
my $help;
my $version="nameClassFiles.pl\tv0.0.2b";
GetOptions(
	'list|tsv:s'=>\$tsv,
	'out:s'=>\$out,
	'ext:s'=>\$ext,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

unless (-e $out){mkdir($out, 0755)};

my %tracker;
open(TSV, "<".$tsv)|| die $!;
while(my $line=<TSV>){
	next if ($line=~ "^#");
	chomp $line;
	next unless $line;
	
	my($old, $new)=split(/\t/, $line);
	if ($ext){
		$old.=".".$ext;
		$new.=".".$ext;
	}
	$tracker{$old}=$new;
}
close TSV;

foreach my $file(keys %tracker){
	my $new=File::Spec->catfile( $out, $tracker{$file} );
	print "Creating:\t$new\n";
	cp($file, $new);
}
