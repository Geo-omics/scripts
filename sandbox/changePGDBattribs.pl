#!/usr/bin/perl

=head1 DESCRIPTION

	changePGDBattribs.pl -- Change the attributes in "data/ptools-local/pgdbs/user/<NAME>/1.0/input/organism.dat" to avoid confusion while opening a PGDB.

=head1 USAGE

	perl changePGDBattribs.pl -list project_names.list

=head2 Options

	-list	<CHAR>	List of the fasta files you used to run MetaPathways.
	-data	<CHAR>	Location of the pathway-tools data folder; default:/opt/packages/pathway-tools/data/
	-prefix	<CHAR>	If you wish to add a prefix to the names, add it here (no spaces allowed)
	-bkp	<BOOLEAN>	Keep the old version (just in case);
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

	Sunit Jain, (Thu Feb 27 15:52:47 EST 2014)
	sunitj [AT] umich [DOT] edu

=head1 License

	This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

	This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;
use File::Basename;
use File::Spec;
use File::Copy;

my $help;
my $version="changePGDBattribs.pl\tv0.0.2b";
my ($list, $keepBkp, $prefix);
my $dataLoc="/opt/packages/pathway-tools/data/";

GetOptions(
	'l|list:s'=>\$list,
	'd|data:s'=>\$dataLoc,
	'p|prefix:s'=>\$prefix,
	'bkp|backup'=>\$keepBkp,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

my @ext=qw(.fasta .fa .fna .faa);
open(LIST, "<".$list)|| die $!;
while(my $line=<LIST>){
	next if ($line=~ /^#/);
	chomp $line;
	next unless $line;

	my $projectName=fileparse($line, @ext);
	my $orgFile=File::Spec->catfile($dataLoc, "ptools-local/pgdbs/user/", $projectName."cyc", "1.0/input/organism.dat");

	if (-e $orgFile){
		print $orgFile."\n";
		&changeAttributes($orgFile);
	}
}
close LIST;

sub changeAttributes{
	my $file=shift;
	my $bkp= $file.".bkp";

	copy($file, $bkp) || die "Could not create backup for $file:".$!."\n";
	
	my $out=$file;
	my ($FH, $OFH);
	open($FH, "<".$bkp)||die $!;
	open($OFH, ">".$out) || die $!;
	my $useThisName;
	while(my $line=<$FH>){
		if($line=~ /^ID/){
			my($tag, $value)=split(/\t/, $line);
			$useThisName=($prefix ? $prefix."_" : "").$value;
			print $OFH $line;
		}
		elsif($line=~ /^NAME/){
			print $OFH "NAME\t$useThisName";
		}
		elsif($line=~ /^ABBREV-NAME/){
			print $OFH "ABBREV-NAME\t$useThisName";
		}
		else{
			print $OFH $line;
		}
	}
	close $FH;
	close $OFH;
	unlink $bkp unless $keepBkp;
}

