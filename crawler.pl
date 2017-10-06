#!/usr/bin/perl

=head1 DESCRIPTION

crawler.pl	-	Checks the 'content.list' file for any new script that may have been created in sandbox. Copies the script to the appropriate directory. Maintains and updates these scripts and README for each category of scripts.

=head1 USAGE

	perl crawler.pl -list ~/src/sandbox/contents.list -path2scripts ~/Github/scripts -path2sandbox ~/src/sandbox

=head2 Options

	-list	-l	<FILE>	A specially formatted file to categorize all the scripts. See the "contents.list" file, available in the same repo, for the description.
	-path2scripts	-p	<PATH> Location of the folder that'll be committed/pushed to Github/BitBucket/Hg etc.
	-path2sandbox	-s	<PATH>	Location of the folder where all the scripts reside.
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message. press q to exit this screen.

=head1 Author

	Sunit Jain, (Fri Sep  6 13:48:16 EDT 2013)
	sunitj [AT] umich [DOT] edu

=cut

use strict;
use Getopt::Long;
use File::Spec;
use File::Path;
use File::Copy;
use FileHandle;

my $version="crawler.pl\tv0.0.6";
my $list="contents.list";
my $path= "./";
my $sandbox="./sandbox";

GetOptions(
	'p|path2scripts:s'=>\$path,
	's|path2sandbox:s'=>\$sandbox,
	'l|list:s'=>\$list,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system('perldoc', $0); exit;},
);
print "\# $version\n";

open(LIST, "<".$list) || die $!;
my $folder="./ERROR/";
mkpath($folder);
my %table_of_contents;
my $MD=FileHandle->new;
while(my $line=<LIST>){
	next if $line=~ /^#/;
	chomp $line;
	next unless $line;

	if ($line=~ /\:$/){
	# Folders
		$line=~ s/\:$//;
		$folder=File::Spec->catfile($path, $line);
		close $MD;
		unless($folder=~ /wrappers/g){
			my $contentsMD=File::Spec->catfile($folder, "README.md");
			$MD=FileHandle->new;
			open($MD, ">".$contentsMD) || die "Can't create table of contents at: $contentsMD\n";
			print $MD "## Script Descriptions\n";
		}

		if (! -d $folder){
			my $dir = eval{ mkpath($folder) };
			die "Failed to create $line: $@\n" unless $dir;
		}
	}
	else{
	# Files
		my($fileName, @comment)=split(/\t/, $line);
		my $file=File::Spec->catfile($sandbox, $fileName);
		unless($file=~ /\./g){
			warn "'$file' must be a folder. Skipping...\n";
			next;
		}

		if(! -e $file){
			warn "Not Found: $file\n";
			next;
		}
		copy($file, $folder) || die "Failed to copy $file: $!\n";

		unless($folder=~ /wrappers/g){
			$fileName=~s/\_/\\\_/g;
			print $MD "* **".$fileName."**"."\t".join("\t", @comment)."\n";
		}
	}
}
close LIST;
