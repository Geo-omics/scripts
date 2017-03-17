#! usr/bin/perl

=head1 DESCRIPTION

	fileChopper.pl - chops a file into a defined number of portions.

=head1 USAGE

	perl fileChopper.pl [-f Filename]

=head2 Optional

	-p OR -parts:	Split the given file  into these many parts.
	-l OR Lines:	Maximum number of lines allowed in a file.
	-avgsize:	defines the average size of a  file. Only used when you don't specify '-p' OR '-l'.

=head1 AUTHOR

	Sunit Jain, Aug, 2011.
	sunitj [AT] umich [DOT] edu.
	 
=cut

use strict;
use Getopt::Long;

my $file;
my $parts;
my $numLines; # max number of Lines allowed in a file
my $avgSize= 5000;

GetOptions(
	'f|file:s'=>\$file,
	'p|parts:i'=>\$parts,
	'l|lines:i'=>\$numLines,
	'avgsize:i'=>\$avgSize,
	'h|help'=>sub{system('perldoc', $0); exit;},
);

my ($totalLines, $filePartSize);
&checkArgs;
print "#Lines per file: \~ $numLines\n";
print "#File chopped into:$parts\n";
my $lastFileName="";
&chopFile;
&checkLastFile;

sub checkArgs{
	if (! $file){
		system('perldoc', $0); exit;
	}

	$totalLines= `wc -l $file`;
	$filePartSize=$totalLines/$avgSize;
	if ((! $parts || $parts == 0 ) && (! $numLines || $numLines == 0)){
		($filePartSize >= 1) ? &bigFile : &smallFile;
	}
	elsif ($parts > 0 && $numLines > 0){
		warn "#[WARNING] You may only use one of the two options (Number of Parts OR Max Line) at a time.\n";
		warn "#[WARNING] Only using the 'Number of Parts' cut-off\n";
		&calcNumLines;
	}
	elsif ($parts > 0 && (! $numLines || $numLines == 0)){
		&calcNumLines;		
	}
	elsif ((! $parts || $parts == 0) && $numLines > 0){
		return;
	}

}

# Get number of Lines to chop the file into; if the user doesn't mention it.
sub bigFile{
	$parts= int($filePartSize + 0.5);
	warn "#[Warning]This file seems like it's larger than my definition of an 'average file'.\n";
	print "#Since you didn't mention a preference I'll split the File  into \~ $parts parts\n";
	print "#You may change the 'average file size' by using the '-avgsize' flag.\nCurrently it is set at $avgSize\n";
	$numLines= int(($totalLines/$parts) + 0.5);
}

sub smallFile{
	my $p= int($filePartSize + 0.5);
	$parts= ($p == 0) ? 1 : $p;
	warn "#[Warning]This file seems like it's smaller than my definition of an 'average file'.\n";
	print "#Since you didn't mention a preference I'll split the File  into \~ $parts parts\n";
	print "#You may change the 'average file size' by using the '-avgsize' flag.\nCurrently it is set at $avgSize\n";
	$numLines= int(($totalLines/$parts) + 0.5);
}

sub calcNumLines{
	$numLines= int(($totalLines/$parts) + 0.5);
}

sub chopFile{

	open (FILE, $file) || die "[error] $! : $file\n";
	my $numOfLines=0;
	my $partCount=1;
	my $totalLines=0;
	my ($fName, @ext)= split(/\./, $file);
	my $restOfName=join(".", @ext);
	my $out=$fName.".0".$partCount.".".$restOfName;
	my $fh;
	open ($fh, ">".$out);
	my %seen=();
	print "\n#List of Files:\n";
	print $out."\n";
	while (my $b = <FILE>) {
		chomp $b;
		next unless $b;
		$numOfLines++;
		$totalLines++;
		if ($numOfLines < $numLines){
			print $fh $b."\n";
		}
		elsif ($numOfLines == $numLines){
			print $fh $b."\n";
			close $fh;

			$numOfLines=0;
			$partCount++;
			my $out=$fName.".0".$partCount.".".$restOfName;
			$lastFileName= $out;
			open ($fh, ">".$out);
			print $out."\n";
		}
	}
}

sub checkLastFile{
	unlink $lastFileName if (-z $lastFileName);
}
