#! /usr/bin/perl

=head1 USAGE

	perl toPhylipAndBack.pl -phylip -f InputFile.fasta -alias file_Prefix
	OR
	perl toPhylipAndBack.pl -original -f InputFile.fasta -alias file_Prefix

=head1 Description

=head2	-phylip

	Converts a regular fasta file to a 10-digit header fasta file and creates an 'aka' file that keeps track of the name changes.

=head2	-original

	Converts the fasta file generated from '-phylip' to it's original self using the 'aka' file.

=head2 Advanced Commands

	-prefix	-p	Add a uniform prefix before the aliased name (default="")
	-num_digits	The alias will be this many digits long. (default=10)
	-num_start Start the aliased numbering from x (default: 1)

=head1 Outputs

	*.aka file, is your alias file. Guard it with your life! I won't be able to convert your fasta file to their original names without this.
	*.aliased.fasta OR *.toOriginal.fasta file, will have the same name as the '.aka' file. This is your converted fasta file.

=head1 Questions/Suggestions/Feedback/Beer

	Sunit Jain, October 2012 (sunitj [AT] umich [DOT] edu)
	Last Updated: July 2016

=cut

use strict;
use Getopt::Long;

my ($file, $aliasFile, $out, $toPhylip, $toOriginal);
my $defPrefix="";
my $num_digits=10;
my $num_start=1;
GetOptions(
	'f|fasta=s'=>\$file,
	'aka|alias=s'=>\$aliasFile,
	'p|prefix:s'=>\$defPrefix,
	'num_digits:i'=>\$num_digits,
	'num_start:f'=>\$num_start,
	'phylip'=>\$toPhylip,
	'original'=>\$toOriginal,
	'h|help'=>sub{system('perldoc', $0); exit;},
);

if (! $file || ! $aliasFile){die "[ERROR: $0] Missing required input. \nFor help, type: '$0 -h'\n";}
if (! $toPhylip && ! $toOriginal){die "[ERROR: $0] Missing conversion type (-phylip ? OR -original ?). \nFor help, type: '$0 -h'\n";}

open(FASTA, $file) || die "[ERROR: $0] Unable to open Fasta File: $! \n";
my %fasta;

$/=">";
while(my $line=<FASTA>){
	chomp $line;
	next unless $line;

	my($name, @seqs)=split("\n", $line);
	my $seq=join("", @seqs);

	$seq=~ s/\n//g;

	$fasta{$name}=$seq;
}
$/="\n";

my $aka;
if($toPhylip){
	$aka=$aliasFile;
	$aliasFile=~ s/\.aka$//;
	$out=$aliasFile.".aliased.fasta";
	$aka=$aliasFile.".aka";
	&createAlias;
}
elsif($toOriginal){
	$aka=$aliasFile;
	$aliasFile=~ s/\.aka$//;
	$out=$aliasFile.".toOriginal.fasta";
	&back2original;
}

sub createAlias{
	open(ALIAS, ">".$aka) || die "[ERROR: $0] Unable to create aliases: $! \n";
	open(OUT, ">".$out) || die "[ERROR: $0] Unable to create aliased from original Fasta: $! \n";

	print ALIAS "# ALIAS\tORIGINAL\n";
	my $i=$num_start-1;
	foreach my $l(keys %fasta){
		$i++;
		my $alias=sprintf("%0${num_digits}d", $i);
		my $defLine= $defPrefix.$alias;
		print ALIAS $defLine."\t".$l."\n";
		print OUT ">".$defLine."\n".$fasta{$l}."\n";
	}
	print "\t".$i-($num_start-1)." new aliases created...\n";
}

sub back2original{
	open(AKA, $aka)|| die "[ERROR: $0] Unable to open alias File: $! \n";
	open(OUT, ">".$out) || die "[ERROR: $0] Unable to create original from aliased Fasta: $! \n";

	my $i=0;
	while(my $line=<AKA>){
		next if $line=~ /^#/;
		chomp $line;
		next unless $line;

		my($alt, $orig)=split(/\t/, $line);
		$i++;
		print OUT ">".$orig."\n".$fasta{$alt}."\n" if ($fasta{$alt});
	}
	print "\t$i new aliases created...\n";
}
