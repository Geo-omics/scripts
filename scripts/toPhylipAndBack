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

=head1 Outputs

	*.aka file, is your alias file. Guard it with your life! I won't be able to convert your fasta file to their original names without this.
	*.toPhylip.fasta OR *.toOriginal.fasta file, will have the same name as the '.aka' file. This is your converted fasta file.

=head1 Questions/Suggestions/Feedback/Beer

	Sunit Jain, October 2012 (sunitj [AT] umich [DOT] edu)

=cut

use strict;
use Getopt::Long;

my ($file, $aliasFile, $out, $toPhylip, $toOriginal);

GetOptions(
	'f|fasta=s'=>\$file,
	'aka|alias=s'=>\$aliasFile,
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
	$out=$aliasFile.".toPhylip.fasta";
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
	my $i=0;
	foreach my $l(keys %fasta){
		$i++;
		my $alias=sprintf("%010d", $i);
		print ALIAS $alias."\t".$l."\n";
		print OUT ">".$alias."\n".$fasta{$l}."\n";
	}
	print "\t$i new aliases created...\n";
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


