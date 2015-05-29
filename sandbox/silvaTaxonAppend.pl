#!/usr/bin/perl

=head1 DESCRIPTION

	silvaTaxonAppend.pl: get the SSU database and make the taxon file.

=head1 USAGE

	# For Blast outpts
	perl silvaTaxonAppend.pl -blast <blast output> -out <output filename>
	OR
	# For Mapper files
	perl silvaTaxonAppend.pl -mapper <mapper output> -out <output filename>

=head2 Options

	-db	[characters]	Silva database fasta; really any fasta that has header in the format ">accession_number[SPACE]description"
						Default: "ssu119" with location set as: /omics/PublicDB/silva/release_119/SILVA_119_SSURef_tax_silva.fasta (SSU version 119)
	-blast	[characters]	Blast output performed against the fasta file.
	-mapper	[characters]	mapper.pl output performed on the blast output.
	-out	[characters]	output file
	
	-h or -help	[boolean]	help; This page.
	-v or -version	[boolean]	script version.

=head1 Author

	Sunit Jain, (Thu Aug  1 15:40:59 EDT 2013)
	sunitj [AT] umich [DOT] edu

=cut

use strict;
use Getopt::Long;

my ($out, $isBlast, $isMapped);
my $fasta;
my %db=(
	"ssu111"=>"/geomicro/data1/COMMON/publicDB/silva/release_111/SSURef_111_NR_tax_silva.fasta",
	"ssu115"=>"/geomicro/data1/COMMON/publicDB/silva/release_115/SSURef_NR99_tax_silva.fasta",
	"ssu119"=>"/omics/PublicDB/silva/release_119/SILVA_119_SSURef_tax_silva.fasta"
);

my $version="silvaTaxonAppend.pl\tv0.1.1";
GetOptions(
	'db=s'=>\$fasta,
	'blast:s'=>\$isBlast,
	'mapper:s'=>\$isMapped,
	'o|out=s'=>\$out,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system('perldoc', $0); exit;},
);

print "#\t".$version."\n";

die "[ERROR $0] Blast or mapped file required\n" if(! $isBlast && !$isMapped);

if (! $fasta){
	$fasta=$db{"ssu119"};
}
elsif($db{lc($fasta)}){
	$fasta=$db{lc($fasta)};
}

if(-s $fasta){
	$fasta=$fasta;
}
else{
	die "[FATAL] Invalid fasta file: $fasta\n";
}

my %dictionary;
$/=">";
open(FASTA, "<".$fasta)|| die $!;
while(my $line=<FASTA>){
	chomp $line;
	next unless $line;

	my($header, @sequence)=split(/\n/, $line);
	my($accNum, @desciption)=split(/\s+/, $header);

	my $seq=join("", @sequence);
	my $desc=join(" ", @desciption);

	$dictionary{$accNum}=$desc;
}
$/="\n";

$isBlast ? &parseBlast : &parseMapper;

sub parseBlast{
	open(BLAST, "<".$isBlast)||die $!;
	open(OUT, ">".$out)|| die $!;
	while(my $line=<BLAST>){
		next if $line=~ /^#/;
		chomp $line;
		next unless $line;

		my ($query, $subject, @restOfData)=split(/\t/, $line);
		my $data=join("\t", @restOfData);

		print OUT $query."\t".$subject."\t".$dictionary{$subject}."\t".$data."\n";
	}
	close BLAST;
	close OUT;
}

sub parseMapper{
	open(MAP, "<".$isMapped)||die $!;
	open(OUT, ">".$out)||die $!;
	while(my $line=<MAP>){
		chomp $line;
		if($.==6){
			my($subject, $query, @restOfData)=split(/\t/, $line);
			my $data=join("\t", @restOfData);

			print OUT $subject."\tDescription\t".$query."\t".$data."\n";
		}
		elsif($line=~ /^#/){
			print OUT $line."\n";
		}
		else{
			my($subject, $query, @restOfData)=split(/\t/, $line);
			my $data=join("\t", @restOfData);

			print OUT $subject."\t".$dictionary{$subject}."\t".$query."\t".$data."\n";
		}
	}
	close OUT;
	close MAP;
}
