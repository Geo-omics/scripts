#!/usr/bin/perl

=head1 DESCRIPTION

	Given a ribopicker tabular output; get distribution by taxa.

=head1 USAGE

	perl ribopicker_taxonDist.pl -tsv ribopicker_output.tsv -out output_filename.txt

=head2 Options

	-tsv	[characters]	.tsv file created by ribopicker along with the fasta files.
	-taxon	[characters]	path to the taxon file, the default shoud work for all cases; default = "/geomicro/data1/COMMON/ribopickerDB/taxons/rrnadb.taxon"
	-out	[characters]	output file name.
	
	-v	[boolean]	script version
	-h	[boolean]	this screen

=head1 Author

	Sunit Jain, (Tue Jul 23 10:40:45 EDT 2013)
	sunitj [AT] umich [DOT] edu

=cut

use strict;
use Getopt::Long;

my $tsv;
my $taxon="/geomicro/data1/COMMON/ribopickerDB/taxons/rrnadb.taxon";
my $out;
my $version="0.0.3";
GetOptions(
	'tsv=s'=>\$tsv,
	'taxon=s'=>\$taxon,
	'out:s'=>\$out,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system('perldoc', $0); exit;},
);

my %refCount;
open(TSV, "<".$tsv)|| die $!;
while(my $line=<TSV>){
	next if $line=~ /^#/;
	chomp $line;
	next unless $line;
	
	my($query_id,$ref_id,$ref_start,$ref_bp_aligned,$query_coverage,$query_identity)=split(/\t/, $line);
	$refCount{$ref_id}++;
}
close TSV;

my %taxaDist;
open(TAXON, "<".$taxon)||die $!;
while(my $line=<TAXON>){
	next if $line=~ /^#/;
	chomp $line;
	next unless $line;
	
	my($ref_id, $ref_length, $domain, $phylum, $db_acc_num)=split(/\t/, $line); # ssr108_1        1437    Bacteria        Thermotogae     A61579,U37021,X91822
	next unless $refCount{$ref_id};
	
	$phylum=lc($phylum);
	$taxaDist{$phylum}+=$refCount{$ref_id};
}
close TAXON;

open(OUT, ">".$out)||die $!;
print OUT "# Taxa\tDistribution\n";
foreach my $phylum(keys %taxaDist){
	print OUT $phylum."\t".$taxaDist{$phylum}."\n";
}
close OUT;
