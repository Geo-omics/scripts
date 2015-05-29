#!/usr/bin/perl

=head1 DESCRIPTION

	Given a ribopicker tabular output; get distribution by taxa.

=head1 USAGE

	perl ribopicker_taxonDist.pl -tsv ribopicker_output.tsv -out output_filename.txt

=head2 Options

	-tsv	[characters]	.tsv file created by ribopicker along with the fasta files.
	-taxon	[characters]	path to the taxon file, the default shoud work for all cases; default = 
			possible choices:
			'rrnadb'	=>	"/geomicro/data1/COMMON/ribopickerDB/taxons/rrnadb.taxon"
			'ssr'	=>	"/geomicro/data1/COMMON/ribopickerDB/taxons/ssr108.taxon"
			'slr'	=>	"/geomicro/data1/COMMON/ribopickerDB/taxons/slr108.taxon"
			OR provide the path to a custom taxon file.
	-out	[characters]	output file name.
	
	-v	[boolean]	script version
	-h	[boolean]	this screen

=head3 Taxon File format

	ribopicker_reference_id <TAB> reference_sequence_length <TAB> Domain <TAB> Phylum <TAB> original_db_accession_numbers (optional)
	EXAMPLE: ssr108_1        1437    Bacteria        Thermotogae     A61579,U37021,X91822
	
=head1 Author

	Sunit Jain, (Tue Jul 23 10:40:45 EDT 2013)
	sunitj [AT] umich [DOT] edu

=cut

use strict;
use Getopt::Long;

my $tsv;
my $taxon="ssr";
my $out;
my $version="0.0.3";
GetOptions(
	'tsv=s'=>\$tsv,
	'taxon:s'=>\$taxon,
	'out:s'=>\$out,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system('perldoc', $0); exit;},
);

if(lc($taxon) eq "ssr"){
	$taxon="/geomicro/data1/COMMON/ribopickerDB/taxons/ssr108.taxon";
}
elsif(lc($taxon) eq "slr"){
	$taxon="/geomicro/data1/COMMON/ribopickerDB/taxons/slr108.taxon";
}
elsif(lc($taxon) eq "rrnadb"){
	$taxon="/geomicro/data1/COMMON/ribopickerDB/taxons/rrnadb.taxon";
}
elsif((-e $taxon)&&(-s $taxon)) {
	print "$taxon set as the taxon db\n";
}
else{
	die "[ERROR]\tYou input '-taxon $taxon' was not valid.
	Please choose either 'ssr','slr' or 'rrnadb' depending on the database that you chose while running ribopicker.
	The script will now exit.\n";
}

my %refCount;
open(TSV, "<".$tsv)|| die "[ERROR]\tCould not open:$tsv:$!\n";
while(my $line=<TSV>){
	next if $line=~ /^#/;
	chomp $line;
	next unless $line;
	
	my($query_id,$ref_id,$ref_start,$ref_bp_aligned,$query_coverage,$query_identity)=split(/\t/, $line);
	$refCount{$ref_id}++;
}
close TSV;

my %taxaDist;
open(TAXON, "<".$taxon)||die "[ERROR]\tCould not open: $taxon : $!\n";
while(my $line=<TAXON>){
	next if $line=~ /^#/;
	chomp $line;
	next unless $line;
	
	my($ref_id, $ref_length, $domain, $phylum, $db_acc_nums)=split(/\t/, $line); # ssr108_1        1437    Bacteria        Thermotogae     A61579,U37021,X91822
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
