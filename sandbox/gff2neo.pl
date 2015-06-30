#!/usr/bin/perl

=head1 DESCRIPTION

gff2neo.pl -- Read a GFF(version 3) file and create nodes and relationship files for upload to a GraphDB. Tested on Neo4j v2.2.3.

=head1 USAGE

perl gff2neo.pl -gff myGFFv3_file.gff

=head2 Options

	-gff		<CHAR>		GFF version 3 file.
	-prefix	-p	<CHAR>		prefix for outputs
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

Sunit Jain, (Tue Jun 30 08:16:35 EDT 2015)
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

use Bio::Tools::GFF;

my $help;
my $version=fileparse($0)."\tv0.0.1b";
my ($gffFile, $prefix);
GetOptions(
	'gff:s'=>\$gffFile,
	'p|prefix:s'=>\$prefix,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

die "[FATAL] A 'GFF v3' File required to run. See '".fileparse($0)." -h' for help." if (! -s $gffFile);
if(!$prefix){
	$prefix=fileparse($gffFile, ".gff");
}

my $parser = new Bio::Tools::GFF->new(-file=> $gffFile, -gff_version => 3);

my (%parents,%outputs);
while( my $result = $parser->next_feature ) {
	my ($id,@junk)= $result->get_tag_values("ID");
	my $type = $result->primary_tag();

	if(!$result){
		last;
	}

	my $seq_id = $result->seq_id();
	my $source = $result->source_tag();
	my $strand = $result->strand();
print $strand. " | ";
next;
	$strand =~ s/-1/-/g;
	$strand =~ s/1/+/g;
	my $end = $result->end();
	my $start = $result->start();
	my @atts = $result->get_all_tags();

	# Open a file handle for each 'type' of feature
	my($NODES, $REL);
	unless($outputs{$type}{"NODES"}){
		$outputs{$type}{"REL"}=FileHandle->new();
		$outputs{$type}{"NODES"}=FileHandle->new();
		$outputs{$type}{"ATTRS"}=\@atts;

		# NODES File
		$NODES = $outputs{$type}{"NODES"};
		open($NODES, ">>", $prefix."_".$type.".nodes") || die $!;
		my $node_header="Scaffold\tSource\tType\tStart\tEnd\tStrand\t";
		$node_header.= $_."\t" foreach (@atts);
		$node_header=~ s/\t$/\n/;
		print $NODES $node_header;

		# RELATIONS File
		$REL = $outputs{$type}{"REL"};
		open($REL, ">>", $prefix."_".$type.".rel") || die $!;
	}
	
	# Redirect output to corresponding file handle
	$NODES=$outputs{$type}{"NODES"};	
	$REL=$outputs{$type}{"REL"};
	
	# Write to Nodes file
	my $line="$seq_id\t$source\t$type\t$start\t$end\t$strand\t";
	foreach my $tag(@{$outputs{$type}{"ATTRS"}}){
		my($value, @junk)=$result->get_tag_values($tag);
		$line.=$value."\t";
	}
	$line=~ s/\t$/\n/;
	print $NODES $line;

	# Write to Relationship file

	next;

	if($type eq "mRNA"){
		my ($parent,@junk)= $result->get_tag_values("Parent");
		$parents{$id} = $parent;
	}
	if($type eq "exon"){
		#find out transcript (parent) and gene for THIS exon
		my ($parent,@junk)= $result->get_tag_values("Parent");
		my $transcript = $parent;
		my $gene = $parents{$transcript};	
		print "$seq_id\t$source\t$type\t$start\t$end\t.\t$strand\t.\tgene_id \"$gene\";transcript_id \"$transcript\";\n";
	}
}

# Close all open file handles.
foreach my $type(keys %outputs){
	close $outputs{$type}{"NODES"};
	close $outputs{$type}{"REL"}
}

exit 0;
