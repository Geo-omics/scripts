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

my $help;
my $version=fileparse($0)."\tv0.0.2b";
my ($gffFile, $prefix);
my $minLen = 200;
GetOptions(
	'gff:s'=>\$gffFile,
	'p|prefix:s'=>\$prefix,
	'l|len:i'=>\$minLen,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

die "[FATAL] A 'GFF v3' File required to run. See '".fileparse($0)." -h' for help." if (! -s $gffFile);
if(!$prefix){
	$prefix=fileparse($gffFile, ".gff");
}

# Can't use BioPerl since it can't deal with the mess that is IMG. Doesn't parse strands from IMG gff outputs.
# use Bio::Tools::GFF;
# my $parser = new Bio::Tools::GFF->new(-file=> $gffFile, -gff_version => 3);

my (%parents,%outputs);
my $GFF=FileHandle->new();
open( $GFF, "<", $gffFile) || die $!;
while(my $line=<$GFF>){
	chomp $line;
	next unless $line;
	next if($line=~ /^#/);

	my($seq_id, $source, $type, $start,$end,$score,$strand,$phase,$attributes)=split(/\t/, $line);
	my @attribs=split(/;/, $attributes);
	my (%atts, $locus_tag);
	foreach my $item(@attribs){
		my($key, $value)=split(/\=/,$item);
		$atts{$key}=$value;
		$locus_tag=$value if (lc($key) eq "locus_tag");
	}
	my @a=sort keys %atts;

	$strand=~ s/-1/-/g;
	$strand=~ s/1/+/g;
	my $length=abs($start-$end)+1;

	next unless ($length >= $minLen);

	# Open a file handle for each 'type' of feature
	my($NODES, $REL);
	unless($outputs{$type}{"NODES"}){
		$outputs{$type}{"REL"}=FileHandle->new();
		$outputs{$type}{"NODES"}=FileHandle->new();
		$outputs{$type}{"ATTRS"}=\@a;

		# NODES File
		$NODES = $outputs{$type}{"NODES"};
		open($NODES, ">>", $prefix."_".$type.".nodes") || die $!;
		my $node_header="Scaffold\tSource\tType\tStart\tEnd\tLength\tStrand\t";
		$node_header.= $_."\t" foreach (@{$outputs{$type}{"ATTRS"}});
		$node_header=~ s/\t$/\n/;
		print $NODES $node_header;

		if($locus_tag){
			# RELATIONS File
			$REL = $outputs{$type}{"REL"};
			open($REL, ">>", $prefix."_".$type.".rel") || die $!;
			print $REL "ID\tTHIS\tTO\tTHAT\n";
		}
	}
	
	# Redirect output to corresponding file handle
	$NODES=$outputs{$type}{"NODES"};	
	$REL=$outputs{$type}{"REL"};
	
	# Write to Nodes file
	my $node;
	foreach my $tag(@{$outputs{$type}{"ATTRS"}}){
		my $value= $atts{$tag};
		$node.=$value."\t";
	}
	$node.="$seq_id\t$source\t$type\t$start\t$end\t$length\t$strand\t";
	$node=~ s/\t$/\n/;
	print $NODES $node;

	# Write to Relationship file
	my $rel=$locus_tag."__".$seq_id."\t".$locus_tag."\tis_".$type."_on\t".$seq_id."\n";
	print $REL $rel;
}
close $GFF;

# Close all open file handles.
foreach my $type(keys %outputs){
	close $outputs{$type}{"NODES"};
	close $outputs{$type}{"REL"}
}

exit 0;
