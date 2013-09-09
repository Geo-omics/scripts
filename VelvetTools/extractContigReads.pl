#!/usr/bin/perl -w
##
##       extractContigReads.pl
##
##       Copyright 2009 Daniel Zerbino <zerbino@ebi.ac.uk>
##
##       This program is free software; you can redistribute it and/or modify
##       it under the terms of the GNU General Public License as published by
##       the Free Software Foundation; either version 2 of the License, or
##       (at your option) any later version.
##
##       This program is distributed in the hope that it will be useful,
##       but WITHOUT ANY WARRANTY; without even the implied warranty of
##       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##       GNU General Public License for more details.
##
##       You should have received a copy of the GNU General Public License
##       along with this program; if not, write to the Free Software
##       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
##       MA 02110-1301, USA.
##
##       A script to split out a single contig's reads from an LastGraph file
##       produced by the Velvet assembler.  The output is in fasta format.
##
##       Usage:  ./extractContigReads.pl <contig number> <directory> 
##
##
##       Where:  <contig number> is the number of the contig of interest.
##		 <directory> is the Velvet directory 
##
######################################################################
=head1 Usage:

	perl extractContigReads.pl -l list_of_Original_contig_Names.txt -dir assembly_directory -o output -r list of reads

	-summary	when you only wish to see the # of reads in each contig and NOT the names.
	-metav	if you've used metavelvet on your data.

=cut


use strict;
use Getopt::Long;

my $version="0.0.3b";
my $readList;
my $metav;
my $trans;
my $list;
my $directory;
my $summary;
my $OUT=$$."cov";

print "Extract Contig Reads	v$version\n";
my $commandLine = join(" ", $0, @ARGV);
print "You executed the following command:\n";
print $commandLine."\n";

GetOptions(
	'l|list:s'=>\$list,
	'r|reads:s'=>\$readList,
	'metav'=>\$metav,
	'trans'=>\$trans,
	'dir=s'=>\$directory,
	'o|out:s'=>\$OUT,
	'summ|summary'=>\$summary,
	'h|help'=>sub{system('perldoc', $0); exit;},
);

die "[$0:ERROR]Folder Name Required\n" if (!$directory);

my ($graphfile,$seqfile,$stats);
if ($metav){
	$graphfile = "$directory/meta-velvetg.LastGraph";
	$seqfile = "$directory/Sequences"; 
}
else{
	$graphfile = "$directory/LastGraph";
	$seqfile = "$directory/Sequences"; 
}

unless(-e $directory){die "$directory does not exist\n"};
unless(-e $graphfile){die "$graphfile does not exist, please re-run Velvet with the reads tracking on\n"};

my @nodes= `cut -d '_' -f 2 $list | sort -g -r`;

#pop @nodes;
chomp @nodes;
my %contigsList;
$contigsList{$_}++ foreach @nodes;
#print $_."\n" foreach keys(%contigsList);
#print "CONTIGS:".keys(%contigsList)."\n";
#exit;
open GRAPH, $graphfile;

$_ = <GRAPH>;
chomp $_;
my @data = split /\t/;
my %reads = ();
my $read_no;
my $readingSeqPath = 0;
my $recordReadIDs = 0;

my $contig_no;


# LastGraph file has Contig Number listed in descending order, if the contig number entered is greater than the very first Contig Number in the LastGraph file, it does not exist.
unless ($data[0] >= $nodes[0]) {die "Contig $nodes[0] does not exist in the LastGraph file.\n"};
$seqfile = "$directory/Sequences"; 
unless (-e $seqfile) {die "$seqfile does not exist, exiting.\n"};
undef @nodes;


while(<GRAPH>) {
	if (/SEQ/) {
		chomp;
		@data = split /\t/;
		$read_no = $data[1];
		$readingSeqPath = 1;
		next;
	}

	if (/NR/) {
		$readingSeqPath = 0;
	}

	if ($readingSeqPath == 1) {
		chomp;
		@data = split /\t/;
		if ($contigsList{abs($data[0])}) {
			$reads{$read_no} = abs($data[0]);
		}
	}

	if (/^NR\t-?(\d+)\t/) {
		$contig_no = $1;
		$recordReadIDs=1;
	} elsif (/^NR/) {
		$recordReadIDs = 0;
	}

	if ($recordReadIDs == 1) {
		chomp;
		@data = split /\t/;
		$read_no = $data[0];
		if ($contigsList{$contig_no}){
			$reads{$read_no} =$contig_no;
		}
	}
}
close GRAPH;
undef %contigsList;

open SEQ, $seqfile;
my %contigs;
while(<SEQ>) {
	if (/>/) {
		chomp;
	# the Sequences file is a fasta format file with a twist; the headers of each sequence have 2 additional columns such that col1 is the sequence name and col2 is a number assigned to the sequence. This number in column 2 is used to make graphs and can be searched for to get the reads used in a contig.
		@data = split /\t/;
		if ($reads{$data[1]}) {
			$data[0]=~ s/^\>//;
			$contigs{$reads{$data[1]}}{$data[0]}++;
#			print STDERR $reads{$data[1]}."\t".$data[0]."\n";
#			exit;
#			print $reads{$data[1]}."\t".$data[1]."\n";
		}
	}
}
undef %reads;
close SEQ;

my %LIST;
if (!$summary){open(LIST, ">".$readList)|| die "Couldn't write to: $readList. $!\n";}
open(OUT, ">".$OUT)|| die "Couldn't write to: $OUT. $!\n";
foreach my $c(keys %contigs){
	print OUT $c."\t".scalar(keys %{$contigs{$c}})."\t";
	if (!$summary){
		foreach my $r(keys %{$contigs{$c}}){
			print LIST $r."\n" unless $LIST{$r};
			print OUT $r."\t";
			$LIST{$r}++;
		}
	}
	print OUT "\n";
}
close OUT;
close LIST if (!$summary);
exit 0;


