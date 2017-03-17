#!/usr/local/bin/perl -w

=head1 DESCRIPTION

genomicFluidity.pl -- Do this.
# aai -- http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1236649/
# Genomic Fluidity -- http://www.biomedcentral.com/1471-2164/12/32 [(conserved1+conserved2)/(Total1+Total2)] * 100

=head1 USAGE

perl genomicFluidity.pl

=head2 Options

	-ref	<CHAR>	Reference Genome
	-blast	<CHAR>	Path to the Post-Blast directory; USE: blastp -outfmt "6 std qcovs"
	-genomes <CHAR>	Path to the protein fasta files
	-ext	<CHAR>	extension for the blast outputs; (Default="blastp")
	-out	<CHAR>	3 column tab delimited output file. GeneName <TAB> ClusterNum <TAB> <Genome>
	-map	<CHAR>	

	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

Sunit Jain, (Wed Aug 20 12:08:14 EDT 2014)
sunitj [AT] umich [DOT] edu

=head1 License

This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;
use File::Basename;
use File::Spec;

my $help;
my $version="genomicFluidity.pl\tv0.2.0";
my ($ref, $postBlast, $map, $genomeDir);
my $setCov=70;
my $setPerc=30;
my $setScore=0;
my $setEval=1e-2;
my $ext="blastp";
my $DEBUG;
my $cluster="tmp.cluster";
GetOptions(
	'r|ref:s'=>\$ref,
	'b|blast:s'=>\$postBlast,
	'g|genomes:s'=>\$genomeDir,
	'c|cov|coverage:i'=>\$setCov,
	'id|perc_id:i'=>\$setPerc,
	's|bit_score:i'=>\$setScore,
	'eval|evalue:s'=>\$setEval,
	'e|ext:s'=>\$ext,
	'o|out:s'=>\$cluster,
	'm|map:s'=>\$map,
	'test'=>\$DEBUG,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

my %index;
if($map){
open(MAP, "<".$map) || die $!;
while(my $line=<MAP>){
	chomp $line;
	next unless $line;
	next if ($line=~ /^#/);

	my($name, $num)=split(/\t/, $line);
	$index{$num}=$name;
}
close MAP;
}

my @files=<$postBlast/*.$ext>;
my @shortlist;
my $refName="";
if($ref){
	$refName=fileparse($ref, (".fasta", ".faa"));
	foreach my $f(@files){
		my $fName=fileparse($f, ".".$ext);
		my ($first, $last)=split(/_vs_/, $fName);
		if(($first == $refName) || ($last == $refName)){
			push(@shortlist, $f);
		}
	}
}
else{
	@shortlist=@files;
}
undef @files;

## SANITY CHECK ##
if($DEBUG){
	print $_."\n" foreach(@shortlist);
}
## ##

my $tmp=$$.".concat.blastp";
my $sortTmp=$$.".concat.sorted.blastp";
my (%aggregates, %geneIndex);
open(TMP, ">".$tmp) || die $!;
foreach my $file(@shortlist){
	my $FH;
	open($FH, "<".$file) || die $!;
	my %seen;
	my $fName=fileparse($file, ".".$ext);
	my ($first, $last)=split(/_vs_/, $fName);
	chomp($first, $last);
	die "[FATAL] Aweful Filename: $fName\n" if (($first eq "") || ($last eq ""));
	while(my $line=<$FH>){
		chomp $line;
		next unless $line;
		next if ($line=~ /^#/);

		my(@splitLine)=split(/\t/, $line);
		die "[FATAL] Blast file [ $file ] has a format that wasn't expected. See 'help'.\n" if (scalar(@splitLine) != 13);
		my $query=$splitLine[0];
		my $subj=$splitLine[1];
		my $pid=$splitLine[2];
		my $eval=$splitLine[10];
		my $score=$splitLine[11];
		my $cov=$splitLine[12];

		next if (($pid < $setPerc) || ($cov < $setCov) || ($score < $setScore) || ($eval > $setEval));
		next if $seen{$query}; # make sure we only get the top hit {for each query}.
		$seen{$query}++;
		$geneIndex{$query}=$first;
		print TMP $line."\n";		

		$aggregates{$first}{$last}{"PID"}+=$pid;
		$aggregates{$first}{$last}{"Count"}++;
	}
	close $FH;
}
close(TMP);

# Order by Increasing Evalue > Decreasing BitScore > Decreasing Coverage
# `sort -k11,11g -k12,12nr -k13,13nr $tmp > $sortTmp`;

# Order by Decreasing BitScore > Increasing Evalue > Decreasing Coverage
`sort -k12,12nr -k11,11g -k13,13nr $tmp > $sortTmp`;
unlink($tmp);

# Generate Gene Families
open(CAT, "<".$sortTmp)|| die $!;
my $clustNum=0;
my (%protClust, %refClust);
while(my $line=<CAT>){
	chomp $line;
	next unless $line;

	my($query, $subj, @etc)=split(/\t/, $line);
	if($protClust{$query}){
		next;
	}
	elsif($protClust{$subj}){
		$protClust{$query}=$protClust{$subj};
	}
	else{
		$clustNum++;
		$protClust{$query}=$clustNum;
		if($geneIndex{$query} eq $refName){
			$refClust{$clustNum}++;
		}
	}
}
close CAT;

open(CLUST, ">".$cluster)|| die $!;
foreach my $query(keys %protClust){
	print CLUST $query."\tc_".$protClust{$query}."\tf_".$geneIndex{$query}."\n";
}
close CLUST;

# Calculate AAI
my (%seen, %numProt);
print "QUERY\tSUBJ\tAAI\tGF\n";
foreach my $first(keys(%aggregates)){
	foreach my $last(keys(%{$aggregates{$first}})){
		next if $seen{$first}{$last};
		next if $seen{$last}{$first};

		my $aai;
		$aai= ($aggregates{$first}{$last}{"PID"} + $aggregates{$last}{$first}{"PID"})/($aggregates{$first}{$last}{"Count"} + $aggregates{$last}{$first}{"Count"});
		$aai=sprintf("%.3f",$aai);
		
		my $queryTotalProt=&getNumProt($first);
		my $subjTotalProt=&getNumProt($last);

		chomp($queryTotalProt, $subjTotalProt);

		my $uniqQuery=$queryTotalProt - $aggregates{$first}{$last}{"Count"};
		my $uniqSubj=$subjTotalProt - $aggregates{$last}{$first}{"Count"};

		my $GF = (($uniqQuery+$uniqSubj)/($queryTotalProt + $subjTotalProt));
		$GF=sprintf("%.3f",$GF);

		if($map){
			print $index{$first}."\t".$index{$last}."\t".$aai."\t".$GF."\n";
		}
		else{
			print "genome_".$first."\tgenome_".$last."\t".$aai."\t".$GF."\n";
		}
		$seen{$first}{$last}++;
		$seen{$last}{$first}++;
	}
}
unlink($sortTmp);
sub getNumProt{
	my $fileName=shift;
	my $file=File::Spec->catfile($genomeDir, $fileName.".fasta");
	if(! $numProt{$file}){
		$numProt{$file}=`grep -c "^>" $file`;
	}
	return $numProt{$file};
}




