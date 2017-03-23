#!/usr/bin/perl

=head1 DESCRIPTION

	gff2tbl.pl	-	Convert JGI's 'GFF3' data to NCBI's ridiculous 'tbl' format.

=head1 USAGE

	perl gff2tbl.pl -fasta scaffolds.fasta -gff jgi_annotated.gff -gene jgi_annotated.gene_product.txt -tbl output.tbl

=head2 Options

	-fasta	[characters]	Original assembled Fasta file [Required]
	-gff	[characters]	JGI's GFF file [Required]
	-tbl	[characters]	Output tbl file [Required]

	-gene	[characters]	gene product file from JGI. Required if GFF doesn't have the tag 'product'

	-aka	[characters]	aliased file; from "toPhylipAndBack.pl" script
	-min	[integers]	minimum sequence length
	-prot_prefix	[characters]	Project specific prefix for CDS regions. example: "abc|unique_lab_name|"

	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message. press q to exit this screen.

=head1 Feedback/Bug Reports

Please use the Github issues page for any issues/bugs with these scripts. Make sure you add the name of the script in the issue header.

=head1 Author

	Sunit Jain, (Thu Oct 10 12:48:37 EDT 2013)
	Last Updated: June 02, 2016

=cut

use strict;
use Getopt::Long;

my ($fasta, $gff, $gene_product);
my ($tbl, $fasta_out);
my $minLen=200;
my $minGeneLen= 300; # just used an arbitary number, to reduce the amount of manual curation required afterwords. This is used to determine 'incompleteness' of a gene.
my $prot_prefix="";
my $aka;
my $help;
my $version="gff2tbl.pl\tv0.2.5";
GetOptions(
	'f|fasta:s'=>\$fasta,
	'gff:s'=>\$gff,
	'gene:s'=>\$gene_product,
	'aka:s'=>\$aka,
	'tbl:s'=>\$tbl,
	'min:i'=>\$minLen,
	'prot_prefix:s'=>\$prot_prefix,
	'o|out:s'=>\$fasta_out,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system('perldoc', $0); exit;},
);
print "\# $version\n";

die("[FATAL] Required parameters not found. See '-h' for help on how to use the script.\n") if ((! $fasta)||(! $gff)||(! $tbl));

my %gene_prod; # gene_prod{locusID - same as one in %annotation}= product name
my %otherInfo;
open(GP, "<".$gene_product); #|| die $!;
while(my $line=<GP>){
	chomp $line;
	next if $line=~ /^#/;
	next unless $line;

	my($img_gene_id, $locusID, $tag, @data)=split(/\t/, $line);
	if(lc($tag) eq "product_name"){
		$gene_prod{$locusID}=join(";", @data);
	}
	else{
		push(@{$otherInfo{$locusID}}, join(",",@data));
	}
}
close GP;

my %alias;
open(ALIAS, "<".$aka);
while(my $line=<ALIAS>){
	chomp $line;
	next if $line=~ /^#/;
	next unless $line;

	my($alt, $orig)=split(/\t/, $line);
	my($name, @desc)=split(/\s+/, $orig);
	$alias{$alt}=$name;
}
close ALIAS;

my %annotation;
open(GFF, "<".$gff)|| die $!;
while(my $line=<GFF>){
	chomp $line;
	next if $line=~ /^#/;
	next unless $line;

	&parseGFF3($line);
}
close GFF;

$/=">";
my %contigLen; # contig = length
open(FASTA, "<".$fasta) || die $!;
open(TBL, ">".$tbl) || die $!;
while(my $line=<FASTA>){
	chomp $line;
	next if $line=~ /^#/;
	next unless $line;

	my($header, @sequence)=split(/\n/, $line);
	my $seq=join("", @sequence);
	my ($name, @desc)=split(/\s+/, $header);
	my $parent;
	if($aka){
		$parent=$alias{$name};
	}
	else{
		$parent=$name;
	}

	next unless ($annotation{$parent});

	&find_Ns($seq, $name);

	my $len=length($seq);

	next if ($len < $minLen);
	print TBL ">Feature ".$name."\n"; #"\tLength:".$len."\n";
#	print ">Feature ".$name."\n"; #"\tLength:".$len."\n";
#	exit
	my $exonNum=0;
	foreach my $locusID(keys %{$annotation{$parent}}){
		my $original_contig_gene_start=$annotation{$parent}{$locusID}{"START"};
		my $original_contig_gene_stop=$annotation{$parent}{$locusID}{"STOP"};
		($original_contig_gene_start, $original_contig_gene_stop)=sort{$a<=>$b} ($original_contig_gene_start, $original_contig_gene_stop);

		# Question: Is the feature incomplete?
		my $incomplete_5="";
		my $incomplete_3="";
		# Workaround: Does the feature start at position 1 **AND** is less than "$minGeneLen" (300 by default). Assume it's incomplete.
		if (($original_contig_gene_start == 1) && (abs($original_contig_gene_stop - $original_contig_gene_start) <= $minGeneLen)){
			$incomplete_5="\<";
		}
		# Workaround: Does the feature stop at the end of the scaffold **AND** is less than "$minGeneLen" (300 by default). Assume it's incomplete.
		if(($original_contig_gene_stop == $len) && (abs($original_contig_gene_stop - $original_contig_gene_start) <= $minGeneLen)){
			$incomplete_3="\>";
		}

		my ($gene_start, $gene_stop);
		if($annotation{$parent}{$locusID}{"STRAND"}=~ /^\-/){
			($gene_start, $gene_stop)=($original_contig_gene_stop, $original_contig_gene_start);
		}
		else{
			($gene_start, $gene_stop)=($original_contig_gene_start, $original_contig_gene_stop);
		}

		print TBL $incomplete_5.$gene_start."\t";
		print TBL $incomplete_3.$gene_stop."\t";
		print TBL "gene\n";
		print TBL "\t\t\tlocus_tag\t$locusID\n";

		print TBL $incomplete_5.$gene_start."\t";
		print TBL $incomplete_3.$gene_stop."\t";
		print TBL $annotation{$parent}{$locusID}{"TYPE"}."\n";
		print TBL "\t\t\tlocus_tag\t$locusID\n";
		print TBL "\t\t\tprotein_id\t$prot_prefix"."$locusID\n";
		print TBL "\t\t\t";

		if($gene_prod{$locusID}){
			if($annotation{$parent}{$locusID}{"TYPE"}=~ /RNA/i){
				print TBL "product\t".$annotation{$parent}{$locusID}{"TYPE"}."-".$gene_prod{$locusID}."\n";
			}
			else{
				print TBL "product\t".$gene_prod{$locusID}."\n";
			}

			foreach my $info(@{$otherInfo{$locusID}}){
				print TBL "\t\t\tnote\t".$info."\n";
			}
		}
		elsif($annotation{$parent}{$locusID}{"TYPE"}=~ /RNA/i){
			print TBL "product\t".$annotation{$parent}{$locusID}{"TYPE"}."\n";
		}
		elsif($annotation{$parent}{$locusID}{"TYPE"} eq "exon"){
			$exonNum++;
			print TBL "number\t".$exonNum."\n";
		}
		elsif($annotation{$parent}{$locusID}{"TYPE"} eq "repeat_region"){
			my($locus, $rpt_type, $unit, $fam)=split(/__/, $locusID);
			print TBL "rpt_type\t".$rpt_type."\n";
			print TBL "\t\t\trpt_unit\t".$unit."\n";
			print TBL "\t\t\trpt_family\t".$fam."\n";
		}
		else{
			my  ($ID, @desc)=split(/\_/, $locusID);
			my $type=join("_", @desc);
			print TBL "note\thypothetical protein\n";
			print TBL "\t\t\tnote\tIMG_locus='".$locusID."'\n";
		}
	}
}
close FASTA;
close TBL;
$/="\n";

sub find_Ns{
	my $seq=shift;
	my $header=shift;
	while($seq=~ /N{10,5000}/ig){
		print STDERR $header."\t".$-[0]."\t".$+[0]."\t".($+[0]-$-[0])."\n";
	}
}

sub parseGFF3{
#http://gmod.org/wiki/GFF
# contig, source, type, start,stop,score,strand, phase,attributes
    my $line=shift;
    my ($contig, $source, $type, $start,$stop,$score,$strand, $phase,$attribs)=split(/\t/, $line);

    my(@attributes)=split(/\;/, $attribs);

    my ($locusID, $ID, $Name,$Alias, $Parent, $Target, $Gap, $Derives_from, $Note, $Dbxref, $Onto, $repeat_type, $repeat_unit, $repeat_fam, $product);
	my $repeatNumber=1;
    foreach my $att(@attributes){
		$locusID=$1 if ($att=~/locus_tag\=(.*)/);
		$ID= $1 if ($att=~/^ID\=(.*)/);
		$Name=$1 if ($att=~/^Name\=(.*)/);
		$Alias=$1 if ($att=~/^Alias\=(.*)/);
		$Parent=$1 if ($att=~/^Parent\=(.*)/);
		$Target=$1 if ($att=~/^Target\=(.*)/);
		$Gap=$1 if ($att=~/^Gap\=(.*)/);
		$Derives_from=$1 if ($att=~/^Derives_from\=(.*)/);
		$Note=$1 if ($att=~/^Note\=(.*)/);
		$Dbxref=$1 if ($att=~/Dbxref\=(.*)/);
		$Onto=$1 if ($att=~/^Ontology.*\=(.*)/);
		$repeat_type=$1 if ($att=~/^rpt_type\=(.*)/);
		$repeat_fam=$1 if ($att=~/^rpt_family\=(.*)/);
		$repeat_unit=$1 if ($att=~/^rpt_unit\=(.*)/);
		$product=$1 if ($att=~/^product\=(.*)/);
    }

	if($locusID){
		$gene_prod{$locusID}=$product unless ($gene_prod{$locusID});
	}
    elsif(! $locusID){
		if ($Parent){
			$locusID=$Parent."__exon"
		}
		elsif($type=~/repeat/){
			# rpt_type=CRISPR;rpt_unit=13023..13055;rpt_family=blah
			$locusID=$ID."__".
			($repeat_type ? "Type_".$repeat_type : "Type_Unknown")."__".
			($repeat_unit ? "Unit_".$repeat_unit : "Unit_Unknown")."__".
			($repeat_fam ? "Family_".$repeat_fam : "Family_Unknown")."__".
			"RepeatNum_".$repeatNumber;
			$repeatNumber++;
        }
        else{
			$locusID=$ID."__".$type;
        }
    }
	$annotation{$contig}{$locusID}{"START"}=$start;
	$annotation{$contig}{$locusID}{"STOP"}=$stop;
	$annotation{$contig}{$locusID}{"TYPE"}=$type;
	$annotation{$contig}{$locusID}{"LEN"}=($stop-$start);
	$annotation{$contig}{$locusID}{"STRAND"}=$strand;
}
