#! usr/bin/perl

=head1 USAGE

	USAGE: perl getSciNames.pl -i <tab-delimited-file> -db [n/p] <nucleotide or protein> -t <tax ids output; optional> -b <if input file is a blast output; mandatory>

=head2 LIMITATION:

	If you tell the script that the input is a blast output; it will only look at the subject cloumn;
	WARNING: will only work if the database followed "gi|gi_number|..." format (NCBI's default fasta format).

=cut

use strict;
use Getopt::Long; 


my ($list, $isBlastOut, $mapped, $db);
my $taxidOut=$$.".taxa";
my $out=$$.".names";
my $names_dmp="/geomicro/data1/COMMON/scripts/NCBI/taxa_dump/names.dmp";
my $dumpFiles="/geomicro/data1/COMMON/scripts/NCBI/taxa_dump/";
my $nucl_dmp="/geomicro/data1/COMMON/scripts/NCBI/taxa_dump/gi_taxid_nucl.dmp";
my $prot_dmp="/geomicro/data1/COMMON/scripts/NCBI/taxa_dump/gi_taxid_prot.dmp";

GetOptions(
	'l:s'=>\$list,
	'b|blast:s'=>\$isBlastOut,
	'm|map:s'=>\$mapped,
	'db:s'=>\$db,
	'taxa:s'=>\$taxidOut,
	'o|out:s'=>\$out,
	'names_dmp:s'=>\$names_dmp,
	'nucl_dmp:s'=>\$nucl_dmp,
	'prot_dmp:s'=>\$prot_dmp,
	'h|help'=> sub{system('perldoc', $0); exit;}
);

my $in;
if ($list){$in=$list}
elsif($isBlastOut){$in=$isBlastOut}
elsif($mapped){$in=$mapped}

if (! $in || ! $db) {
	print "Tab-Delimited File: $in\n";
	print "Database Used: $db\n";
	sub{system('perldoc', $0); exit;}
}

#ARGV[0]: blast output
open (LIST, $in) || die "[err] $in:\n$!\n"; 
my %idList;
while (my $line=<LIST>){
	next if $line=~ m/^#/;
	chomp $line;
	$line=~ s/\r//;
	next unless $line;

	if ($isBlastOut){
		my(@stuff)=split(/\t/, $line);
		my ($x,$sgi, @blah)=split(/\|/, $stuff[1]);
		$idList{$sgi}++;
	}
	elsif($list){
		$idList{$line}++;
	}
	elsif($mapped){
		my(@stuff)=split(/\t/, $line);
		my ($x,$sgi, @blah)=split(/\|/, $stuff[0]);
		$idList{$sgi}++;
	}
}
close LIST;

my ($fType, $dbType, $ref, $giList);
if (lc($db) eq "p"){
	$ref=$prot_dmp;
	$giList="giListP.tmp";
}
elsif (lc($db) eq "n"){
	$ref=$nucl_dmp;
	$giList="giListN.tmp";
}
else{
	warn "[err]: Invalid Database Type. Enter either \"p \(for Proteins\)\" OR \"n\(for Nucleotides\)\".\nYou Entered ".$db."\n";
	print "The Script will now exit.\n";
	exit;
}

open (INDEX, $ref) || die "[err] $ref:\n$!\n";
open (TAXID, ">".$taxidOut) if $taxidOut;
my (%index, %checkList);
while(my $line=<INDEX>){
	my($gi, $taxa)=split(/\t/, $line);
	chomp($gi);
	chomp($taxa);
	if ($idList{$gi}){
		$index{$gi}=$taxa;
		print TAXID $gi."\t".$taxa."\n" if $taxidOut;
		$checkList{$taxa}++;
	}
}
close INDEX;

my %taxa;
open (NAMES, $names_dmp) || die $!;
while (my $line2=<NAMES>){
	my(@thisLine)=split(/\t\|\t/, $line2); # tax_id, name, uniqName, commonName
	chomp(@thisLine);
	my $all=join("\t", @thisLine);
	$taxa{$thisLine[0]}=$all if ($checkList{$thisLine[0]});
}
undef %checkList;

#open (TD, ">sgi.log");
while (my($k, $v)=each(%index)){
	#print TD "$k\n";
	my @taxaD=split(/\t/, $taxa{$v});
	my $sciName=$taxaD[1];
	$index{$k}=$sciName
}
#close TD;
undef %taxa;

open (BOUT, $in) || die "[err] $in:\n$!\n"; 
open (OUT, ">".$out);
while (my $line=<BOUT>){
	unless($line=~ m/^#/){
		chomp $line;
		$line=~ s/\r//;
		next unless $line;
	
		if ($isBlastOut){
			my(@stuff3)=split(/\t/, $line);
			my ($x3,$sgi3, @blah3)=split(/\|/, $stuff3[1]);
#			$stuff3[1]=$index{$sgi3};
#			my $putBack=join ("\t", @stuff3);
			my $putBack=join ("\t", $index{$sgi3}, $line);

			print OUT $putBack."\n";
		}
		elsif($list){
			print OUT "$line\t$index{$line}\n";		
		}
		elsif($mapped){
			my(@stuff3)=split(/\t/, $line);
			my ($x, $sgi3, @blah3)=split(/\|/, $stuff3[0]);
#			$stuff3[0]=$index{$sgi3};
#			my $putBack=join ("\t", @stuff3);
			my $putBack=join ("\t", $stuff3[0],$index{$sgi3}, @stuff3[1..$#stuff3]);
			print OUT $putBack."\n";
		}
	}
	elsif($.==6 && $mapped){
		my @headers=split(/\t/, $line);
		my $newLine=join("\t", $headers[0],"Description", @headers[1..$#headers]);
		print OUT $newLine;
	}
	else{
		print OUT $line;
	}
}
close BOUT;
close OUT;
