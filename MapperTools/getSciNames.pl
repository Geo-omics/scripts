#! usr/bin/perl

=head1 USAGE

	perl getSciNames.pl -list <gi-list-file> -dbtype <nucl/prot>
	OR
	perl getSciNames.pl -blast <tabular-blast-output> -dbtype <nucl/prot>
	OR
	perl getSciNames.pl -map <tabular-mapper-output> -dbtype <nucl/prot>

=head2 Options

=head3 Input_type_(choose_one_only)

	-l	-list	:	if your input is a list of GI numbers
	-b	-blast	:	if your input is a blast output
	-m	-map	:	if your input file is a mapper output

=head3 Other_Options

	-dbtype		:	database type <"prot" or "nucl">; Required
	-o	-out	:	output file name <default=process_id.names>

	-taxa		:	prints Taxa ids to this file <optional>
	-dump		:	location of your *.dmp file downloaded from NCBI's Taxonomy FTP site; <default="/geomicro/data1/COMMON/scripts/Ungit/NCBI/taxa_dump/">
	-h	-help	:	This message; press "q" to exit this screen.
	
=head2 Limitation:

	If the input is a blast output, the script will only look at the subject cloumn and expect it to be in the NCBI FASTA header format (gi|gi_number|...)
	
=head1 Author

	Sunit Jain, May 2010
	sunitj-at-umich-dot-edu
	sunitjain-dot-com

=cut

use strict;
use Getopt::Long; 

my $version="1.1.0";
my ($list, $isBlastOut, $mapped, $db);
my $taxidOut;
my $out=$$.".names";
my $dumpFiles="/geomicro/data1/COMMON/scripts/Ungit/NCBI/taxa_dump/";

GetOptions(
	'l|list:s'=>\$list,
	'b|blast:s'=>\$isBlastOut,
	'm|map:s'=>\$mapped,
	'dbtype:s'=>\$db,
	'taxa:s'=>\$taxidOut,
	'o|out:s'=>\$out,
	'dump:s'=>\$dumpFiles,
	'h|help'=> sub{system('perldoc', $0); exit;},
	'v|version'=>sub{print "Version: $0\t v$version\n"; exit;},
);

$|++;

my $names_dmp=$dumpFiles."names.dmp";
my $nucl_dmp=$dumpFiles."gi_taxid_nucl.dmp";
my $nucl_diff_dmp=$dumpFiles."gi_taxid_nucl_diff.dmp";
my $prot_dmp=$dumpFiles."gi_taxid_prot.dmp";
my $prot_diff_dmp=$dumpFiles."gi_taxid_prot_diff.dmp";

my $in;
if ($list){$in=$list}
elsif($isBlastOut){$in=$isBlastOut}
elsif($mapped){$in=$mapped}

if (! $in || ! $db) {
	print "# Tab-Delimited File: $in\n";
	print "# Database Used: $db\n";
	sub{system('perldoc', $0); exit;}
}

print "# Parsing input...";
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
print "\tDONE\n";

my ($fType, $dbType, $ref, $giList, $diff);
if (lc($db) eq "prot"){
	$ref=$prot_dmp;
	$diff=$prot_diff_dmp;
	$giList="giListP.tmp";
}
elsif (lc($db) eq "nucl"){
	$ref=$nucl_dmp;
	$diff=$nucl_diff_dmp;	
	$giList="giListN.tmp";
}
else{
	warn "[err]: Invalid Database Type. Enter either \"p \(for Proteins\)\" OR \"n\(for Nucleotides\)\".\nYou Entered ".$db."\n";
	print "The Script will now exit.\n";
	exit;
}

print "# Checking Dump...";
open (INDEX, $ref) || die "[err] $ref:\n$!\n";
open (TAXID, ">".$taxidOut) if $taxidOut;
my (%index, %checkList);
while(my $line=<INDEX>){
	chomp $line;
	my($gi, $taxa)=split(/\t/, $line);
	if ($idList{$gi} && $taxa){
		$index{$gi}=$taxa;
		print TAXID $gi."\t".$taxa."\n" if $taxidOut;
		$checkList{$taxa}++;
	}
}
close INDEX;
print "\tDONE\n";

if (scalar(keys(%index)) < scalar(keys(%idList))){
	print "#\tOnly [ ".scalar(keys(%index))." ] GI numbers were current...\n";
	print "# Trying again...";
	open (DIFF, $diff) || die "[err] $diff:\n$!\n";
	while(my $line=<DIFF>){
		chomp $line;
		my($gi, $taxa)=split(/\t/, $line);
		if ($idList{$gi} && ! $index{$gi} && $taxa){
			$index{$gi}=$taxa;
			print TAXID $gi."\t".$taxa."\n" if $taxidOut;
			$checkList{$taxa}++;
		}
	}
	close DIFF;
	close TAXA;
	print "\tDONE\n";
	print "#\t[ ".scalar(keys(%index))." ] out of a possible [ ".scalar(keys(%idList))." ] GI numbers found...\n";
	print "#\tGiving Up...\n" if (scalar(keys(%index)) != scalar(keys(%idList)));
}

print "# Getting Taxonomy Names...";
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
print "\tDONE\n";

print "# Writing to the output file...";
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
			my $putBack=join ("\t", $index{$sgi3}, $line);

			print OUT $putBack."\n";
		}
		elsif($list){
			print OUT "$line\t$index{$line}\n";		
		}
		elsif($mapped){
			my(@stuff3)=split(/\t/, $line);
			my ($x, $sgi3, @blah3)=split(/\|/, $stuff3[0]);
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
print "\tDONE\n";
