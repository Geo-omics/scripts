#!/usr/bin/perl

=head1 DESCRIPTION

	extractTranslationsFromGbk.pl -- Do this.

=head1 USAGE

	perl extractTranslationsFromGbk.pl -list prophage_tbl.txt.list -gbk corresponding_genbank.gbk

=head2 Options

    -list
    -gbk
    -meta
    -fasta
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

	Sunit Jain, (Tue Nov 26 16:49:24 EST 2013)
	sunitj [AT] umich [DOT] edu

=head1 License

	This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

	This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;

my ($gbk,$meta, $list, $fasta);
my $help;
my $version="extractTranslationsFromGbk.pl\tv2.0.5";
GetOptions(
	'list:s'=>\$list,
	'gbk:s'=>\$gbk,
	'o|out|meta:s'=>\$meta,
	'f|fasta:s'=>\$fasta,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

open(LIST, $list)||die $!;
my %index;
while(my $line=<LIST>){
	$line=strip($line);
	$index{uc($line)}++;
}
close LIST;

sub strip{
	my $line=shift;
	chomp $line;
	$line=~ s/^\s+//;
	$line=~ s/\s+$//;
	return $line;
}

open(GBK, "<".$gbk)||die $!;
open(FASTA, ">".$fasta)||die $!;
my $found=0;
my ($accNum, $taxa, $desc, $inFeatures, $isPlasmid, $genomeLength, $bioproject);
$isPlasmid=0;
my %DB;
while(my $line=<GBK>){
	chomp $line;
	
	if($line=~ /^LOCUS/){
	    my ($TAG, @items)=split(/\s+/, $line);
	    foreach my $i(@items){
#	        print $i."\n";
	        next unless ($i=~ /^\d{3,10}$/);
	        $genomeLength=$i;
	    }
	    next;
	}
	
	if($line=~ /BioProject: (.*)/){
	    $bioproject=$1;
	    next;
	}
	
	if ($line=~ /^DEFINITION/){
		if($line=~ /plasmid/i){
		    $isPlasmid++;
		}
		$inFeatures=0;
		next;
	}
	
	if ($line=~ /^ACCESSION/){
		(my $TAG, $accNum)=split(/\s+/, $line);
		$inFeatures=0;
		next;
	}
	
	if ($line=~ /^SOURCE/){
		while($line!~ /ORGANISM/){
			$line=<GBK>;
		}
		my $organism=$line;
		$taxa=<GBK>;
		$organism=strip($organism);
		$taxa=strip($taxa);
		(my $org, $desc)=split(/  /, $organism);
		next;
	}
	
	$inFeatures++ if($line=~ /^FEATURES/);
	next unless ($inFeatures > 0);
	$line=strip($line);
#	print $line."\n";
	if ($line=~ m/^CDS/){
		my($geneLocus, $location, $proteinID, $gene, $product, $transTable, $start, $stop, $strand);
		$found=0;
		$geneLocus="";
		$gene="";
		$proteinID="";
		$product="";
		$strand=0;
		(my $TAG, $location)=split(/\s+/, $line);
        if($location=~ /(\d+)\.\.(\d+)/){
#            print "$1\t$2\n";
            $start=$1;
            $stop=$2;
        }
        if($location=~ /complement/){
            $strand="-1";
        }
        else{
            $strand="+1";
        }
        
		$line=<GBK>;
#     CDS             complement(1627..2319)
		while($line!~ m/^\ {5}\w+/){
            $line=strip($line);
			if($line=~ m#/locus_tag="(.*)"#){
				if ($index{uc($1)}){
#					print $1."\n";
					$found++;
					$geneLocus=$1;
					$DB{$accNum}{$geneLocus}{"Start"}=$start;
					$DB{$accNum}{$geneLocus}{"Stop"}=$stop;
					$DB{$accNum}{$geneLocus}{"Strand"}=$strand;
					$DB{$accNum}{$geneLocus}{"Taxa"}=$taxa;
					$DB{$accNum}{$geneLocus}{"Org"}=$desc;
				}
				else{
					last;
				}
			}
			elsif($line=~ m#/plasmid=(.*)#){
			    $isPlasmid++;
			}
			elsif($line=~ m#/gene=\"(.*)\"#){
				$DB{$accNum}{$geneLocus}{"Gene"}=$1;
			}
			elsif($line=~ m#/protein_id=\"(.*)\"#){
				$DB{$accNum}{$geneLocus}{"protID"}=$1;			
			}
			elsif($line=~ m#/product=\"(.*)\"#){
				$DB{$accNum}{$geneLocus}{"product"}=$1;			
			}
			elsif($line=~ m#/transl_table=(.*)#){
				$DB{$accNum}{$geneLocus}{"transTable"}=$1;
			}
			elsif (($line=~ m#/translation=\"(.*)#) && ($found > 0)){
				print FASTA ">".$geneLocus."\n";
#				print $1;
				my $seq=$1;
				while($line!~ m/\"$/){
					$line=<GBK>;
					$line=strip($line);
					$seq.=$line;
				}
				$seq=~ s/\"$//g;
				print FASTA $seq."\n";
				$found=0;
				$geneLocus="";
			}
			
		    if($isPlasmid > 0){
                $DB{$accNum}{$geneLocus}{"isPlasmid"}="Plasmid";		
	    	}
			$line=<GBK>;
		}
	}
}
close GBK;
close FASTA;

my @colOrder=qw(Start Stop Strand Org Taxa protID Gene product transTable, isPlasmid);
open(META, ">".$meta)||die $!;
print META "# BioProject\tContig\tContigLength\tGeneLocus\tGeneStart\tGeneStop\tStrand\tOrganism\tTaxonomy\tProteinID\tGene_Symbol\tProduct\tTranslationTable\tisPlasmid\n";
foreach my $acc(keys %DB){
	foreach my $geneLocus(keys %{$DB{$acc}}){
		my $line=$bioproject."\t".$acc."\t".$genomeLength."\t".$geneLocus;
		next unless ($DB{$acc}{$geneLocus}{"Strand"});
		foreach my $col(@colOrder){
			if($DB{$acc}{$geneLocus}{$col}){
				$line.="\t".$DB{$acc}{$geneLocus}{$col};
			}
			else{
				$line.="\t0";
			}
		}
		print META $line."\n";
	}
}
close META

