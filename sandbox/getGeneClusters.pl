#!/usr/local/bin/perl

=head1 DESCRIPTION

getGeneClusters.pl -- Do this.

=head1 USAGE

perl getGeneClusters.pl

=head2 Options


	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

Sunit Jain, (Fri Mar 20 11:04:27 EDT 2015)
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

use Bio::SeqIO;

my $help;
my $version=fileparse($0)."\tv0.0.4";
my ($gbkFile, $geneClustFile, $prefix);
GetOptions(
        'gbk:s'=>\$gbkFile,
        'gc|gene_clust:s'=>\$geneClustFile,
        'prefix:s'=>\$prefix,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

my %gbk;
parseGBK($gbkFile);

my $LIST=FileHandle->new();
my %index;
open( $LIST, "<", $geneClustFile) || die $!;
my $clustNum=0;
while(my $line=<$LIST>){
    chomp $line;
    next unless $line;
    next if ($line=~/^#/);
    my($shortName, $fullName, $clustProduct, $clustLocii, $acc)=split(/\t/, $line);
    
    $clustNum++;
    my $product=trim($clustProduct);
    my $fileName=$prefix."_product_".
                ($product ? $product : "Unknown").
                "_cnum_".$clustNum;
    my $faaFile=$fileName.".faa";
    my $metaFile=$fileName.".meta.txt";
    
    my $FAA=FileHandle->new();
    open($FAA, ">", $faaFile) || die $!;
    my $META=FileHandle->new();
    open($META, ">", $metaFile) || die $!;
    print $META "#Seq_Name\tCluster_Product\tORF_Start\tORF_Stop\t".
                "Locus_Tag\tMotif_Start\tMotif_Stop\tProduct\t".
                "PFam_Num\tPFamA_Hit\tScore\tEvalue\tDomain_Range\tNotes\n";
    foreach (split(/\;/, $clustLocii)){
        print $FAA ">".$_."\n";
        print $FAA $gbk{uc($_)}{"SEQ"}."\n";
        
        print $META $fullName."\t".$clustProduct."\t".
                    $gbk{uc($_)}{"GENE_START"}."\t".
                    $gbk{uc($_)}{"GENE_END"}."\t".
                    $_."\t".
                    $gbk{uc($_)}{"MOTIF_START"}."\t".
                    $gbk{uc($_)}{"MOTIF_END"}."\t".
                    $gbk{uc($_)}{"PRODUCT"}."\t".
                    $gbk{uc($_)}{"PFAM"}."\t".
                    $gbk{uc($_)}{"PFAMA"}."\t".
                    $gbk{uc($_)}{"SCORE"}."\t".
                    $gbk{uc($_)}{"EVALUE"}."\t".
                    $gbk{uc($_)}{"DR"}."\t".
                    $gbk{uc($_)}{"NOTES"}."\n";
    }
    close $FAA;
    close $META;
}
close $LIST;

##############################
######   Sub-Routines   ######  
##############################

sub trim{
    my $line=shift;
    chomp $line;
    $line=~s/^\s+//;
    $line=~s/\s+$//;
    $line=~s/\W+/_/g;
    return $line;
}

sub parseGBK{
    my $genbank = Bio::SeqIO->new(
                                -format=> 'genbank',
                                -file=> shift,
                                );
    
    while (my $seqObj=$genbank->next_seq) {
        foreach my $featObj ($seqObj->get_SeqFeatures){
            if($featObj->primary_tag eq "CDS"){
                my @headers=$featObj->get_tag_values('locus_tag');
                
                my @sequences=$featObj->get_tag_values('translation');
                my $start=$featObj->location->start;
                my $end=$featObj->location->end;
                
                $gbk{uc($headers[0])}{"GENE_START"}=$start;
                $gbk{uc($headers[0])}{"GENE_END"}=$end;
                $gbk{uc($headers[0])}{"SEQ"}=join("",@sequences);
            }
        
            if($featObj->primary_tag eq "CDS_motif"){
                my @printMeta;
                my @headers=eval{$featObj->get_tag_values('locus_tag')};
                next if $@;
                
                push(@printMeta, $headers[0]);
                
                my $start=$featObj->location->start;
                my $end=$featObj->location->end;
                
                my @notes=$featObj->get_tag_values('note');
                my $product=shift(@notes);
                
                $gbk{uc($headers[0])}{"MOTIF_START"}=$start;
                $gbk{uc($headers[0])}{"MOTIF_END"}=$end;
                $gbk{uc($headers[0])}{"PRODUCT"}=$product;
                
# Pfam-A.hmm-Hit: Hydrolase_like. Score: 44.5. E-value: 9.1e-12. Domain range: 2..73.;
                foreach (@notes){
                    if ($_=~/PFAM-Id: (.*)/) {
                        $gbk{uc($headers[0])}{"PFAM"}=$1;
                    }
                    elsif($_=~/Pfam-A.hmm-Hit: (.*). Score: (.*). E-value: (.*). Domain range: (.*)./){
                        $gbk{uc($headers[0])}{"PFAMA"}=$1;
                        $gbk{uc($headers[0])}{"SCORE"}=$2;
                        $gbk{uc($headers[0])}{"EVALUE"}=$3;
                        $gbk{uc($headers[0])}{"DR"}=$4;
                        $gbk{uc($headers[0])}{"NOTES"}.=$_.";";
                    }
                    else{
                        $gbk{uc($headers[0])}{"NOTES"}.=$_.";";
                    }
                }
            }
        }
    }
    return;
}
