#!/usr/bin/perl

=head1 DESCRIPTION

extractGenbankMetaData -- Extract Meta data and protein translations from a Genbank file. This script was created to be used as a part of the VirDB pipeline.
    
=head2 DEPENDENCIES

    BioPerl module Bio::SeqIO is required to run this script.

=head1 USAGE

    perl extractGenbankMetadata.pl -gbk in.gbk -faa out.faa -meta out.meta -list geneTags.list

=head2 Options

    -gbk    -g  <CHAR>  Genbank file (input)
    -list -l <CHAR>  List of Gene locus tags that you're interested in. (optional input; default: get everything)
    -meta   -m  <CHAR>  Metadata file (output 1)  
    -faa    -f  <CHAR>  Translated CDS (output 2)
    -null       <CHAR>  Set null value as (default: NA)
    -version -v     <BOOLEAN>       version of the current script
    -help   -h      <BOOLEAN>       This message.

=head1 Author

Sunit Jain, (Fri Jan 3 12:27:13 EST 2015)
sunitj [AT] umich [DOT] edu

=head1 License

This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

# Core Modules
use strict;
#use warnings;
use Getopt::Long;
use File::Basename;
use FileHandle;

# External Modules
use Bio::SeqIO;

my $self=fileparse($0);
my ($gbkFile,$metaFile,$locusFile,$faaFile);
my $nullValue="NA";
my $help;
my $version="$self\tv0.1.5";
GetOptions(
    'g|gbk:s'=>\$gbkFile,
    'f|faa:s'=>\$faaFile,
    'm|meta:s'=>\$metaFile,
    'l|list:s'=>\$locusFile,
    'null:s'=>\$nullValue,
    'v|version'=>sub{print $version."\n"; exit;},
    'h|help'=>sub{system("perldoc $0 | cat"); exit;},
);
print "\# $version\n";

my $LOCI=FileHandle->new();
my %locusTags;
if ($locusFile) {   
    open( $LOCI, "<", $locusFile) || die $!;
    while (my $line=<$LOCI>) {
        chomp $line;
        next if ($line=~ /^#/);
        next unless $line;
        
        $locusTags{uc($line)}++;
    }
    close $LOCI;
}

my $FAA=FileHandle->new();
open($FAA, ">", $faaFile) || die $!;
my $META=FileHandle->new();
open($META, ">", $metaFile) || die $!;
# Meta-Data Headers                
print $META "# BioProject\tContig\tContigLength\tisPlasmid\tContigDesc\tDBsources\t";
print $META "DB-xRef\tGeneLocus\tGeneStart\tGeneStop\tStrand\tProtein_Length\t";
print $META "Species\tTaxonomy\tProteinID\tProduct\tTranslation_table\tNotes\n";
parseGBK($gbkFile);
close $FAA;
close $META;

##############################
######   Sub-Routines   ######  
##############################

sub parseGBK{
    my $genbank = Bio::SeqIO->new(
                                -format=> 'genbank',
                                -file=> shift,
                                );
    
    while (my $seqObj=$genbank->next_seq) {
        # Contig Level Details that can be overwritten at the CDS level
        my($definition, $contigIsPlasmid, $LocusID, $Length, $species, $taxa, $bioProject, $xLinks);
        
        # Get definition line and check if the word palsmid is mentioned
        $definition=$seqObj->desc;
#        print $definition."\n";

        if ($definition=~ /plasmid/igo) {
            $contigIsPlasmid++;
        }
        
        # LOCUS Identifier
        $LocusID=$seqObj->display_id;
        
        # Contig/Scaffold Length
        $Length=$seqObj->length;
        
        eval { 
            # Organism
            $species = $seqObj->species->node_name;
            my @classification = $seqObj->species->classification;
            $taxa=join(";", @classification);
            
            if ($taxa eq ".") {
                $taxa=$definition;
            }
        
            # DB Cross links
            my @dblinks=$seqObj->get_Annotations('dblink');
            foreach(@dblinks){
                $xLinks.=$_->display_text."," ;
                if($_->display_text=~ /Project:(.*)/){
                    $bioProject=$1;
                }
            } 
            $xLinks=~ s/\,$//;
        }; warn $@ if $@;
        # Features for all CDSs in file.
        foreach my $featObj ($seqObj->get_SeqFeatures){
            if($featObj->primary_tag eq "CDS"){
                # CDS Level details that should be specific to the current CDS. Can overwrite Contig level variables as well.
                my($locus,$product, $start, $stop, $strand, $transTable, $note, $protID,$aa_length, $xRef, $translation, $markPlasmid);
                my $isPlasmid="";
                $start=$featObj->location->start;
                $stop=$featObj->location->end;
                $strand=$featObj->location->strand;
                
                # Get Locus Tag
                if ($featObj->has_tag('locus_tag')) {
                    my @locii=$featObj->get_tag_values("locus_tag");
                    $locus=join(",", @locii);
                    if ($locusFile) {
                        next unless $locusTags{uc($locus)};
                    }
                }
                else{
                    $locus=$nullValue;
                }
                
                # Check if Plasmid
                if($contigIsPlasmid){
                    $isPlasmid="Yes: Contig definition";
                    $markPlasmid++;
                }
                
                if ($featObj->has_tag('plasmid')) {
                    my @plasmids=$featObj->get_tag_values("plasmid");
                    $isPlasmid.=join(",", @plasmids);
                    $markPlasmid++;
                }
                elsif(! $contigIsPlasmid){
                    $isPlasmid="No";
                }
                
                # Get product description, if present
                if ($featObj->has_tag('product')) {
                    my @products=$featObj->get_tag_values("product");
                    $product=join(",", @products);
                }
                else{
                    $product=$nullValue;
                }
                
                # Get Translation Table, if present
                if ($featObj->has_tag('transl_table')) {
                    my @transTables=$featObj->get_tag_values("transl_table");
                    $transTable=join(",", @transTables);
                }
                else{
                    $transTable=$nullValue;
                }
                
                # Get Protein IDs, if present
                if ($featObj->has_tag('protein_id')) {
                    my @protIDs=$featObj->get_tag_values("protein_id");
                    $protID=join(",", @protIDs);
                    
                }
                else{
                    $protID=$nullValue;
                }
                
                # Get notes, if present
                if ($featObj->has_tag('note')) {
                    my @notes=$featObj->get_tag_values("note");
                    $note=join(";", @notes);
                }
                else{
                    $note=$nullValue;
                }
                
                # Get translation (protein seq), if present
                if ($featObj->has_tag('translation')) {
                    my @translations=$featObj->get_tag_values("translation");
                    $translation=join("", @translations);
                    $aa_length=length($translation);
                }
                else{
                    $translation="";
                    $aa_length=0;
                }
                
                # Get DB XREFs, if present
                if ($featObj->has_tag('db_xref')) {
                    my @xRefs=$featObj->get_tag_values("db_xref");
                    $xRef=join(",", @xRefs);
                }
                else{
                    $xRef=$nullValue;
                }
                
# Meta-Data Headers                
# print $META "# BioProject\tContig\tContigLength\tisPlasmid\tContigDesc\tDBsources\t";
# print $META "DB-xRef\tGeneLocus\tGeneStart\tGeneStop\tStrand\tProtein_Length\t";
# print $META "Species\tTaxonomy\tProteinID\tProduct\tTranslation_table\tNotes\n";
                
                if ($aa_length > 0) {   
                    # Write Meta-Data
                    ## Contig/Scaffold
                    print $META ($bioProject ? $bioProject : $nullValue)."\t".
                                $LocusID."\t".$Length."\t".$isPlasmid."\t".$definition."\t".
                                ($xLinks ? $xLinks : $nullValue)."\t";
                    ## Gene/Locus
                    print $META $xRef."\t".$locus."\t".$start."\t".$stop."\t".$strand."\t".$aa_length."\t";
                    print $META ($species ? $species : $nullValue)."\t".
                                ($taxa ? $taxa : $nullValue)."\t".
                                $protID."\t".$product."\t".$transTable."\t".$note."\n";
                    
                    # Write to protein fasta file
                    print $FAA ">".$locus.($LocusID ? "__".$LocusID : "")."\t".$protID."\t".$product."\t".$xRef.($markPlasmid ? "\tPLASMID" : "")."\n";
                    print $FAA $translation."\n";
                }
            }
        }
    }
    return;
}