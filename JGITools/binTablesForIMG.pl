#!/usr/bin/perl

=head1 DESCRIPTION

binTablesForIMG.pl -- Given a list of Contigs reformat it to upload to IMG Scaffold workspace.

=head1 USAGE

perl binTablesForIMG.pl

=head2 Options

        -map    -m   <CHAR>  names.map.txt
        -bins   -b   <CHAR>  Scaffolds_By_Sample_by_Bin/LOYAL_60 directory
        -outdir -o  <CHAR>  Bins_For_IMG/LOYAL_60
        -map_bins   <BOOLEAN>   Add bin numbers to project data downloaded from IMG in the map file.
        -version    -v  <BOOLEAN>   version of the current script
        -help   -h  <BOOLEAN>   This message.

=head1 Author

Sunit Jain, (Wed Feb 18 15:26:27 EST 2015)
sunitj [AT] umich [DOT] edu

=head1 License

This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;
use FileHandle;
use File::Spec qw(catdir catfile);
use File::Basename qw( fileparse );
use File::Path qw( make_path );

my $help;
my $version="fileparse($0)\tv0.1.0";
my $projectMap="names_map.txt";
my $outDir="Bins_For_IMG/LOYAL_60";
my $binsDir="Scaffolds_By_Sample_by_Bin/LOYAL_60";
my $mapBinData;

GetOptions(
        'm|map:s'=>\$projectMap,
        'b|bins:s'=>\$binsDir,
        'o|outdir:s'=>\$outDir,
        'map_bins'=>\$mapBinData,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

my %projects;
my $MAP=FileHandle->new();
open( $MAP, "<", $projectMap) || die $!;
while(my $line=<$MAP>){
    chomp $line;
    next unless ($line=~/^\d+/);
    next unless $line;
    my($imgName, $sampleNum, $desc)=split(/\t/, $line);
    $projects{$imgName}=$sampleNum;
    print "$imgName\t$sampleNum\t$desc\n";
}
close $MAP;


#foreach $project
#    read *.a.map.txt for originalScaffold <=> imgScaffold associations
#    look in corresponding $sampleNum $binDir for list files
#    foreach list file in $sampleNum $binDir
#        add column 1 = img Project name
#        replace originalScaffold with "assembled imgScaffold"
#        5 * \t
#        coverage
#   write new *.a.map.txt with column 3 = bin number

foreach my $proj(keys %projects){
    
    # Make an output directory if it doesn't already exist.
    my $projDirs=File::Spec->catdir($outDir,"IMG_".$proj."__".$projects{$proj});
    if (! -d $projDirs) {
        make_path $projDirs || die "Failed to create path: $projDirs\n";
    }

    # Read Map file into memory to relate original and IMG scaffold names.
    my $scafMapFile=File::Spec->catfile($proj, $proj.".a.map.txt");
    my (%scafMap, %bins);
    my $SCAF=FileHandle->new();
    open( $SCAF, "<", $scafMapFile) || die $!;
    while(my $line=<$SCAF>){
        next if ($line=~/^#/);
        chomp $line;
        next unless $line;
        my($origName, $imgName)=split(/\t/, $line);
        $scafMap{$origName}=$imgName;
    }
    close $SCAF;
    
    # Read IMG Scaffold file into Memory to seprate into Bins
    my %binSeq;
    my $projScafFasta=File::Spec->catfile($proj, $proj.".a.fna");
    my $FASTA=FileHandle->new();
    open($FASTA, "<",$projScafFasta) || die $!;
    $/=">";
    while (my $line=<$FASTA>) {
	chomp $line;
	next unless $line;
	my($header,@seqs)=split(/\n/,$line);
	$binSeq{$header}=join("",@seqs);
    }
    $/="\n";
    close $FASTA;
    
    if ($mapBinData){
	# Create a modified MAP file to use for consolidation by Bin.
	my $newScafMapFile=File::Spec->catfile($proj, $proj.".a.map.orig.txt");
	unless(-e $newScafMapFile){
	    rename($scafMapFile,$newScafMapFile);    
	}
	else{
	    print "[WARNING] File: $newScafMapFile already exists, NOT overwritten.\n";
	}
    }
    # Get Bins and their average coverage.
    my @binList=glob "$binsDir/$projects{$proj}/*.bin.list";
    next unless (scalar(@binList) > 0);
    # Create an IMG uploadable file for each bin and its corresponding fasta file.
    foreach my $binFile(@binList){
        my $binFileName=fileparse($binFile, ".bin.list");
        my $imgFile=File::Spec->catfile($projDirs, $binFileName.".tsv");
        my $binFasta=File::Spec->catfile($projDirs, $binFileName.".fna");
	
        my $BIN=FileHandle->new();
        my $IMG=FileHandle->new();
	my $BFASTA=FileHandle->new();
	
	open( $BFASTA, ">", $binFasta) || die $!;
        open( $IMG, ">", $imgFile) || die $!;
        print $IMG "Scaffold ID	Scaffold Name	Genome	Gene Count	Sequence Length	GC Content	Read Depth	Lineage Domain	Lineage Phylum	Lineage Class	Lineage Order	Lineage Family	Lineage Genus	Lineage Species	Lineage Percentage\n";
        open( $BIN, "<", $binFile) || die $!;

        while(my $line=<$BIN>){
            chomp $line;
            next if ($line=~/^#/);
            next unless $line;
            my($origName, $coverage)=split(/\s+/, $line);
	    my $imgContigName=$scafMap{$origName};
            print $IMG $proj." assembled ".$imgContigName."\t\t\t\t\t\t".$coverage."\n";
	    print $BFASTA ">".$imgContigName."\n".$binSeq{$imgContigName}."\n";
            if($mapBinData){
                if($bins{$origName}){
                    $bins{$origName}.="\t".$binFileName; 
                }
                else{
                    $bins{$origName}=$binFileName;
                }
            }
        }
        close $SCAF;
        close $IMG;
    }
    if ($mapBinData) {
        my $binScafMapFile=$scafMapFile;  # File::Spec->catfile($proj, $proj.".a.bins.map.txt");
        my $BINSCAF=FileHandle->new();
        open( $BINSCAF, ">", $binScafMapFile) || die $!;
        foreach my $scafName(keys %scafMap){
            print $BINSCAF $scafName."\t".$scafMap{$scafName};
            print $BINSCAF ($bins{$scafName} ? "\t".$bins{$scafName} : "");
            print $BINSCAF "\n";
        }
        close $BINSCAF;
    }
}

sub readFileWith2columns{
    my $file=shift;
    my %fileHash;
    
    return %fileHash;
}
